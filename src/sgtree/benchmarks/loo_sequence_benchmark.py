"""Small held-out-source sequence benchmark for US-009.

The benchmark simulates aligned proteins, infers every marker and species tree
with VeryFastTree, and keeps production LOO pruning disabled. It is an
engineering and safety gate, not a biological performance estimate.
"""

from __future__ import annotations

import argparse
import csv
import functools
import hashlib
import json
import math
import random
import shutil
import time
import warnings
from collections import defaultdict
from itertools import combinations
from pathlib import Path

warnings.simplefilter("ignore", SyntaxWarning)

from ete3 import Tree

from sgtree.benchmarks.loo_tree_fixtures import (
    _balanced,
    _nearest_genome_set,
    _opaque,
    _pectinate,
    _score_cmtv_current,
    cmtv_rf_quality_weights,
    wilson_interval,
)
from sgtree.benchmarks.sequence_evidence import assign_contig_gene_split_votes
from sgtree.id_schema import build_sequence_id, parse_sequence_id
from sgtree.marker_selection import (
    _rf_distance_between,
    build_singleton_output_tree,
    choose_singleton_prune,
    select_singleton_proposals,
)
from sgtree.marker_selection.contig_evidence import contig_gene_vote_gate
from sgtree.marker_selection.loo_profile import REVIEW_ROBUST_Z, score_loo_profiles
from sgtree.phylogeny import run_fasttree


REGIMES = {
    "alpha_like": {
        "topology": "balanced",
        "sentinel": (8, 0),
    },
    "bacteroidota_like": {
        "topology": "pectinate",
        "sentinel": (9, 3),
    },
}
SOURCE_DISTANCE_STRATA = {
    ("alpha_like", 607): "near",
    ("alpha_like", 809): "far",
    ("bacteroidota_like", 607): "far",
    ("bacteroidota_like", 809): "near",
}
DEVELOPMENT_SEEDS = (101, 401, 503)
SEEDS = (607, 809)
EVENT_CLASSES = ("clean", "gene_rich_contaminant", "solo_marker_contaminant")
TAXA_COUNT = 12
MARKER_COUNT = 8
NON_MARKER_GENE_COUNT = 3
NATIVE_CONTIG_GENE_COUNT = 3
GREEDY_MARKER_COUNT = 16
GREEDY_PANEL_TAG = "greedy16"
SWEEP_DONOR_GENE_COUNTS = (0, 1, 2, 3, 4, 5, 7, 10)
# Review-warning margin floor, calibrated on the v2 development instrument:
# truth informative-vote margin means were 31.1 and 50.5 bits, the highest
# false mean was 25.5. Two truth events only, so this is calibration, not
# validation. Donor-origin genes show a sharp attachment/background affinity
# cliff; native genes sit on a smooth within-panel divergence gradient.
REVIEW_MIN_VOTE_MARGIN = 28.0
MARKER_LENGTH = 400
NON_MARKER_LENGTH = 220
MAX_MARKER_TREE_BUILDS = 96
EXPECTED_DISTINCT_MARKER_TREE_BUILDS = 40
EXPECTED_MARKER_TREE_INSTANCES = 96
MAX_THREADS = 4
MAX_RUNTIME_SECONDS = 60.0
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
TREE_CACHE_SCHEMA = 1
SCORER_SIMPLICITY_ORDER = ("loo", "cmtv_weighted", "recipient_consensus")


def _seed(*parts: object) -> int:
    payload = "\x1f".join(str(part) for part in parts).encode()
    return int.from_bytes(hashlib.blake2s(payload, digest_size=8).digest(), "big")


def _attach_sister(
    topology: str | tuple,
    target: str,
    sister: str,
) -> str | tuple:
    if isinstance(topology, str):
        return (target, sister) if topology == target else topology
    return tuple(_attach_sister(child, target, sister) for child in topology)


def _leaf_path(
    topology: str | tuple,
    target: str,
    path: tuple[int, ...] = (),
) -> tuple[int, ...] | None:
    if isinstance(topology, str):
        return path if topology == target else None
    for index, child in enumerate(topology):
        result = _leaf_path(child, target, path + (index,))
        if result is not None:
            return result
    return None


def _topology_edge_distance(topology: str | tuple, left: str, right: str) -> int:
    left_path = _leaf_path(topology, left)
    right_path = _leaf_path(topology, right)
    if left_path is None or right_path is None:
        raise ValueError("source-distance taxon is absent from the fixed topology")
    shared = 0
    for left_step, right_step in zip(left_path, right_path):
        if left_step != right_step:
            break
        shared += 1
    return (len(left_path) - shared) + (len(right_path) - shared)


def _source_anchor(
    topology: str | tuple,
    genomes: tuple[str, ...],
    *,
    recipient: str,
    regime: str,
    seed: int,
) -> tuple[str, str, int]:
    stratum = SOURCE_DISTANCE_STRATA[(regime, seed)]
    distances = {
        genome: _topology_edge_distance(topology, recipient, genome)
        for genome in genomes
        if genome != recipient
    }
    target_distance = (
        min(distances.values())
        if stratum == "near"
        else max(distances.values())
    )
    choices = sorted(
        genome
        for genome, distance in distances.items()
        if distance == target_distance
    )
    choice = choices[_seed("source_anchor", regime, seed) % len(choices)]
    return choice, stratum, target_distance


def _mutate(sequence: str, *, seed_parts: tuple[object, ...], count: int) -> str:
    rng = random.Random(_seed(*seed_parts))
    result = list(sequence)
    for position in rng.sample(range(len(result)), count):
        result[position] = rng.choice(AMINO_ACIDS.replace(result[position], ""))
    return "".join(result)


def _sequence_difference_count(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def _simulate_family(
    topology: str | tuple,
    *,
    regime: str,
    family_key: tuple[object, ...],
    length: int,
) -> dict[str, str]:
    rng = random.Random(_seed(*family_key, "root"))
    root = [rng.choice(AMINO_ACIDS) for _ in range(length)]
    sequences: dict[str, str] = {}

    def descend(node: str | tuple, parent: list[str], depth: int, path: str) -> None:
        if isinstance(node, str):
            sequences[node] = "".join(parent)
            return
        for child_index, child in enumerate(node):
            child_sequence = list(parent)
            if REGIMES[regime]["topology"] == "balanced":
                mutation_count = 10 if depth < 2 else 6
            else:
                mutation_count = max(3, 8 - (depth // 3))
            branch_rng = random.Random(
                _seed(*family_key, "branch", path, child_index)
            )
            for position in branch_rng.sample(range(length), mutation_count):
                child_sequence[position] = branch_rng.choice(
                    AMINO_ACIDS.replace(child_sequence[position], "")
                )
            descend(
                child,
                child_sequence,
                depth + 1,
                f"{path}.{child_index}",
            )

    descend(topology, root, 0, "root")
    return sequences


def _event_record(
    fixture_key: tuple[object, ...],
    *,
    recipient: str,
    marker: str,
    kind: str,
) -> tuple[str, str]:
    contig = _opaque("c", *fixture_key, kind, "contig")
    gene = _opaque("q", *fixture_key, kind, marker)
    return contig, build_sequence_id(recipient, contig, gene)


def build_panel(
    regime: str,
    seed: int,
    *,
    marker_count: int = MARKER_COUNT,
    tag: str | None = None,
) -> dict:
    """Build one predeclared panel without calling a scorer or tree builder."""
    if regime not in REGIMES or seed not in SEEDS:
        raise ValueError("unknown fixed sequence-panel coordinate")

    panel_id = (
        _opaque("p", regime, seed)
        if tag is None
        else _opaque("p", regime, seed, tag)
    )
    fixture_key = (regime, seed) if tag is None else (regime, seed, tag)
    genomes = tuple(_opaque("x", panel_id, index) for index in range(TAXA_COUNT))
    markers = tuple(_opaque("m", panel_id, index) for index in range(marker_count))
    source_id = _opaque("s", panel_id, "held_out_source")
    main_recipient = genomes[2]
    main_marker = markers[0]
    topology = (
        _balanced(list(genomes))
        if REGIMES[regime]["topology"] == "balanced"
        else _pectinate(list(genomes))
    )
    main_anchor, distance_stratum, source_tree_edge_distance = _source_anchor(
        topology,
        genomes,
        recipient=main_recipient,
        regime=regime,
        seed=seed,
    )
    source_topology = _attach_sister(topology, main_anchor, source_id)

    marker_sequences: dict[str, dict[str, str]] = {}
    marker_records: dict[str, dict[str, str]] = {}
    source_marker_sequence = None
    for marker_index, marker in enumerate(markers):
        family_sequences = _simulate_family(
            source_topology if marker == main_marker else topology,
            regime=regime,
            family_key=(panel_id, "marker", marker_index),
            length=MARKER_LENGTH,
        )
        marker_sequences[marker] = {
            genome: family_sequences[genome]
            for genome in genomes
        }
        if marker == main_marker:
            source_marker_sequence = family_sequences[source_id]
        marker_records[marker] = {
            genome: build_sequence_id(
                genome,
                _opaque("c", panel_id, marker, genome),
                _opaque("q", panel_id, marker, genome),
            )
            for genome in genomes
        }

    non_marker_sequences: dict[str, dict[str, str]] = {genome: {} for genome in genomes}
    non_marker_records: dict[str, dict[str, str]] = {genome: {} for genome in genomes}
    source_non_marker_sequences: dict[str, str] = {}
    for gene_index in range(NON_MARKER_GENE_COUNT):
        family = _opaque("n", panel_id, gene_index)
        family_sequences = _simulate_family(
            source_topology,
            regime=regime,
            family_key=(panel_id, "non_marker", gene_index),
            length=NON_MARKER_LENGTH,
        )
        for genome in genomes:
            record_id = build_sequence_id(
                genome,
                _opaque("c", panel_id, "native_evidence", genome),
                _opaque("q", panel_id, "native_evidence", gene_index, genome),
            )
            non_marker_records[genome][family] = record_id
            non_marker_sequences[genome][family] = family_sequences[genome]
        source_non_marker_sequences[family] = family_sequences[source_id]

    # Native gene neighborhoods: every native marker contig carries the gate
    # minimum of non-marker genes, so contig audits on native records must be
    # decided by vote direction rather than gene scarcity.
    native_contig_genes: dict[str, dict[str, dict[str, str]]] = {
        genome: {} for genome in genomes
    }
    sweep_donor_gene_sequences: list[str] = [
        source_non_marker_sequences[family]
        for family in sorted(source_non_marker_sequences)
    ]
    for marker_index, marker in enumerate(markers):
        for slot in range(NATIVE_CONTIG_GENE_COUNT):
            family_sequences = _simulate_family(
                source_topology,
                regime=regime,
                family_key=(panel_id, "native_contig", marker_index, slot),
                length=NON_MARKER_LENGTH,
            )
            for genome in genomes:
                contig = _opaque("c", panel_id, marker, genome)
                record_id = build_sequence_id(
                    genome,
                    contig,
                    _opaque(
                        "q", panel_id, "native_contig", marker_index, slot, genome
                    ),
                )
                native_contig_genes[genome].setdefault(contig, {})[record_id] = (
                    family_sequences[genome]
                )
            sweep_donor_gene_sequences.append(family_sequences[source_id])

    main_contig, main_record = _event_record(
        fixture_key,
        recipient=main_recipient,
        marker=main_marker,
        kind="source_replacement",
    )
    if source_marker_sequence is None:
        raise RuntimeError("held-out source marker was not simulated")
    main_sequence = source_marker_sequence
    main_genes = {
        build_sequence_id(
            main_recipient,
            main_contig,
            _opaque("q", *fixture_key, "source", gene_index),
        ): source_non_marker_sequences[family]
        for gene_index, family in enumerate(sorted(non_marker_sequences[main_anchor]))
    }
    source_anchor_gene_differences = [
        _sequence_difference_count(
            source_non_marker_sequences[family],
            non_marker_sequences[main_anchor][family],
        )
        for family in sorted(source_non_marker_sequences)
    ]

    sentinel_recipient_index, sentinel_anchor_index = REGIMES[regime]["sentinel"]
    sentinel_recipient = genomes[sentinel_recipient_index]
    sentinel_anchor = genomes[sentinel_anchor_index]
    sentinel_marker = markers[1]
    sentinel_contig, sentinel_record = _event_record(
        fixture_key,
        recipient=sentinel_recipient,
        marker=sentinel_marker,
        kind="native_contig_sentinel",
    )
    sentinel_sequence = _mutate(
        marker_sequences[sentinel_marker][sentinel_anchor],
        seed_parts=(*fixture_key, "sentinel", "marker"),
        count=2,
    )
    sentinel_genes = {
        build_sequence_id(
            sentinel_recipient,
            sentinel_contig,
            _opaque("q", *fixture_key, "sentinel", gene_index),
        ): non_marker_sequences[sentinel_recipient][family]
        for gene_index, family in enumerate(
            sorted(non_marker_sequences[sentinel_recipient])
        )
    }

    return {
        "panel_id": panel_id,
        "regime": regime,
        "seed": seed,
        "tag": tag,
        "topology_name": REGIMES[regime]["topology"],
        "genomes": genomes,
        "markers": markers,
        "source_id": source_id,
        "topology": topology,
        "marker_sequences": marker_sequences,
        "marker_records": marker_records,
        "non_marker_sequences": non_marker_sequences,
        "non_marker_records": non_marker_records,
        "native_contig_genes": native_contig_genes,
        "sweep_donor_gene_sequences": sweep_donor_gene_sequences,
        "main_event": {
            "event_id": _opaque("e", *fixture_key, "source_replacement"),
            "event_kind": "source_replacement",
            "is_contaminant": True,
            "recipient": main_recipient,
            "source_id": source_id,
            "source_generation": "joint_13_taxon_sister_simulation",
            "anchor": main_anchor,
            "distance_stratum": distance_stratum,
            "source_tree_edge_distance": source_tree_edge_distance,
            "source_anchor_rule": "truth_topology_extreme_then_blake2_tiebreak",
            "marker": main_marker,
            "native_record_id": marker_records[main_marker][main_recipient],
            "observed_record_id": main_record,
            "contig_id": main_contig,
            "sequence": main_sequence,
            "source_anchor_marker_difference_count": _sequence_difference_count(
                main_sequence,
                marker_sequences[main_marker][main_anchor],
            ),
            "source_anchor_gene_difference_counts": source_anchor_gene_differences,
            "candidate_genes": main_genes,
        },
        "sentinel_event": {
            "event_id": _opaque("e", *fixture_key, "native_contig_sentinel"),
            "event_kind": "native_contig_sentinel",
            "is_contaminant": False,
            "recipient": sentinel_recipient,
            "anchor": sentinel_anchor,
            "marker": sentinel_marker,
            "native_record_id": marker_records[sentinel_marker][sentinel_recipient],
            "observed_record_id": sentinel_record,
            "contig_id": sentinel_contig,
            "sequence": sentinel_sequence,
            "candidate_genes": sentinel_genes,
        },
    }


def build_panels() -> list[dict]:
    return [build_panel(regime, seed) for regime in REGIMES for seed in SEEDS]


def _sequence_hashes(records: dict[str, str]) -> dict[str, str]:
    return {
        record_id: hashlib.sha256(sequence.encode()).hexdigest()
        for record_id, sequence in sorted(records.items())
    }


def _panel_sequence_hash(panel: dict) -> str:
    """Hash every identifier and sequence that can affect a benchmark score."""
    payload = {
        "marker_records": panel["marker_records"],
        "marker_sequences": panel["marker_sequences"],
        "non_marker_records": panel["non_marker_records"],
        "non_marker_sequences": panel["non_marker_sequences"],
        "native_contig_genes": panel["native_contig_genes"],
        "sweep_donor_gene_sequences": panel["sweep_donor_gene_sequences"],
        "events": {
            key: {
                "observed_record_id": panel[key]["observed_record_id"],
                "sequence": panel[key]["sequence"],
                "candidate_genes": panel[key]["candidate_genes"],
            }
            for key in ("main_event", "sentinel_event")
        },
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()


def build_manifest(panels: list[dict] | None = None) -> dict:
    """Return the fixed design manifest written before any held-out score."""
    panels = panels or build_panels()
    panel_rows = []
    case_rows = []
    for panel in panels:
        events = []
        for event_key in ("main_event", "sentinel_event"):
            event = panel[event_key]
            event_row = {
                key: event[key]
                for key in (
                    "event_id",
                    "event_kind",
                    "is_contaminant",
                    "recipient",
                    "anchor",
                    "marker",
                    "native_record_id",
                    "observed_record_id",
                    "contig_id",
                )
            }
            event_row["candidate_gene_hashes"] = _sequence_hashes(
                event["candidate_genes"]
            )
            for key in (
                "source_id",
                "source_generation",
                "source_anchor_marker_difference_count",
                "source_anchor_gene_difference_counts",
                "source_tree_edge_distance",
                "source_anchor_rule",
            ):
                if key in event:
                    event_row[key] = event[key]
            events.append(event_row)
        panel_rows.append(
            {
                "panel_id": panel["panel_id"],
                "regime": panel["regime"],
                "seed": panel["seed"],
                "topology": panel["topology_name"],
                "genomes": list(panel["genomes"]),
                "markers": list(panel["markers"]),
                "source_id": panel["source_id"],
                "source_excluded_from_panel": panel["source_id"] not in panel["genomes"],
                "distance_stratum": panel["main_event"]["distance_stratum"],
                "sequence_design_sha256": _panel_sequence_hash(panel),
                "events": events,
            }
        )
        for event_class in EVENT_CLASSES:
            case_rows.append(
                {
                    "case_id": _case_opaque(panel, event_class),
                    "panel_id": panel["panel_id"],
                    "regime": panel["regime"],
                    "seed": panel["seed"],
                    "event_class": event_class,
                    "distance_stratum": panel["main_event"]["distance_stratum"],
                }
            )
    marker_counts = sorted({len(panel["markers"]) for panel in panels})
    return {
        "benchmark_kind": "synthetic_aligned_protein_inference",
        "primary_estimand": "truth_topology_far_source_single_copy_contamination",
        "near_source_role": "identifiability_control_expected_abstention",
        "development_seeds": list(DEVELOPMENT_SEEDS),
        "held_out_seeds": list(SEEDS),
        "production_pruning_enabled": False,
        "scorer_simplicity_order": list(SCORER_SIMPLICITY_ORDER),
        "taxa_per_panel": TAXA_COUNT,
        "markers_per_case": (
            marker_counts[0] if len(marker_counts) == 1 else marker_counts
        ),
        "panel_realism": "native_contig_genes_v2",
        "native_contig_gene_count": NATIVE_CONTIG_GENE_COUNT,
        "panels": panel_rows,
        "cases": case_rows,
    }


def _case_state(panel: dict, event_class: str) -> dict:
    records = {
        marker: dict(per_genome)
        for marker, per_genome in panel["marker_records"].items()
    }
    sequences = {
        marker: dict(per_genome)
        for marker, per_genome in panel["marker_sequences"].items()
    }
    contig_genes: dict[str, dict[str, str]] = {}
    for genome in panel["genomes"]:
        for contig, genes in panel["native_contig_genes"][genome].items():
            contig_genes[contig] = dict(genes)
    events = []
    if event_class in {"gene_rich_contaminant", "solo_marker_contaminant"}:
        event = panel["main_event"]
        records[event["marker"]][event["recipient"]] = event["observed_record_id"]
        sequences[event["marker"]][event["recipient"]] = event["sequence"]
        events.append(event)
        if event_class == "gene_rich_contaminant":
            contig_genes[event["contig_id"]] = dict(event["candidate_genes"])
            sentinel = panel["sentinel_event"]
            records[sentinel["marker"]][sentinel["recipient"]] = sentinel[
                "observed_record_id"
            ]
            sequences[sentinel["marker"]][sentinel["recipient"]] = sentinel[
                "sequence"
            ]
            contig_genes[sentinel["contig_id"]] = dict(
                sentinel["candidate_genes"]
            )
            events.append(sentinel)
    return {
        "records": records,
        "sequences": sequences,
        "contig_genes": contig_genes,
        "events": events,
    }


def _reference_proteomes(panel: dict) -> dict[str, dict[str, str]]:
    result = {genome: {} for genome in panel["genomes"]}
    for marker in panel["markers"]:
        for genome in panel["genomes"]:
            result[genome][panel["marker_records"][marker][genome]] = panel[
                "marker_sequences"
            ][marker][genome]
    for genome in panel["genomes"]:
        for family in sorted(panel["non_marker_records"][genome]):
            result[genome][panel["non_marker_records"][genome][family]] = panel[
                "non_marker_sequences"
            ][genome][family]
    for genome in panel["genomes"]:
        for contig in sorted(panel["native_contig_genes"][genome]):
            result[genome].update(panel["native_contig_genes"][genome][contig])
    return result


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{record_id}\n{sequence}\n" for record_id, sequence in records))


@functools.cache
def _tree_builder_signature() -> str:
    executable = next(
        (
            resolved
            for candidate in ("VeryFastTree", "veryfasttree")
            if (resolved := shutil.which(candidate)) is not None
        ),
        None,
    )
    if executable is None:
        raise FileNotFoundError("Could not find VeryFastTree in PATH")
    digest = hashlib.sha256()
    with Path(executable).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tree_cache_digest(records: list[tuple[str, str]], threads: int) -> str:
    payload = {
        "schema": TREE_CACHE_SCHEMA,
        "builder_sha256": _tree_builder_signature(),
        "threads": threads,
        "records": records,
    }
    serialized = json.dumps(payload, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()


def _tree_has_expected_leaves(tree_path: Path, expected: list[str]) -> bool:
    try:
        observed = [str(leaf.name) for leaf in Tree(str(tree_path), format=1).iter_leaves()]
    except Exception:
        return False
    return len(observed) == len(expected) and set(observed) == set(expected)


def _ensure_tree(
    alignment_path: Path,
    tree_path: Path,
    records: list[tuple[str, str]],
    *,
    threads: int,
) -> bool:
    cache_path = tree_path.with_suffix(tree_path.suffix + ".cache.json")
    expected_leaves = [record_id for record_id, _sequence in records]
    input_digest = _tree_cache_digest(records, threads)
    if tree_path.exists() and cache_path.exists() and tree_path.stat().st_size > 0:
        try:
            cached = json.loads(cache_path.read_text())
        except (OSError, json.JSONDecodeError):
            cached = {}
        if (
            cached.get("input_sha256") == input_digest
            and _tree_has_expected_leaves(tree_path, expected_leaves)
        ):
            return False
    _write_fasta(alignment_path, records)
    tree_path.parent.mkdir(parents=True, exist_ok=True)
    run_fasttree(str(alignment_path), str(tree_path), threads)
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        raise RuntimeError(f"VeryFastTree did not create {tree_path}")
    if not _tree_has_expected_leaves(tree_path, expected_leaves):
        raise RuntimeError(f"VeryFastTree tree has the wrong leaves: {tree_path}")
    cache_path.write_text(
        json.dumps(
            {
                "schema": TREE_CACHE_SCHEMA,
                "input_sha256": input_digest,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return True


def _variant_alignment(panel: dict, marker: str, variant: str) -> list[tuple[str, str]]:
    records = dict(panel["marker_records"][marker])
    sequences = dict(panel["marker_sequences"][marker])
    if variant == "source":
        event = panel["main_event"]
    elif variant == "sentinel":
        event = panel["sentinel_event"]
    elif variant == "clean":
        event = None
    else:
        raise ValueError(f"unknown marker-tree variant: {variant}")
    if event is not None:
        if event["marker"] != marker:
            raise ValueError("event variant requested for the wrong marker")
        records[event["recipient"]] = event["observed_record_id"]
        sequences[event["recipient"]] = event["sequence"]
    return [(records[genome], sequences[genome]) for genome in panel["genomes"]]


def _build_marker_trees(panel: dict, outdir: Path, threads: int) -> tuple[dict, int]:
    panel_dir = outdir / "marker_cache" / panel["panel_id"]
    paths: dict[tuple[str, str], Path] = {}
    built = 0
    for marker in panel["markers"]:
        variants = ["clean"]
        if marker == panel["main_event"]["marker"]:
            variants.append("source")
        if marker == panel["sentinel_event"]["marker"]:
            variants.append("sentinel")
        for variant in variants:
            stem = f"{variant}_{marker}"
            alignment_path = panel_dir / "alignments" / f"{stem}.faa"
            tree_path = panel_dir / "trees" / f"{stem}.nwk"
            built += int(
                _ensure_tree(
                    alignment_path,
                    tree_path,
                    _variant_alignment(panel, marker, variant),
                    threads=threads,
                )
            )
            paths[(marker, variant)] = tree_path
    return paths, built


def _case_tree_paths(panel: dict, event_class: str, cache: dict) -> dict[str, Path]:
    result = {}
    for marker in panel["markers"]:
        variant = "clean"
        if (
            event_class in {"gene_rich_contaminant", "solo_marker_contaminant"}
            and marker == panel["main_event"]["marker"]
        ):
            variant = "source"
        if (
            event_class == "gene_rich_contaminant"
            and marker == panel["sentinel_event"]["marker"]
        ):
            variant = "sentinel"
        result[marker] = cache[(marker, variant)]
    return result


def _species_alignment(
    panel: dict,
    state: dict,
    removals: set[tuple[str, str]],
) -> list[tuple[str, str]]:
    records = []
    for genome in panel["genomes"]:
        blocks = []
        for marker in panel["markers"]:
            record_id = state["records"][marker][genome]
            sequence = state["sequences"][marker][genome]
            blocks.append(
                "-" * len(sequence)
                if (marker, record_id) in removals
                else sequence
            )
        records.append((genome, "".join(blocks)))
    return records


def _species_tree(
    panel: dict,
    state: dict,
    *,
    case_id: str,
    arm: str,
    removals: set[tuple[str, str]],
    outdir: Path,
    threads: int,
) -> tuple[Path, bool]:
    digest = hashlib.blake2s(
        json.dumps(sorted(removals), separators=(",", ":")).encode(),
        digest_size=6,
    ).hexdigest()
    stem = f"{case_id}_{arm}_{digest}"
    alignment_path = outdir / "species" / "alignments" / f"{stem}.faa"
    tree_path = outdir / "species" / "trees" / f"{stem}.nwk"
    built = _ensure_tree(
        alignment_path,
        tree_path,
        _species_alignment(panel, state, removals),
        threads=threads,
    )
    return tree_path, built


def _rf_norm(reference_path: Path, observed_path: Path) -> float:
    reference = Tree(str(reference_path))
    observed = Tree(str(observed_path))
    shared = sorted(
        {str(leaf.name) for leaf in reference.iter_leaves()}
        & {str(leaf.name) for leaf in observed.iter_leaves()}
    )
    if len(shared) < 4:
        return 0.0
    reference.prune(shared)
    observed.prune(shared)
    rf, maxrf, *_ = reference.robinson_foulds(
        observed,
        unrooted_trees=True,
    )
    return rf / maxrf if maxrf else 0.0


def _patristic_nrmse(reference_path: Path, observed_path: Path) -> float:
    reference = Tree(str(reference_path))
    observed = Tree(str(observed_path))
    reference_leaves = {
        str(leaf.name): leaf
        for leaf in reference.iter_leaves()
    }
    observed_leaves = {
        str(leaf.name): leaf
        for leaf in observed.iter_leaves()
    }
    shared = sorted(reference_leaves.keys() & observed_leaves.keys())
    pairs = list(combinations(shared, 2))
    if not pairs:
        return 0.0
    reference_distances = [
        float(reference_leaves[left].get_distance(reference_leaves[right]))
        for left, right in pairs
    ]
    observed_distances = [
        float(observed_leaves[left].get_distance(observed_leaves[right]))
        for left, right in pairs
    ]
    scale = sum(reference_distances) / len(reference_distances)
    if scale <= 0.0:
        return 0.0
    rmse = math.sqrt(
        sum(
            (observed_distance - reference_distance) ** 2
            for reference_distance, observed_distance in zip(
                reference_distances,
                observed_distances,
                strict=True,
            )
        )
        / len(pairs)
    )
    return rmse / scale


def _marker_trees(paths: dict[str, Path]) -> dict[str, Tree]:
    return {
        marker: Tree(str(path), format=1)
        for marker, path in sorted(paths.items())
    }


def _all_marker_record_ids(panel: dict) -> set[str]:
    result = {
        record_id
        for records in panel["marker_records"].values()
        for record_id in records.values()
    }
    result.add(panel["main_event"]["observed_record_id"])
    result.add(panel["sentinel_event"]["observed_record_id"])
    return result


def _case_opaque(panel: dict, event_class: str) -> str:
    tag = panel.get("tag")
    if tag is None:
        return _opaque("f", panel["regime"], event_class, panel["seed"])
    return _opaque("f", panel["regime"], event_class, panel["seed"], tag)


def _select_with_budget(proposals: list[dict], panel: dict) -> list[dict]:
    return select_singleton_proposals(
        proposals,
        genome_marker_counts={
            genome: len(panel["markers"]) for genome in panel["genomes"]
        },
        min_markers_per_genome=1,
        max_prunes_per_genome=1,
    )


def _informative_margin_stats(
    votes: list[dict],
) -> tuple[float | None, float | None]:
    margins = [
        float(vote["score_margin"])
        for vote in votes
        if vote.get("informative") and vote.get("score_margin") is not None
    ]
    if not margins:
        return None, None
    return sum(margins) / len(margins), min(margins)


def _loo_proposals(rows: list[dict]) -> list[dict]:
    return [
        {
            "marker_name": str(row["marker_name"]),
            "leaf_name": str(row["leaf_name"]),
            "genome": str(row["genome"]),
            "score": float(row["loo_score"]),
            "delta_rf": 0.0,
            "topoknn_score": 0.0,
        }
        for row in rows
        if row["loo_class"] == "discordant_marker"
    ]


def _common_attachment_audits(
    panel: dict,
    state: dict,
    trees: dict[str, Tree],
    rows_by_key: dict[tuple[str, str], dict],
    candidate_keys: set[tuple[str, str]],
) -> dict[tuple[str, str], dict]:
    reference_proteomes = _reference_proteomes(panel)
    marker_record_ids = _all_marker_record_ids(panel)
    audits: dict[tuple[str, str], dict] = {}
    for key in sorted(candidate_keys):
        row = rows_by_key[key]
        genome, contig_id, _gene_id = parse_sequence_id(str(row["leaf_name"]))
        candidate_genes = state["contig_genes"].get(contig_id, {})
        attachment = set(row["loo_attachment_taxa"])
        background = {
            str(leaf.name).split("|", 1)[0]
            for leaf in trees[str(row["marker_name"])].iter_leaves()
        } - {genome} - attachment
        evidence = assign_contig_gene_split_votes(
            candidate_genes,
            reference_proteomes,
            recipient_genome=genome,
            candidate_contig_id=contig_id,
            marker_record_ids=marker_record_ids,
            attachment_taxa=attachment,
            background_taxa=background,
        )
        votes = evidence["votes"]
        gate = contig_gene_vote_gate(votes, row["loo_attachment_clade"])
        audits[key] = {"evidence": evidence, "gate": gate}
    return audits


def _guard_proposals(
    proposals: list[dict],
    *,
    marker_paths: dict[str, Path],
    species_tree_path: Path,
    mode: str,
) -> tuple[set[tuple[str, str]], list[dict]]:
    by_marker: dict[str, list[dict]] = defaultdict(list)
    for proposal in proposals:
        by_marker[str(proposal["marker_name"])].append(proposal)

    removals: set[tuple[str, str]] = set()
    decisions = []
    species_tree = Tree(str(species_tree_path))
    for marker, marker_proposals in sorted(by_marker.items()):
        marker_path = marker_paths[marker]
        original = Tree(str(marker_path), format=1)
        accepted_leaf_names = sorted(
            str(proposal["leaf_name"])
            for proposal in marker_proposals
        )
        accepted_leaf_set = set(accepted_leaf_names)
        remaining = [
            str(leaf.name)
            for leaf in original.iter_leaves()
            if str(leaf.name) not in accepted_leaf_set
        ]
        candidate = original.copy(method="deepcopy")
        candidate.prune(remaining)
        rf_before = _rf_distance_between(species_tree, original)
        rf_after = _rf_distance_between(species_tree, candidate)
        chosen, decision = build_singleton_output_tree(
            marker_tree_path=str(marker_path),
            species_tree_path=str(species_tree_path),
            accepted_leaf_names=accepted_leaf_names,
            mode=mode,
        )
        chosen_leaves = {str(leaf.name) for leaf in chosen.iter_leaves()}
        marker_removed = {
            (marker, leaf_name)
            for leaf_name in accepted_leaf_names
            if leaf_name not in chosen_leaves
        }
        removals.update(marker_removed)
        decisions.append(
            {
                "marker_name": marker,
                "proposed_leaf_names": accepted_leaf_names,
                "removed_leaf_names": sorted(
                    leaf_name for _marker, leaf_name in marker_removed
                ),
                "rf_before": rf_before,
                "rf_after": rf_after,
                "rf_guard_decision": decision,
            }
        )
    return removals, decisions


def _cmtv_proposals(cmtv_rows) -> list[dict]:
    return [
        {
            "marker_name": str(row.marker_name),
            "leaf_name": str(row.leaf_name),
            "genome": str(row.genome),
            "score": float(row.cmtv_combined),
            "delta_rf": 0.0,
            "topoknn_score": 0.0,
        }
        for row in cmtv_rows.itertuples(index=False)
        if row.cmtv_class == "contamination_candidate"
    ]


def _recipient_consensus_proposals(
    trees: dict[str, Tree],
    marker_paths: dict[str, Path],
    species_tree_path: Path,
) -> list[dict]:
    species_tree = Tree(str(species_tree_path))
    proposals = []
    for marker, tree in sorted(trees.items()):
        marker_path = marker_paths[marker]
        alignment_path = (
            marker_path.parents[1] / "alignments" / f"{marker_path.stem}.faa"
        )
        chosen = choose_singleton_prune(
            species_tree,
            tree,
            mode="recipient_consensus",
            k=min(10, TAXA_COUNT - 1),
            alignment_path=str(alignment_path),
        )
        if chosen is not None:
            proposals.append({**chosen, "marker_name": marker})
    return proposals


def _comparison_rows(
    *,
    case_id: str,
    panel: dict,
    event_class: str,
    events: list[dict],
    rows_by_key: dict[tuple[str, str], dict],
    audits: dict[tuple[str, str], dict],
    pipelines: dict[str, dict],
) -> list[dict]:
    events_by_key = {
        (event["marker"], event["observed_record_id"]): event
        for event in events
    }
    raw_keys = {
        scorer: {
            (str(proposal["marker_name"]), str(proposal["leaf_name"]))
            for proposal in pipeline["raw"]
        }
        for scorer, pipeline in pipelines.items()
    }
    budgeted_keys = {
        scorer: {
            (str(proposal["marker_name"]), str(proposal["leaf_name"]))
            for proposal in pipeline["budgeted"]
        }
        for scorer, pipeline in pipelines.items()
    }
    raw_scores = {
        scorer: {
            (str(proposal["marker_name"]), str(proposal["leaf_name"])): float(
                proposal["score"]
            )
            for proposal in pipeline["raw"]
        }
        for scorer, pipeline in pipelines.items()
    }
    guard_decisions = {}
    for scorer, pipeline in pipelines.items():
        guard_decisions[scorer] = {
            (str(decision["marker_name"]), str(leaf_name)): str(
                decision["rf_guard_decision"]
            )
            for decision in pipeline["guards"]
            for leaf_name in decision["proposed_leaf_names"]
        }

    keys = set(events_by_key)
    for keys_for_scorer in raw_keys.values():
        keys.update(keys_for_scorer)

    result = []
    for key in sorted(keys):
        marker, record_id = key
        row = rows_by_key[key]
        event = events_by_key.get(key)
        audit = audits.get(key, {})
        gate = audit.get("gate", {})
        evidence = audit.get("evidence", {})
        for scorer in pipelines:
            raw_call = key in raw_keys[scorer]
            contig_gate_pass = (
                bool(gate.get("contig_gate_pass"))
                if raw_call
                else None
            )
            budget_pass = key in budgeted_keys[scorer]
            final_removed = key in pipelines[scorer]["removals"]
            rf_decision = guard_decisions[scorer].get(key)
            if not raw_call:
                reason = "no_raw_call"
                gate_result = "not_evaluated"
            elif not contig_gate_pass:
                gate_reason = gate.get("contig_abstention_reason")
                reason = (
                    "attachment_abstention"
                    if gate_reason == "missing_proposed_clade"
                    else "evidence_abstention"
                )
                gate_result = "failed"
            elif not budget_pass:
                reason = "budget_rejected"
                gate_result = "passed"
            elif not final_removed:
                reason = "rf_guard_rejected"
                gate_result = "passed"
            else:
                reason = "removed"
                gate_result = "passed"
            result.append(
                {
                    "case_id": case_id,
                    "panel_id": panel["panel_id"],
                    "regime": panel["regime"],
                    "seed": panel["seed"],
                    "event_class": event_class,
                    "distance_stratum": panel["main_event"]["distance_stratum"],
                    "event_id": event["event_id"] if event else "",
                    "event_kind": event["event_kind"] if event else "unsupported_call",
                    "is_declared_event": event is not None,
                    "is_contaminant": bool(event and event["is_contaminant"]),
                    "marker_name": marker,
                    "record_id": record_id,
                    "recipient": (
                        event["recipient"]
                        if event
                        else parse_sequence_id(record_id)[0]
                    ),
                    "scorer": scorer,
                    "raw_call": raw_call,
                    "raw_score": raw_scores[scorer].get(key),
                    "attachment_clade": row["loo_attachment_clade"],
                    "evidence_status": (
                        evidence.get("input_status")
                        if raw_call
                        else None
                    ),
                    "contig_gate_pass": contig_gate_pass,
                    "gate_result": gate_result,
                    "gate_reason": (
                        gate.get("contig_abstention_reason")
                        if raw_call
                        else None
                    ),
                    "budget_pass": budget_pass,
                    "rf_guard_decision": rf_decision,
                    "final_action": "removed" if final_removed else "kept",
                    "reason": reason,
                }
            )
    return result


def _score_case(
    panel: dict,
    event_class: str,
    *,
    marker_cache: dict,
    outdir: Path,
    threads: int,
    clean_species_path: Path,
) -> tuple[dict, int]:
    case_id = _case_opaque(panel, event_class)
    state = _case_state(panel, event_class)
    marker_paths = _case_tree_paths(panel, event_class, marker_cache)
    trees = _marker_trees(marker_paths)
    main = panel["main_event"]
    sentinel = panel["sentinel_event"]
    main_key = (main["marker"], main["observed_record_id"])
    sentinel_key = (sentinel["marker"], sentinel["observed_record_id"])
    truth_reference_removals = (
        {main_key}
        if event_class in {"gene_rich_contaminant", "solo_marker_contaminant"}
        else set()
    )
    if truth_reference_removals:
        truth_reference_path, truth_reference_built = _species_tree(
            panel,
            state,
            case_id=case_id,
            arm="truth_reference",
            removals=truth_reference_removals,
            outdir=outdir,
            threads=threads,
        )
    else:
        truth_reference_path = clean_species_path
        truth_reference_built = False
    initial_species_path, initial_built = _species_tree(
        panel,
        state,
        case_id=case_id,
        arm="initial",
        removals=set(),
        outdir=outdir,
        threads=threads,
    )

    loo_rows = score_loo_profiles(trees)
    rows_by_key = {
        (str(row["marker_name"]), str(row["leaf_name"])): row
        for row in loo_rows
    }
    initial_species_tree = Tree(str(initial_species_path))
    cmtv_marker_rf_distances = {
        marker: _rf_distance_between(initial_species_tree, tree)
        for marker, tree in trees.items()
    }
    cmtv_marker_rf_weights = cmtv_rf_quality_weights(
        cmtv_marker_rf_distances
    )
    cmtv_rows = _score_cmtv_current(
        trees,
        loo_rows,
        marker_rf_weights=cmtv_marker_rf_weights,
    )
    raw_proposals = {
        "loo": _loo_proposals(loo_rows),
        "cmtv_weighted": _cmtv_proposals(cmtv_rows),
        "recipient_consensus": _recipient_consensus_proposals(
            trees,
            marker_paths,
            initial_species_path,
        ),
    }
    raw_candidate_keys = {
        (str(proposal["marker_name"]), str(proposal["leaf_name"]))
        for proposals in raw_proposals.values()
        for proposal in proposals
    }
    review_keys = {
        (str(row["marker_name"]), str(row["leaf_name"]))
        for row in loo_rows
        if row["loo_review_candidate"]
    }
    common_audits = _common_attachment_audits(
        panel,
        state,
        trees,
        rows_by_key,
        raw_candidate_keys | review_keys,
    )
    pipelines = {}
    for scorer, proposals in raw_proposals.items():
        evidence_passed = [
            proposal
            for proposal in proposals
            if common_audits[
                (
                    str(proposal["marker_name"]),
                    str(proposal["leaf_name"]),
                )
            ]["gate"]["contig_gate_pass"]
        ]
        budgeted = _select_with_budget(evidence_passed, panel)
        removals, guards = _guard_proposals(
            budgeted,
            marker_paths=marker_paths,
            species_tree_path=initial_species_path,
            mode="loo_profile",
        )
        pipelines[scorer] = {
            "raw": proposals,
            "budgeted": budgeted,
            "removals": removals,
            "guards": guards,
        }

    species_builds = int(initial_built) + int(truth_reference_built)
    final_paths = {}
    for arm, pipeline in pipelines.items():
        removals = pipeline["removals"]
        if removals and removals == truth_reference_removals:
            final_paths[arm] = truth_reference_path
        elif removals:
            path, built = _species_tree(
                panel,
                state,
                case_id=case_id,
                arm=arm,
                removals=removals,
                outdir=outdir,
                threads=threads,
            )
            species_builds += int(built)
            final_paths[arm] = path
        else:
            final_paths[arm] = initial_species_path

    truth_keys = (
        {main_key}
        if event_class in {"gene_rich_contaminant", "solo_marker_contaminant"}
        else set()
    )
    gene_rich_truth = (
        {main_key}
        if event_class == "gene_rich_contaminant"
        else set()
    )
    review_margin_stats = {
        key: _informative_margin_stats(common_audits[key]["evidence"]["votes"])
        for key in review_keys
    }
    review_warned_keys = {
        key
        for key in review_keys
        if common_audits[key]["gate"]["contig_gate_pass"]
        and review_margin_stats[key][0] is not None
        and review_margin_stats[key][0] >= REVIEW_MIN_VOTE_MARGIN
    }
    review_rows = []
    for key in sorted(review_keys):
        row = rows_by_key[key]
        gate = common_audits[key]["gate"]
        evidence = common_audits[key]["evidence"]
        margin_mean, margin_min = review_margin_stats[key]
        review_rows.append(
            {
                "case_id": case_id,
                "panel_id": panel["panel_id"],
                "regime": panel["regime"],
                "seed": panel["seed"],
                "event_class": event_class,
                "distance_stratum": main["distance_stratum"],
                "marker_name": key[0],
                "record_id": key[1],
                "genome": str(row["genome"]),
                "loo_robust_z": row["loo_robust_z"],
                "loo_target_discordance": row["loo_target_discordance"],
                "loo_voter_upper": row["loo_voter_upper"],
                "loo_target_support": row["loo_target_support"],
                "loo_voter_search_mode": row["loo_voter_search_mode"],
                "evidence_status": evidence.get("input_status"),
                "informative_gene_count": gate.get("informative_gene_count", 0),
                "contig_gate_pass": gate.get("contig_gate_pass", False),
                "contig_abstention_reason": gate.get("contig_abstention_reason"),
                "vote_margin_mean": margin_mean,
                "vote_margin_min": margin_min,
                "margin_pass": (
                    margin_mean is not None
                    and margin_mean >= REVIEW_MIN_VOTE_MARGIN
                ),
                "review_warning": key in review_warned_keys,
                "is_truth_event": key in truth_keys,
                "is_gene_rich_truth": key in gene_rich_truth,
            }
        )
    method_counts = {}
    for method, pipeline in pipelines.items():
        removals = pipeline["removals"]
        method_counts.update(
            {
                f"{method}_truth_positive_count": len(removals & truth_keys),
                f"{method}_truth_false_positive_count": len(removals - truth_keys),
                f"{method}_gene_rich_true_positive_count": len(
                    removals & gene_rich_truth
                ),
                f"{method}_gene_rich_false_positive_count": (
                    len(removals - gene_rich_truth)
                    if event_class == "gene_rich_contaminant"
                    else 0
                ),
                f"{method}_unsupported_action_count": (
                    len(removals - gene_rich_truth)
                    if event_class == "gene_rich_contaminant"
                    else len(removals)
                ),
                f"{method}_sentinel_removed": sentinel_key in removals,
            }
        )
    observed_record_keys = {
        (marker, record_id)
        for marker, per_genome in state["records"].items()
        for record_id in per_genome.values()
    }
    loo_removals = pipelines["loo"]["removals"]
    unaffected_record_keys = observed_record_keys - truth_keys
    loo_retained_record_keys = observed_record_keys - loo_removals
    loo_unaffected_losses = sorted(
        unaffected_record_keys - loo_retained_record_keys
    )
    loo_remaining_markers_by_genome = {
        genome: sum(
            (marker, state["records"][marker][genome]) not in loo_removals
            for marker in panel["markers"]
        )
        for genome in panel["genomes"]
    }
    final_taxa = {
        arm: sorted(str(leaf.name) for leaf in Tree(str(path)).iter_leaves())
        for arm, path in final_paths.items()
    }

    clean_target = panel["marker_records"][main["marker"]][main["recipient"]]
    cmtv_k5_preserved = _nearest_genome_set(
        Tree(str(marker_cache[(main["marker"], "clean")]), format=1),
        clean_target,
        5,
    ) == _nearest_genome_set(
        Tree(str(marker_cache[(main["marker"], "source")]), format=1),
        main["observed_record_id"],
        5,
    )

    comparison_rows = _comparison_rows(
        case_id=case_id,
        panel=panel,
        event_class=event_class,
        events=state["events"],
        rows_by_key=rows_by_key,
        audits=common_audits,
        pipelines=pipelines,
    )
    comparison_by_event = {
        (row["marker_name"], row["record_id"], row["scorer"]): row
        for row in comparison_rows
        if row["is_declared_event"]
    }
    events = []
    for event in state["events"]:
        key = (event["marker"], event["observed_record_id"])
        row = rows_by_key[key]
        audit = common_audits.get(key, {})
        event_result = {
            "event_id": event["event_id"],
            "event_kind": event["event_kind"],
            "is_contaminant": event["is_contaminant"],
            "marker_name": event["marker"],
            "record_id": event["observed_record_id"],
            "recipient": event["recipient"],
            "anchor": event["anchor"],
            "loo_class": row["loo_class"],
            "loo_abstention_reason": row["loo_abstention_reason"],
            "loo_score": row["loo_score"],
            "loo_target_support": row["loo_target_support"],
            "loo_attachment_clade": row["loo_attachment_clade"],
            "sequence_evidence_status": audit.get("evidence", {}).get(
                "input_status"
            ),
            "sequence_informative_vote_count": audit.get("evidence", {}).get(
                "informative_vote_count", 0
            ),
            "contig_gate_pass": audit.get("gate", {}).get(
                "contig_gate_pass", False
            ),
            "contig_abstention_reason": audit.get("gate", {}).get(
                "contig_abstention_reason"
            ),
            "loo_robust_z": row["loo_robust_z"],
            "loo_review_candidate": row["loo_review_candidate"],
            "loo_review_warning": key in review_warned_keys,
            "loo_voter_search_mode": row["loo_voter_search_mode"],
        }
        for scorer in pipelines:
            comparison = comparison_by_event[(*key, scorer)]
            event_result.update(
                {
                    f"{scorer}_raw_call": comparison["raw_call"],
                    f"{scorer}_budget_pass": comparison["budget_pass"],
                    f"{scorer}_removed": comparison["final_action"] == "removed",
                    f"{scorer}_reason": comparison["reason"],
                }
            )
        events.append(event_result)
    if event_class == "gene_rich_contaminant":
        evaluation_role = (
            "primary_far_source_positive"
            if main["distance_stratum"] == "far"
            else "near_source_identifiability_control"
        )
        primary_positive = main["distance_stratum"] == "far"
    elif event_class == "clean":
        evaluation_role = "clean_negative_control"
        primary_positive = False
    else:
        evaluation_role = "solo_marker_abstention_control"
        primary_positive = False
    return {
        "case_id": case_id,
        "panel_id": panel["panel_id"],
        "regime": panel["regime"],
        "seed": panel["seed"],
        "event_class": event_class,
        "distance_stratum": main["distance_stratum"],
        "evaluation_role": evaluation_role,
        "primary_positive": primary_positive,
        "taxa_count": TAXA_COUNT,
        "marker_count": len(panel["markers"]),
        "truth_reference_removed_count": len(truth_reference_removals),
        "truth_reference_preserves_sentinel": sentinel_key not in truth_reference_removals,
        "cmtv_k5_neighbors_preserved_posthoc": cmtv_k5_preserved,
        "cmtv_marker_rf_distances": cmtv_marker_rf_distances,
        "cmtv_marker_rf_weights": cmtv_marker_rf_weights,
        "loo_review_candidate_count": len(review_keys),
        "loo_review_warning_count": len(review_warned_keys),
        "loo_review_warning_truth_count": len(review_warned_keys & truth_keys),
        "loo_review_warning_false_count": len(review_warned_keys - truth_keys),
        "_review_rows": review_rows,
        **{
            f"{scorer}_raw_call_count": len(pipeline["raw"])
            for scorer, pipeline in pipelines.items()
        },
        **{
            f"{scorer}_removed_count": len(pipeline["removals"])
            for scorer, pipeline in pipelines.items()
        },
        **method_counts,
        "loo_unaffected_record_loss_count": len(loo_unaffected_losses),
        "loo_unaffected_record_losses": [
            {"marker_name": marker, "record_id": record_id}
            for marker, record_id in loo_unaffected_losses
        ],
        "loo_unknown_removal_count": len(loo_removals - observed_record_keys),
        "loo_min_remaining_marker_count": min(
            loo_remaining_markers_by_genome.values()
        ),
        "loo_final_taxa": final_taxa["loo"],
        "loo_final_taxa_match": final_taxa["loo"] == sorted(panel["genomes"]),
        "cmtv_weighted_final_taxa_match": final_taxa["cmtv_weighted"] == sorted(
            panel["genomes"]
        ),
        "recipient_consensus_final_taxa_match": final_taxa[
            "recipient_consensus"
        ] == sorted(panel["genomes"]),
        "initial_rf_norm": _rf_norm(truth_reference_path, initial_species_path),
        "loo_rf_norm": _rf_norm(truth_reference_path, final_paths["loo"]),
        "cmtv_weighted_rf_norm": _rf_norm(
            truth_reference_path,
            final_paths["cmtv_weighted"],
        ),
        "recipient_consensus_rf_norm": _rf_norm(
            truth_reference_path,
            final_paths["recipient_consensus"],
        ),
        "initial_patristic_nrmse": _patristic_nrmse(
            truth_reference_path,
            initial_species_path,
        ),
        "loo_patristic_nrmse": _patristic_nrmse(
            truth_reference_path,
            final_paths["loo"],
        ),
        "cmtv_weighted_patristic_nrmse": _patristic_nrmse(
            truth_reference_path,
            final_paths["cmtv_weighted"],
        ),
        "recipient_consensus_patristic_nrmse": _patristic_nrmse(
            truth_reference_path,
            final_paths["recipient_consensus"],
        ),
        **{
            f"{scorer}_guard_decisions": pipeline["guards"]
            for scorer, pipeline in pipelines.items()
        },
        "comparison_rows": comparison_rows,
        "events": events,
    }, species_builds


def _write_frozen_manifest(outdir: Path, manifest: dict) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "benchmark_manifest.json"
    if path.exists():
        existing = json.loads(path.read_text())
        if existing != manifest:
            raise ValueError(
                f"existing benchmark manifest differs from the fixed design: {path}"
            )
    else:
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return path


def _interval(successes: int, total: int, unit: str) -> dict:
    return {**wilson_interval(successes, total), "unit": unit}


def _f1(true_positives: int, false_positives: int, false_negatives: int) -> float:
    denominator = (2 * true_positives) + false_positives + false_negatives
    return (2 * true_positives / denominator) if denominator else 0.0


def _aggregate_method_metrics(
    cases: list[dict],
    *,
    label: str,
    field_prefix: str,
) -> dict:
    gene_cases = [
        case for case in cases if case["event_class"] == "gene_rich_contaminant"
    ]
    far_gene_cases = [
        case for case in gene_cases if case["distance_stratum"] == "far"
    ]
    near_gene_cases = [
        case for case in gene_cases if case["distance_stratum"] == "near"
    ]
    clean_cases = [case for case in cases if case["event_class"] == "clean"]
    solo_cases = [
        case for case in cases if case["event_class"] == "solo_marker_contaminant"
    ]
    contaminated_cases = gene_cases + solo_cases
    removed_key = f"{field_prefix}_removed_count"
    gene_tp = sum(
        case[f"{field_prefix}_gene_rich_true_positive_count"]
        for case in gene_cases
    )
    gene_fp = sum(
        case[f"{field_prefix}_gene_rich_false_positive_count"]
        for case in gene_cases
    )
    truth_tp = sum(
        case[f"{field_prefix}_truth_positive_count"]
        for case in cases
    )
    truth_fp = sum(
        case[f"{field_prefix}_truth_false_positive_count"]
        for case in cases
    )
    removals = sum(case[removed_key] for case in cases)
    unsupported = sum(
        case[f"{field_prefix}_unsupported_action_count"]
        for case in cases
    )
    return {
        f"{label}_gene_rich_recall": _interval(
            gene_tp,
            len(gene_cases),
            "held_out_regime_seed_panel",
        ),
        f"{label}_gene_rich_marker_precision": _interval(
            gene_tp,
            gene_tp + gene_fp,
            "guarded_gene_rich_removal",
        ),
        f"{label}_far_source_recall": _interval(
            sum(
                case[f"{field_prefix}_gene_rich_true_positive_count"]
                for case in far_gene_cases
            ),
            len(far_gene_cases),
            "truth_topology_far_source_panel",
        ),
        f"{label}_near_source_abstention": _interval(
            sum(
                case[f"{field_prefix}_gene_rich_true_positive_count"] == 0
                for case in near_gene_cases
            ),
            len(near_gene_cases),
            "truth_topology_near_source_event",
        ),
        f"{label}_truth_marker_precision": _interval(
            truth_tp,
            removals,
            "guarded_removal_against_event_truth",
        ),
        f"{label}_truth_marker_recall": _interval(
            truth_tp,
            len(contaminated_cases),
            "held_out_panel_context",
        ),
        f"{label}_clean_safety": _interval(
            sum(case[removed_key] == 0 for case in clean_cases),
            len(clean_cases),
            "held_out_clean_panel",
        ),
        f"{label}_solo_abstention": _interval(
            sum(case[removed_key] == 0 for case in solo_cases),
            len(solo_cases),
            "held_out_solo_panel",
        ),
        f"{label}_sentinel_safety": _interval(
            sum(
                not case[f"{field_prefix}_sentinel_removed"]
                for case in gene_cases
            ),
            len(gene_cases),
            "held_out_native_contig_sentinel",
        ),
        f"{label}_gene_rich_f1": _f1(
            gene_tp,
            gene_fp,
            len(gene_cases) - gene_tp,
        ),
        f"{label}_gene_rich_true_positive_count": gene_tp,
        f"{label}_gene_rich_false_positive_count": gene_fp,
        f"{label}_truth_positive_count": truth_tp,
        f"{label}_truth_false_positive_count": truth_fp,
        f"{label}_unsupported_action_count": unsupported,
        f"{label}_removal_count": removals,
    }


def _raw_discrimination_metrics(comparison_rows: list[dict]) -> dict:
    metrics = {}
    for scorer in ("loo", "cmtv_weighted", "recipient_consensus"):
        rows = [
            row
            for row in comparison_rows
            if row["scorer"] == scorer
            and row["event_class"] == "gene_rich_contaminant"
        ]
        positives = sum(
            row["event_kind"] == "source_replacement"
            and row["is_contaminant"]
            for row in rows
        )
        true_positives = sum(
            row["raw_call"]
            and row["event_kind"] == "source_replacement"
            and row["is_contaminant"]
            for row in rows
        )
        false_positives = sum(
            row["raw_call"]
            and not (
                row["event_kind"] == "source_replacement"
                and row["is_contaminant"]
            )
            for row in rows
        )
        false_negatives = positives - true_positives
        metrics[scorer] = {
            "true_positive_count": true_positives,
            "false_positive_count": false_positives,
            "false_negative_count": false_negatives,
            "precision": (
                true_positives / (true_positives + false_positives)
                if true_positives + false_positives
                else 0.0
            ),
            "recall": (
                true_positives / positives
                if positives
                else 0.0
            ),
            "f1": _f1(true_positives, false_positives, false_negatives),
        }
    return metrics


def _wp1_decision(raw_metrics: dict) -> dict:
    scores = {
        scorer: float(values["f1"])
        for scorer, values in raw_metrics.items()
    }
    best_score = max(scores.values())
    tied = sorted(
        scorer
        for scorer, score in scores.items()
        if abs(score - best_score) <= 1e-12
    )
    if scores["loo"] < best_score - 1e-12:
        return {
            "status": "stop_loo_worse",
            "selected_scorer": None,
            "empirical_confirmation_ready": False,
            "criterion": "gene_rich_raw_call_f1",
            "tie_break_order": list(SCORER_SIMPLICITY_ORDER),
            "scores": scores,
            "tied_scorers": tied,
            "rationale": "LOO has worse raw candidate discrimination than another scorer.",
        }
    if tied == ["loo"]:
        return {
            "status": "retain_loo",
            "selected_scorer": "loo",
            "empirical_confirmation_ready": True,
            "criterion": "gene_rich_raw_call_f1",
            "tie_break_order": list(SCORER_SIMPLICITY_ORDER),
            "scores": scores,
            "tied_scorers": tied,
            "rationale": "LOO has the best raw candidate discrimination.",
        }
    selected = next(
        scorer for scorer in SCORER_SIMPLICITY_ORDER if scorer in tied
    )
    return {
        "status": "select_simpler_tied_scorer",
        "selected_scorer": selected,
        "empirical_confirmation_ready": True,
        "criterion": "gene_rich_raw_call_f1",
        "tie_break_order": list(SCORER_SIMPLICITY_ORDER),
        "scores": scores,
        "tied_scorers": tied,
        "rationale": (
            "LOO ties for best raw candidate discrimination; "
            f"{selected} is first in the predeclared simplicity order "
            "(marker-tree-only LOO, weighted CMTV, then alignment-aware "
            "recipient_consensus)."
        ),
    }


def _write_comparison_table(outdir: Path, rows: list[dict]) -> Path:
    path = outdir / "per_event_comparison.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


REVIEW_TIER_FIELDS = [
    "case_id",
    "panel_id",
    "regime",
    "seed",
    "event_class",
    "distance_stratum",
    "marker_name",
    "record_id",
    "genome",
    "loo_robust_z",
    "loo_target_discordance",
    "loo_voter_upper",
    "loo_target_support",
    "loo_voter_search_mode",
    "evidence_status",
    "informative_gene_count",
    "contig_gate_pass",
    "contig_abstention_reason",
    "vote_margin_mean",
    "vote_margin_min",
    "margin_pass",
    "review_warning",
    "is_truth_event",
    "is_gene_rich_truth",
]


def _write_review_tier_table(outdir: Path, rows: list[dict]) -> Path:
    path = outdir / "review_tier.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REVIEW_TIER_FIELDS,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


SWEEP_FIELDS = [
    "regime",
    "seed",
    "distance_stratum",
    "donor_gene_count",
    "evidence_status",
    "informative_gene_count",
    "agreement_count",
    "strongest_conflict_count",
    "contig_gate_pass",
    "contig_abstention_reason",
]


def run_donor_gene_sweep(
    outdir: Path,
    *,
    threads: int = 1,
    panel_filter: set[tuple[str, int]] | None = None,
) -> dict:
    """Trace the contig-gate operating curve over donor non-marker gene counts.

    Gate-stage only: marker trees come from the shared cache and the truth
    row is scored once per panel; only the donor candidate gene set varies.
    """
    if threads < 1 or threads > MAX_THREADS:
        raise ValueError(f"threads must be between 1 and {MAX_THREADS}")
    rows = []
    for panel in build_panels():
        coordinate = (panel["regime"], panel["seed"])
        if panel_filter is not None and coordinate not in panel_filter:
            continue
        cache, _built = _build_marker_trees(panel, outdir, threads)
        trees = _marker_trees(_case_tree_paths(panel, "gene_rich_contaminant", cache))
        loo_rows = score_loo_profiles(trees)
        main = panel["main_event"]
        main_key = (main["marker"], main["observed_record_id"])
        row = next(
            item
            for item in loo_rows
            if (str(item["marker_name"]), str(item["leaf_name"])) == main_key
        )
        reference_proteomes = _reference_proteomes(panel)
        marker_record_ids = _all_marker_record_ids(panel)
        attachment = set(row["loo_attachment_taxa"])
        background = {
            str(leaf.name).split("|", 1)[0]
            for leaf in trees[str(row["marker_name"])].iter_leaves()
        } - {main["recipient"]} - attachment
        donor_genes = panel["sweep_donor_gene_sequences"]
        stratum = SOURCE_DISTANCE_STRATA[coordinate]
        for count in SWEEP_DONOR_GENE_COUNTS:
            candidate_genes = {
                build_sequence_id(
                    main["recipient"],
                    main["contig_id"],
                    _opaque("q", panel["panel_id"], "sweep", index),
                ): sequence
                for index, sequence in enumerate(donor_genes[:count])
            }
            evidence = assign_contig_gene_split_votes(
                candidate_genes,
                reference_proteomes,
                recipient_genome=main["recipient"],
                candidate_contig_id=main["contig_id"],
                marker_record_ids=marker_record_ids,
                attachment_taxa=attachment,
                background_taxa=background,
            )
            gate = contig_gene_vote_gate(evidence["votes"], row["loo_attachment_clade"])
            rows.append(
                {
                    "regime": panel["regime"],
                    "seed": panel["seed"],
                    "distance_stratum": stratum,
                    "donor_gene_count": count,
                    "evidence_status": evidence.get("input_status"),
                    "informative_gene_count": gate.get("informative_gene_count", 0),
                    "agreement_count": gate.get("agreement_count", 0),
                    "strongest_conflict_count": gate.get(
                        "strongest_conflict_count", 0
                    ),
                    "contig_gate_pass": gate.get("contig_gate_pass", False),
                    "contig_abstention_reason": gate.get("contig_abstention_reason"),
                }
            )
    path = outdir / "donor_gene_sweep.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SWEEP_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return {"rows": rows, "path": str(path)}


def run_benchmark(outdir: Path, *, threads: int = 1) -> dict:
    if threads < 1 or threads > MAX_THREADS:
        raise ValueError(f"threads must be between 1 and {MAX_THREADS}")
    started = time.monotonic()
    panels = build_panels()
    manifest = build_manifest(panels)
    manifest_path = _write_frozen_manifest(outdir, manifest)

    marker_caches = {}
    marker_tree_builds_this_run = 0
    distinct_marker_tree_count = 0
    for panel in panels:
        cache, built = _build_marker_trees(panel, outdir, threads)
        marker_caches[panel["panel_id"]] = cache
        marker_tree_builds_this_run += built
        distinct_marker_tree_count += len(cache)

    reference_species_paths = {}
    species_tree_builds_this_run = 0
    for panel in panels:
        clean_case_id = _case_opaque(panel, "clean")
        clean_path, built = _species_tree(
            panel,
            _case_state(panel, "clean"),
            case_id=clean_case_id,
            arm="initial",
            removals=set(),
            outdir=outdir,
            threads=threads,
        )
        reference_species_paths[panel["panel_id"]] = clean_path
        species_tree_builds_this_run += int(built)

    cases = []
    for panel in panels:
        for event_class in EVENT_CLASSES:
            case, builds = _score_case(
                panel,
                event_class,
                marker_cache=marker_caches[panel["panel_id"]],
                outdir=outdir,
                threads=threads,
                clean_species_path=reference_species_paths[panel["panel_id"]],
            )
            species_tree_builds_this_run += builds
            cases.append(case)

    comparison_rows = [
        row
        for case in cases
        for row in case["comparison_rows"]
    ]
    comparison_path = _write_comparison_table(outdir, comparison_rows)
    # Greedy voter-selection tier: one far-source 16-marker panel exercising
    # the production greedy path. Scored outside the frozen 12-case matrix
    # with its own tree accounting.
    greedy_panel = build_panel(
        "alpha_like",
        809,
        marker_count=GREEDY_MARKER_COUNT,
        tag=GREEDY_PANEL_TAG,
    )
    greedy_outdir = outdir / "greedy_tier"
    greedy_cache, greedy_marker_builds = _build_marker_trees(
        greedy_panel, greedy_outdir, threads
    )
    greedy_clean_path, greedy_clean_built = _species_tree(
        greedy_panel,
        _case_state(greedy_panel, "clean"),
        case_id=_case_opaque(greedy_panel, "clean"),
        arm="initial",
        removals=set(),
        outdir=greedy_outdir,
        threads=threads,
    )
    greedy_case, greedy_species_builds = _score_case(
        greedy_panel,
        "gene_rich_contaminant",
        marker_cache=greedy_cache,
        outdir=greedy_outdir,
        threads=threads,
        clean_species_path=greedy_clean_path,
    )
    greedy_case.pop("_review_rows")
    greedy_tier = {
        "panel_tag": GREEDY_PANEL_TAG,
        "marker_count": GREEDY_MARKER_COUNT,
        "distance_stratum": greedy_case["distance_stratum"],
        "marker_tree_builds": greedy_marker_builds,
        "distinct_marker_trees": len(greedy_cache),
        "species_tree_builds": greedy_species_builds + int(greedy_clean_built),
        "voter_search_modes": sorted(
            {
                str(event["loo_voter_search_mode"])
                for event in greedy_case["events"]
                if event.get("loo_voter_search_mode") is not None
            }
        ),
        "case": greedy_case,
    }

    review_rows = [row for case in cases for row in case.pop("_review_rows")]
    review_path = _write_review_tier_table(outdir, review_rows)
    review_tier = {
        "review_robust_z": REVIEW_ROBUST_Z,
        "review_min_vote_margin": REVIEW_MIN_VOTE_MARGIN,
        "candidate_count": len(review_rows),
        "gate_only_warning_count": sum(
            row["contig_gate_pass"] for row in review_rows
        ),
        "gate_only_false_count": sum(
            row["contig_gate_pass"] and not row["is_truth_event"]
            for row in review_rows
        ),
        "warning_count": sum(row["review_warning"] for row in review_rows),
        "warning_truth_count": sum(
            row["review_warning"] and row["is_truth_event"] for row in review_rows
        ),
        "warning_false_count": sum(
            row["review_warning"] and not row["is_truth_event"]
            for row in review_rows
        ),
        "near_source_gene_rich_warned": sum(
            row["review_warning"]
            and row["is_gene_rich_truth"]
            and row["distance_stratum"] == "near"
            for row in review_rows
        ),
        "near_source_gene_rich_total": sum(
            case["event_class"] == "gene_rich_contaminant"
            and case["distance_stratum"] == "near"
            for case in cases
        ),
    }
    raw_discrimination = _raw_discrimination_metrics(comparison_rows)
    decision = _wp1_decision(raw_discrimination)
    metrics = {}
    for label, field_prefix in (
        ("loo", "loo"),
        ("cmtv_weighted", "cmtv_weighted"),
        ("recipient_consensus", "recipient_consensus"),
    ):
        metrics.update(
            _aggregate_method_metrics(
                cases,
                label=label,
                field_prefix=field_prefix,
            )
        )
    gene_cases = [
        case for case in cases if case["event_class"] == "gene_rich_contaminant"
    ]
    far_gene_cases = [
        case for case in gene_cases if case["distance_stratum"] == "far"
    ]
    guarded_loo_removals = [
        decision
        for case in cases
        for decision in case["loo_guard_decisions"]
        if decision["removed_leaf_names"]
    ]
    metrics.update(
        {
            "species_topology_rf_nonzero_case_count": sum(
                any(
                    case[key] > 0.0
                    for key in (
                        "initial_rf_norm",
                        "loo_rf_norm",
                        "cmtv_weighted_rf_norm",
                        "recipient_consensus_rf_norm",
                    )
                )
                for case in cases
            ),
            "gene_rich_initial_patristic_nonzero_count": sum(
                case["initial_patristic_nrmse"] > 0.0
                for case in gene_cases
            ),
            "loo_gene_rich_patristic_nonworsening": _interval(
                sum(
                    case["loo_patristic_nrmse"]
                    <= case["initial_patristic_nrmse"] + 1e-12
                    for case in gene_cases
                ),
                len(gene_cases),
                "held_out_regime_seed_panel",
            ),
            "loo_gene_rich_patristic_improvement": _interval(
                sum(
                    case["loo_patristic_nrmse"]
                    < case["initial_patristic_nrmse"] - 1e-12
                    for case in gene_cases
                ),
                len(gene_cases),
                "held_out_regime_seed_panel",
            ),
            "loo_far_source_patristic_improvement": _interval(
                sum(
                    case["loo_patristic_nrmse"]
                    < case["initial_patristic_nrmse"] - 1e-12
                    for case in far_gene_cases
                ),
                len(far_gene_cases),
                "truth_topology_far_source_panel",
            ),
            "loo_marker_rf_guard_improvement": _interval(
                sum(
                    decision["rf_after"] < decision["rf_before"]
                    for decision in guarded_loo_removals
                ),
                len(guarded_loo_removals),
                "guarded_marker_removal",
            ),
        }
    )
    return {
        "benchmark_kind": manifest["benchmark_kind"],
        "primary_estimand": manifest["primary_estimand"],
        "near_source_role": manifest["near_source_role"],
        "manifest_path": str(manifest_path),
        "production_pruning_enabled": False,
        "manifest_written_before_tree_inference": True,
        "fair_pipeline_stages": [
            "raw_call",
            "common_attachment_and_contig_audit",
            "candidate_budget",
            "rf_guard",
        ],
        "compared_scorers": [
            "loo",
            "cmtv_weighted",
            "recipient_consensus",
        ],
        "cmtv_baseline": "CMTV-current with positive RF-quality voter weights",
        "raw_candidate_discrimination": raw_discrimination,
        "wp1_decision": decision,
        "review_tier": review_tier,
        "review_tier_path": str(review_path),
        "greedy_tier": greedy_tier,
        "per_event_comparison_path": str(comparison_path),
        "per_event_comparison_rows": comparison_rows,
        "case_count": len(cases),
        "panel_count": len(panels),
        "marker_tree_instance_count": len(cases) * MARKER_COUNT,
        "distinct_marker_tree_build_count": distinct_marker_tree_count,
        "marker_tree_builds_this_run": marker_tree_builds_this_run,
        "species_tree_builds_this_run": species_tree_builds_this_run,
        "threads": threads,
        "runtime_seconds": time.monotonic() - started,
        "metrics": metrics,
        "limitations": [
            "simulated aligned proteins do not test gene calling, HMM search, or alignment",
            "four held-out positives form an engineering gate, not a significance test",
            "phmmer split votes are benchmark-only; production LOO remains report-only",
            "CMTV voter weights are positive one-minus-normalized-RF values from each case's initial species tree",
            "source-distance strata and anchors use only the truth topology; CMTV k=5 behavior is analyzed post hoc",
            "near-source contamination that stays within clean-marker dispersion is an identifiability control and remains unpruned",
            "gene-rich and solo contexts reuse four source events and are not independent biological samples",
            "species-tree topology RF may be insensitive at this scale; normalized patristic RMSE supplies a branch-length-sensitive secondary check",
        ],
        "cases": cases,
    }


def check_benchmark(report: dict) -> list[str]:
    errors = []
    expected_coordinates = {
        (regime, seed, event_class)
        for regime in REGIMES
        for seed in SEEDS
        for event_class in EVENT_CLASSES
    }
    coordinates = {
        (case["regime"], case["seed"], case["event_class"])
        for case in report["cases"]
    }
    if coordinates != expected_coordinates or report["case_count"] != 12:
        errors.append("case matrix is not exactly 2 regimes x 2 seeds x 3 classes")
    if report["marker_tree_instance_count"] != EXPECTED_MARKER_TREE_INSTANCES:
        errors.append("case matrix does not contain exactly 96 marker-tree instances")
    if report["distinct_marker_tree_build_count"] != EXPECTED_DISTINCT_MARKER_TREE_BUILDS:
        errors.append("shared cache does not contain exactly 40 distinct marker trees")
    if report["distinct_marker_tree_build_count"] > MAX_MARKER_TREE_BUILDS:
        errors.append("marker-tree build ceiling exceeded")
    if report["threads"] > MAX_THREADS:
        errors.append("thread ceiling exceeded")
    if report["production_pruning_enabled"]:
        errors.append("benchmark enabled production pruning")
    if report.get("compared_scorers") != [
        "loo",
        "cmtv_weighted",
        "recipient_consensus",
    ]:
        errors.append("fair scorer set is incomplete or reordered")
    if report.get("fair_pipeline_stages") != [
        "raw_call",
        "common_attachment_and_contig_audit",
        "candidate_budget",
        "rf_guard",
    ]:
        errors.append("fair pipeline stages are not recorded in the required order")

    gene_cases = [
        case for case in report["cases"]
        if case["event_class"] == "gene_rich_contaminant"
    ]
    for case in gene_cases:
        event_by_kind = {event["event_kind"]: event for event in case["events"]}
        main = event_by_kind.get("source_replacement", {})
        sentinel = event_by_kind.get("native_contig_sentinel", {})
        if main.get("loo_removed") and (
            main.get("loo_class") != "discordant_marker"
            or main.get("sequence_informative_vote_count", 0) < 3
            or not main.get("contig_gate_pass")
        ):
            errors.append(f"source removal bypassed a required gate in {case['case_id']}")
        if sentinel.get("loo_class") != "discordant_marker":
            errors.append(f"LOO missed native-contig sentinel topology in {case['case_id']}")
        if sentinel.get("contig_gate_pass"):
            errors.append(f"native-contig sentinel passed evidence gate in {case['case_id']}")
    for case in report["cases"]:
        if len(case["cmtv_marker_rf_weights"]) != MARKER_COUNT or any(
            weight <= 0.0
            for weight in case["cmtv_marker_rf_weights"].values()
        ):
            errors.append(f"CMTV weights are missing or nonpositive in {case['case_id']}")
        for scorer in report["compared_scorers"]:
            if not case[f"{scorer}_final_taxa_match"]:
                errors.append(f"{scorer} lost a final taxon in {case['case_id']}")

    allowed_reasons = {
        "no_raw_call",
        "attachment_abstention",
        "evidence_abstention",
        "budget_rejected",
        "rf_guard_rejected",
        "removed",
    }
    comparison_rows = report.get("per_event_comparison_rows", [])
    comparison_keys = {
        (
            row["case_id"],
            row["marker_name"],
            row["record_id"],
            row["scorer"],
        )
        for row in comparison_rows
    }
    if len(comparison_keys) != len(comparison_rows):
        errors.append("per-event comparison table contains duplicate scorer rows")
    if {row["scorer"] for row in comparison_rows} != set(report["compared_scorers"]):
        errors.append("per-event comparison table omits a scorer")
    for row in comparison_rows:
        if row["reason"] not in allowed_reasons:
            errors.append(f"unknown terminal reason in {row['case_id']}")
        if row["reason"] == "no_raw_call" and row["raw_call"]:
            errors.append(f"raw call has no-call reason in {row['case_id']}")
        if row["budget_pass"] and (
            not row["raw_call"] or not row["contig_gate_pass"]
        ):
            errors.append(f"candidate bypassed evidence before budget in {row['case_id']}")
        if row["rf_guard_decision"] is not None and not row["budget_pass"]:
            errors.append(f"candidate bypassed budget before RF guard in {row['case_id']}")
        if row["final_action"] == "removed" and (
            row["reason"] != "removed"
            or not row["raw_call"]
            or not row["contig_gate_pass"]
            or not row["budget_pass"]
            or row["rf_guard_decision"] != "pruned"
        ):
            errors.append(f"removal bypassed the fair gate chain in {row['case_id']}")

    if report["metrics"]["loo_truth_marker_precision"]["rate"] is None or report["metrics"]["loo_truth_marker_precision"]["rate"] < 0.90:
        errors.append("LOO guarded marker precision is below 0.90")
    if report["metrics"]["loo_far_source_recall"]["rate"] < 0.75:
        errors.append("LOO far-source recall is below 0.75")
    if report["metrics"]["loo_near_source_abstention"]["rate"] != 1.0:
        errors.append("LOO did not abstain on a near-source identifiability control")
    for metric in ("loo_clean_safety", "loo_solo_abstention", "loo_sentinel_safety"):
        if report["metrics"][metric]["rate"] != 1.0:
            errors.append(f"safety gate failed: {metric}")
    if report["metrics"]["loo_marker_rf_guard_improvement"]["rate"] != 1.0:
        errors.append("LOO marker RF guard did not improve every accepted marker tree")
    if report["metrics"]["gene_rich_initial_patristic_nonzero_count"] != len(gene_cases):
        errors.append("a gene-rich case lacks a branch-length-sensitive perturbation")
    if report["metrics"]["loo_gene_rich_patristic_nonworsening"]["rate"] != 1.0:
        errors.append("LOO worsened species-tree patristic error in a gene-rich case")
    if report["metrics"]["loo_far_source_patristic_improvement"]["rate"] != 1.0:
        errors.append("LOO did not improve patristic error for every far-source case")
    expected_decision = _wp1_decision(report["raw_candidate_discrimination"])
    if report.get("wp1_decision") != expected_decision:
        errors.append("WP1 scorer decision does not match the raw comparison rows")
    for scorer in report["compared_scorers"]:
        for metric in ("clean_safety", "solo_abstention", "sentinel_safety"):
            if report["metrics"][f"{scorer}_{metric}"]["rate"] != 1.0:
                errors.append(f"safety gate failed: {scorer}_{metric}")
    if any(
        not case["loo_final_taxa_match"]
        or case["loo_unaffected_record_loss_count"]
        or case["loo_unknown_removal_count"]
        or case["loo_min_remaining_marker_count"] < 1
        for case in report["cases"]
    ):
        errors.append("LOO lost a taxon or unaffected marker record")
    if report["runtime_seconds"] > MAX_RUNTIME_SECONDS:
        errors.append(
            f"runtime {report['runtime_seconds']:.3f}s exceeds {MAX_RUNTIME_SECONDS:.0f}s"
        )
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="runs/us009_sequence")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--sweep-donor-genes", action="store_true")
    args = parser.parse_args(argv)
    if args.sweep_donor_genes:
        sweep = run_donor_gene_sweep(Path(args.outdir), threads=args.threads)
        print(json.dumps(sweep["rows"], indent=2, sort_keys=True))
        print(f"WROTE: {sweep['path']}")
        return 0
    report = run_benchmark(Path(args.outdir), threads=args.threads)
    print(json.dumps(report, indent=2, sort_keys=True))
    if not args.check:
        return 0
    errors = check_benchmark(report)
    if errors:
        print("FAIL: " + "; ".join(errors))
        return 1
    print(
        f"PASS: {report['case_count']} cases, "
        f"{report['marker_tree_instance_count']} marker-tree instances, "
        f"{report['distinct_marker_tree_build_count']} distinct builds, "
        f"{report['runtime_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
