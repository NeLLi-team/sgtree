"""Bounded in-memory tree-fixture screen for US-009.

The screen keeps 24 small mechanism cases and adds eight 50/100-leaf scale
cases without starting an SGTree job. It fixes events and opaque record IDs
before either scorer runs. It runs the LOO and baseline CMTV scoring cores in
memory. Synthetic contig votes and RF checks exercise the decision gates, but
never enable production pruning.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.simplefilter("ignore", SyntaxWarning)

from ete3 import Tree

from sgtree.id_schema import build_sequence_id
from sgtree.marker_selection import _rf_distance_between, choose_tree_by_rf
from sgtree.marker_selection.contig_evidence import contig_gene_vote_gate
from sgtree.marker_selection.loo_profile import score_loo_profiles


SHAPES = {"balanced": 10, "pectinate": 12}
EVENT_CLASSES = (
    "clean",
    "gene_rich_contaminant",
    "solo_marker_contaminant",
    "safety_controls",
)
SEEDS = (101, 202, 303)
SCALE_TAXA_COUNTS = (50, 100)
SCALE_EVENT_CLASSES = ("scale_clean", "scale_far_source")
SCALE_SEED = 101
MARKER_COUNT = 8
MAX_RUNTIME_SECONDS = 30.0
MIN_CMTV_RF_WEIGHT = 1e-6


def cmtv_rf_quality_weights(
    marker_rf_distances: dict[str, float],
) -> dict[str, float]:
    """Convert normalized RF distances to positive CMTV voter weights."""
    weights = {}
    for marker, value in marker_rf_distances.items():
        rf_distance = float(value)
        if not np.isfinite(rf_distance) or not 0.0 <= rf_distance <= 1.0:
            raise ValueError(
                f"invalid normalized RF distance for marker {marker!r}: {value!r}"
            )
        weights[str(marker)] = max(MIN_CMTV_RF_WEIGHT, 1.0 - rf_distance)
    return weights


def load_cmtv_rf_weights(path: str | Path) -> dict[str, float]:
    """Load the four-column SGTree RF report as per-marker quality weights."""
    rf_path = Path(path)
    if not rf_path.exists():
        return {}

    with rf_path.open() as handle:
        header = next(handle, "").strip().split()
        expected = ["ProteinID", "MarkerGene", "RFdistance", "Status"]
        if header != expected:
            raise ValueError(
                f"malformed CMTV RF header in {rf_path}: expected {' '.join(expected)}"
            )
        marker_rf_values: dict[str, list[float]] = {}
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) != len(expected):
                raise ValueError(
                    f"malformed CMTV RF row {line_number} in {rf_path}"
                )
            marker = parts[1]
            try:
                rf_distance = float(parts[2])
            except ValueError as exc:
                raise ValueError(
                    f"invalid RFdistance on row {line_number} in {rf_path}"
                ) from exc
            marker_rf_values.setdefault(marker, []).append(rf_distance)

    if not marker_rf_values:
        raise ValueError(f"CMTV RF report contains no data rows: {rf_path}")
    marker_distances = {}
    for marker, values in marker_rf_values.items():
        if max(values) - min(values) > 1e-12:
            raise ValueError(
                f"inconsistent RFdistance values for marker {marker!r} in {rf_path}"
            )
        marker_distances[marker] = values[0]
    return cmtv_rf_quality_weights(marker_distances)


def _opaque(prefix: str, *parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode()
    return prefix + hashlib.blake2s(payload, digest_size=6).hexdigest()


def _balanced(items: list[str]) -> str | tuple:
    if len(items) == 1:
        return items[0]
    midpoint = len(items) // 2
    return _balanced(items[:midpoint]), _balanced(items[midpoint:])


def _pectinate(items: list[str]) -> str | tuple:
    if len(items) == 2:
        return items[0], items[1]
    return items[0], _pectinate(items[1:])


def _newick(node: str | tuple, branch_length: float) -> str:
    if isinstance(node, str):
        return f"{node}:{branch_length:.6f}"
    left, right = node
    return (
        f"({_newick(left, branch_length)},{_newick(right, branch_length)})"
        f"0.95:{branch_length:.6f}"
    )


def _perturb_branch_lengths(tree: Tree, panel_id: str, marker_name: str) -> None:
    """Add small frozen marker-specific variation without scorer feedback."""
    for node in tree.traverse():
        if node.is_root():
            continue
        descendants = ",".join(sorted(str(leaf.name) for leaf in node.iter_leaves()))
        payload = f"{panel_id}\x1f{marker_name}\x1f{descendants}".encode()
        value = int.from_bytes(
            hashlib.blake2s(payload, digest_size=4).digest(),
            byteorder="big",
        ) / (2**32 - 1)
        node.dist *= 0.97 + (0.06 * value)


def _regraft_target(
    tree: Tree,
    target_record: str,
    donor_record: str,
    replacement_record: str,
    *,
    support: float,
) -> Tree:
    """Move one leaf next to a donor while preserving other pair distances."""
    result = tree.copy(method="deepcopy")
    target = result & target_record
    old_parent = target.up
    if old_parent is None or old_parent.is_root() or len(old_parent.children) != 2:
        raise ValueError("fixture target must have a collapsible binary parent")
    sibling = next(child for child in old_parent.children if child is not target)
    grandparent = old_parent.up
    old_parent_distance = float(old_parent.dist)
    sibling_distance = float(sibling.dist)

    target.detach()
    sibling.detach()
    old_parent.detach()
    sibling.dist = sibling_distance + old_parent_distance
    grandparent.add_child(child=sibling)

    donor = result & donor_record
    donor_parent = donor.up
    if donor_parent is None:
        raise ValueError("fixture donor cannot be the root")
    donor_distance = float(donor.dist)
    donor.detach()
    attachment = donor_parent.add_child(
        name=f"{support:.2f}",
        dist=donor_distance / 2.0,
    )
    donor.dist = donor_distance / 2.0
    attachment.add_child(child=donor)
    target.name = replacement_record
    target.dist = donor_distance / 2.0
    attachment.add_child(child=target)
    return result


def _event(
    fixture_key: tuple[object, ...],
    kind: str,
    marker_name: str,
    native_record_id: str,
    observed_record_id: str,
    donor_genome: str,
    *,
    is_contaminant: bool,
    expected_loo: bool,
    expected_screen_candidate: bool,
) -> dict:
    contaminant_id = observed_record_id if is_contaminant else None
    return {
        "event_id": _opaque("e", *fixture_key, kind),
        "event_kind": kind,
        "marker_name": marker_name,
        "native_record_id": native_record_id,
        "contaminant_record_id": contaminant_id,
        "observed_record_id": observed_record_id,
        "donor_genome": donor_genome,
        "is_contaminant": is_contaminant,
        "expected_loo": expected_loo,
        "expected_screen_candidate": expected_screen_candidate,
    }


def _gene_votes(
    fixture_key: tuple[object, ...],
    event_kind: str,
    clades: list[str],
) -> list[dict]:
    return [
        {
            "gene_id": _opaque("q", *fixture_key, event_kind, index),
            "assigned_clade": clade,
            "informative": True,
        }
        for index, clade in enumerate(clades)
    ]


def _nearest_genome_set(tree: Tree, record_id: str, k: int) -> set[str]:
    target = tree & record_id
    ranked = sorted(
        (
            float(target.get_distance(other)),
            str(other.name).split("|", 1)[0],
        )
        for other in tree.iter_leaves()
        if other is not target
    )
    return {genome for _distance, genome in ranked[:k]}


def build_fixture(
    shape: str,
    event_class: str,
    seed: int,
    *,
    scale_taxa_count: int | None = None,
) -> dict:
    """Build one fixed tree panel; event choices never inspect model output."""
    if scale_taxa_count is None:
        if shape not in SHAPES or event_class not in EVENT_CLASSES or seed not in SEEDS:
            raise ValueError("unknown fixed fixture coordinate")
        fixture_tier = "mechanism"
        fixture_key = (shape, event_class, seed)
        panel_id = _opaque("p", shape, seed)
        taxa_count = SHAPES[shape]
    else:
        if (
            shape not in SHAPES
            or event_class not in SCALE_EVENT_CLASSES
            or seed != SCALE_SEED
            or scale_taxa_count not in SCALE_TAXA_COUNTS
        ):
            raise ValueError("unknown fixed scale fixture coordinate")
        fixture_tier = "scale"
        fixture_key = (
            fixture_tier,
            scale_taxa_count,
            shape,
            event_class,
            seed,
        )
        panel_id = _opaque("p", fixture_tier, scale_taxa_count, shape, seed)
        taxa_count = scale_taxa_count

    fixture_id = _opaque("f", *fixture_key)
    genomes = tuple(_opaque("x", panel_id, "genome", index) for index in range(taxa_count))
    markers = tuple(_opaque("m", panel_id, "marker", index) for index in range(MARKER_COUNT))
    seed_rank = SEEDS.index(seed)
    recipient_indices = (2 + seed_rank, 3 + seed_rank, 4 + seed_rank)
    donor_indices = (taxa_count - 1, taxa_count - 2, taxa_count - 3)

    records: dict[str, dict[str, str]] = {}
    trees: dict[str, Tree] = {}
    genome_topology = _balanced(list(genomes)) if shape == "balanced" else _pectinate(list(genomes))
    reference_tree = Tree(_newick(genome_topology, 1.0) + ";", format=1)
    for marker_index, marker_name in enumerate(markers):
        marker_records = {
            genome: build_sequence_id(
                genome,
                _opaque("c", panel_id, marker_name, genome),
                _opaque("q", panel_id, marker_name, genome),
            )
            for genome in genomes
        }
        records[marker_name] = marker_records
        ordered_records = [marker_records[genome] for genome in genomes]
        topology = _balanced(ordered_records) if shape == "balanced" else _pectinate(ordered_records)
        scale = 1.0 + (0.04 * marker_index) + (0.01 * seed_rank)
        trees[marker_name] = Tree(_newick(topology, scale) + ";", format=1)
        _perturb_branch_lengths(trees[marker_name], panel_id, marker_name)

    events: list[dict] = []
    votes_by_record: dict[str, list[dict]] = {}

    def add_regraft(
        *,
        marker_index: int,
        recipient_index: int,
        donor_index: int,
        kind: str,
        support: float,
        replacement: bool,
        expected_loo: bool,
        expected_screen_candidate: bool,
        vote_clades: list[str],
    ) -> None:
        marker_name = markers[marker_index]
        recipient = genomes[recipient_index]
        donor = genomes[donor_index]
        native_record = records[marker_name][recipient]
        donor_record = records[marker_name][donor]
        observed_record = (
            build_sequence_id(
                recipient,
                _opaque("c", fixture_id, kind, "contig"),
                _opaque("q", fixture_id, kind, "marker"),
            )
            if replacement
            else native_record
        )
        event = _event(
            fixture_key,
            kind,
            marker_name,
            native_record,
            observed_record,
            donor,
            is_contaminant=replacement,
            expected_loo=expected_loo,
            expected_screen_candidate=expected_screen_candidate,
        )
        neighbors_before = _nearest_genome_set(
            trees[marker_name],
            native_record,
            5,
        )
        trees[marker_name] = _regraft_target(
            trees[marker_name],
            native_record,
            donor_record,
            observed_record,
            support=support,
        )
        event["cmtv_k5_neighbors_preserved"] = (
            neighbors_before
            == _nearest_genome_set(trees[marker_name], observed_record, 5)
        )
        events.append(event)
        votes_by_record[observed_record] = _gene_votes(
            fixture_key,
            kind,
            vote_clades,
        )

    if event_class == "scale_far_source":
        donor = genomes[-1]
        add_regraft(
            marker_index=0,
            recipient_index=2,
            donor_index=taxa_count - 1,
            kind="scale_far_source",
            support=0.95,
            replacement=True,
            expected_loo=True,
            expected_screen_candidate=True,
            vote_clades=[donor] * 3,
        )
    elif event_class == "gene_rich_contaminant":
        donor = genomes[0]
        add_regraft(
            marker_index=0,
            recipient_index=2,
            donor_index=0,
            kind="replacement",
            support=0.95,
            replacement=True,
            expected_loo=True,
            expected_screen_candidate=True,
            vote_clades=[donor] * 3,
        )
    elif event_class == "solo_marker_contaminant":
        add_regraft(
            marker_index=0,
            recipient_index=2,
            donor_index=0,
            kind="solo",
            support=0.95,
            replacement=True,
            expected_loo=True,
            expected_screen_candidate=False,
            vote_clades=[],
        )
    elif event_class == "safety_controls":
        native_recipient = genomes[recipient_indices[0]]
        add_regraft(
            marker_index=0,
            recipient_index=recipient_indices[0],
            donor_index=donor_indices[0],
            kind="native",
            support=0.95,
            replacement=False,
            expected_loo=True,
            expected_screen_candidate=False,
            vote_clades=[native_recipient] * 3,
        )
        weak_donor = genomes[donor_indices[1]]
        add_regraft(
            marker_index=1,
            recipient_index=recipient_indices[1],
            donor_index=donor_indices[1],
            kind="weak_support",
            support=0.60,
            replacement=True,
            expected_loo=False,
            expected_screen_candidate=False,
            vote_clades=[weak_donor] * 3,
        )
        conflict_donor = genomes[donor_indices[2]]
        conflict_recipient = genomes[recipient_indices[2]]
        add_regraft(
            marker_index=2,
            recipient_index=recipient_indices[2],
            donor_index=donor_indices[2],
            kind="conflicting_votes",
            support=0.95,
            replacement=True,
            expected_loo=True,
            expected_screen_candidate=False,
            vote_clades=(
                [conflict_donor] * 4
                + [conflict_recipient] * 2
            ),
        )

    if event_class in {"gene_rich_contaminant", "solo_marker_contaminant"}:
        if not all(event["cmtv_k5_neighbors_preserved"] for event in events):
            raise AssertionError("fixed 2-to-0 regraft changed the CMTV k=5 neighbor set")

    return {
        "fixture_id": fixture_id,
        "fixture_tier": fixture_tier,
        "shape": shape,
        "event_class": event_class,
        "seed": seed,
        "genomes": genomes,
        "markers": markers,
        "reference_tree": reference_tree,
        "trees": trees,
        "events": events,
        "votes_by_record": votes_by_record,
    }


def build_fixtures() -> list[dict]:
    return [
        build_fixture(shape, event_class, seed)
        for shape in SHAPES
        for event_class in EVENT_CLASSES
        for seed in SEEDS
    ]


def build_scale_fixtures() -> list[dict]:
    return [
        build_fixture(
            shape,
            event_class,
            SCALE_SEED,
            scale_taxa_count=taxa_count,
        )
        for taxa_count in SCALE_TAXA_COUNTS
        for shape in SHAPES
        for event_class in SCALE_EVENT_CLASSES
    ]


def _score_cmtv_current(
    trees: dict[str, Tree],
    rows: list[dict],
    marker_rf_weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Run the tracked in-memory equivalent of CMTV's scoring core."""
    if marker_rf_weights is None:
        marker_rf_weights = {marker: 1.0 for marker in trees}
    missing_weights = sorted(set(trees) - set(marker_rf_weights))
    if missing_weights:
        raise ValueError(
            "CMTV voter weights are missing markers: "
            + ", ".join(missing_weights)
        )
    invalid_weights = sorted(
        marker
        for marker in trees
        if not np.isfinite(float(marker_rf_weights[marker]))
        or float(marker_rf_weights[marker]) <= 0.0
    )
    if invalid_weights:
        raise ValueError(
            "CMTV voter weights must be positive and finite: "
            + ", ".join(invalid_weights)
        )
    features = pd.DataFrame(
        {
            "genome": row["genome"],
            "marker_name": row["marker_name"],
            "leaf_name": row["leaf_name"],
        }
        for row in rows
    )
    marker_neighbors: dict[str, dict[str, set[str]]] = {}
    for marker_name in sorted(trees):
        tree = trees[marker_name]
        per_genome = {}
        for target in tree.iter_leaves():
            genome = str(target.name).split("|", 1)[0]
            ranked = sorted(
                (
                    float(target.get_distance(other)),
                    str(other.name).split("|", 1)[0],
                )
                for other in tree.iter_leaves()
                if str(other.name).split("|", 1)[0] != genome
            )
            per_genome[genome] = {other_genome for _distance, other_genome in ranked[:5]}
        marker_neighbors[marker_name] = per_genome

    agreement_rows = []
    genomes = sorted(
        {
            genome
            for per_genome in marker_neighbors.values()
            for genome in per_genome
        }
    )
    for genome in genomes:
        genome_marker_sets = {
            marker_name: per_genome[genome]
            for marker_name, per_genome in marker_neighbors.items()
            if genome in per_genome
        }
        for target_marker, target_set in genome_marker_sets.items():
            jaccards = []
            disagree_count = 0
            total_weight = 0.0
            for other_marker, other_set in genome_marker_sets.items():
                if other_marker == target_marker:
                    continue
                union = target_set | other_set
                jaccard = len(target_set & other_set) / len(union) if union else 1.0
                weight = float(marker_rf_weights[other_marker])
                jaccards.append(jaccard * weight)
                total_weight += weight
                if jaccard < 0.4:
                    disagree_count += 1
            voter_count = len(genome_marker_sets) - 1
            agreement_rows.append(
                {
                    "genome": genome,
                    "marker_name": target_marker,
                    "cmtv_agreement": (
                        sum(jaccards) / total_weight if total_weight else 1.0
                    ),
                    "cmtv_disagreement_frac": (
                        disagree_count / voter_count if voter_count else 0.0
                    ),
                    "cmtv_n_voters": voter_count,
                }
            )

    scored = features.merge(
        pd.DataFrame(agreement_rows),
        on=["genome", "marker_name"],
        how="left",
    )
    scored["cmtv_agreement"] = scored["cmtv_agreement"].fillna(1.0)
    scored["cmtv_disagreement_frac"] = scored["cmtv_disagreement_frac"].fillna(0.0)
    scored["cmtv_raw_score"] = 1.0 - scored["cmtv_agreement"]
    scored["cmtv_marker_z"] = scored.groupby("marker_name")["cmtv_raw_score"].transform(
        lambda values: (
            (values - values.mean()) / values.std()
            if values.std() > 0
            else 0.0
        )
    )
    scored["cmtv_genome_z"] = scored.groupby("genome")["cmtv_raw_score"].transform(
        lambda values: (
            (values - values.mean()) / values.std()
            if values.std() > 0
            else 0.0
        )
    )
    scored["cmtv_genome_pct"] = scored.groupby("genome")["cmtv_raw_score"].rank(
        pct=True
    )
    scored["cmtv_marker_z_pct"] = scored["cmtv_marker_z"].rank(pct=True)
    scored["cmtv_combined"] = np.sqrt(
        scored["cmtv_genome_pct"] * scored["cmtv_marker_z_pct"]
    )
    scored["cmtv_class"] = "clean"
    for _genome, group in scored.groupby("genome"):
        if len(group) <= 2:
            scored.loc[group.index, "cmtv_class"] = "ambiguous"
            continue
        top_index = group["cmtv_combined"].idxmax()
        top_disagree = scored.loc[top_index, "cmtv_disagreement_frac"]
        top_combined = scored.loc[top_index, "cmtv_combined"]
        top_marker_z = scored.loc[top_index, "cmtv_marker_z"]
        top_genome_z = scored.loc[top_index, "cmtv_genome_z"]
        if (
            top_disagree >= 0.6
            and top_combined >= 0.80
            and (top_genome_z >= 1.0 or top_marker_z >= 1.5)
        ):
            scored.loc[top_index, "cmtv_class"] = "contamination_candidate"
        elif top_disagree >= 0.4 and top_combined >= 0.65 and top_genome_z >= 0.5:
            scored.loc[top_index, "cmtv_class"] = "ambiguous"
    return scored


def score_fixture(fixture: dict) -> dict:
    loo_rows = score_loo_profiles(fixture["trees"])
    cmtv_rows = _score_cmtv_current(fixture["trees"], loo_rows)
    cmtv_candidates = cmtv_rows[
        cmtv_rows["cmtv_class"] == "contamination_candidate"
    ]
    cmtv_calls = set(
        zip(cmtv_candidates["marker_name"], cmtv_candidates["leaf_name"])
    )
    event_keys = [
        (event["marker_name"], event["observed_record_id"])
        for event in fixture["events"]
    ]
    if len(event_keys) != len(set(event_keys)):
        raise ValueError("duplicate marker and observed-record truth key")
    events_by_key = {
        key: event
        for key, event in zip(event_keys, fixture["events"])
    }

    loo_calls: set[tuple[str, str]] = set()
    screen_candidates: set[tuple[str, str]] = set()
    scored_by_key: dict[tuple[str, str], dict] = {}
    for row in loo_rows:
        record_id = row["leaf_name"]
        key = (row["marker_name"], record_id)
        if key in scored_by_key:
            raise ValueError("duplicate marker and scored-record key")
        gate = contig_gene_vote_gate(
            fixture["votes_by_record"].get(record_id, []),
            row["loo_attachment_clade"],
        )
        rf_before = None
        rf_after = None
        rf_guard_pass = False
        if row["loo_class"] == "discordant_marker":
            loo_calls.add(key)
            if gate["contig_gate_pass"]:
                original_tree = fixture["trees"][row["marker_name"]]
                remaining = [
                    str(leaf.name)
                    for leaf in original_tree.iter_leaves()
                    if str(leaf.name) != record_id
                ]
                candidate_tree = original_tree.copy(method="deepcopy")
                candidate_tree.prune(remaining)
                rf_before = _rf_distance_between(
                    fixture["reference_tree"],
                    original_tree,
                )
                rf_after = _rf_distance_between(
                    fixture["reference_tree"],
                    candidate_tree,
                )
                chosen_tree = choose_tree_by_rf(
                    fixture["reference_tree"],
                    original_tree,
                    candidate_tree,
                )
                rf_guard_pass = chosen_tree is candidate_tree
                if rf_guard_pass:
                    screen_candidates.add(key)
        scored_by_key[key] = {
            **row,
            **gate,
            "rf_before": rf_before,
            "rf_after": rf_after,
            "rf_guard_pass": rf_guard_pass,
        }

    event_results = []
    for key, event in sorted(events_by_key.items()):
        marker_name, record_id = key
        scored = scored_by_key[key]
        event_results.append(
            {
                "event_id": event["event_id"],
                "event_kind": event["event_kind"],
                "marker_name": marker_name,
                "record_id": record_id,
                "is_contaminant": event["is_contaminant"],
                "expected_loo": event["expected_loo"],
                "expected_screen_candidate": event["expected_screen_candidate"],
                "cmtv_k5_neighbors_preserved": event[
                    "cmtv_k5_neighbors_preserved"
                ],
                "loo_class": scored["loo_class"],
                "loo_abstention_reason": scored["loo_abstention_reason"],
                "contig_gate_pass": scored["contig_gate_pass"],
                "contig_abstention_reason": scored["contig_abstention_reason"],
                "loo_target_support": scored["loo_target_support"],
                "loo_score": scored["loo_score"],
                "loo_voter_mad": scored["loo_voter_mad"],
                "loo_voter_upper": scored["loo_voter_upper"],
                "loo_conflict_margin": scored["loo_conflict_margin"],
                "loo_marker_rank": scored["loo_marker_rank"],
                "loo_marker_margin": scored["loo_marker_margin"],
                "rf_before": scored["rf_before"],
                "rf_after": scored["rf_after"],
                "rf_guard_pass": scored["rf_guard_pass"],
                "screen_candidate": key in screen_candidates,
                "cmtv_called": key in cmtv_calls,
            }
        )

    truth_keys = {
        (event["marker_name"], event["observed_record_id"])
        for event in fixture["events"]
        if event["is_contaminant"]
    }
    expected_loo_keys = {
        (event["marker_name"], event["observed_record_id"])
        for event in fixture["events"]
        if event["expected_loo"]
    }
    expected_screen_keys = {
        (event["marker_name"], event["observed_record_id"])
        for event in fixture["events"]
        if event["expected_screen_candidate"]
    }
    case_pass = (
        loo_calls == expected_loo_keys
        and screen_candidates == expected_screen_keys
        and screen_candidates <= truth_keys
    )
    voter_uppers = [
        float(row["loo_voter_upper"])
        for row in loo_rows
        if row["loo_voter_upper"] is not None
    ]
    return {
        "fixture_id": fixture["fixture_id"],
        "fixture_tier": fixture["fixture_tier"],
        "shape": fixture["shape"],
        "event_class": fixture["event_class"],
        "seed": fixture["seed"],
        "taxa_count": len(fixture["genomes"]),
        "marker_count": len(fixture["markers"]),
        "truth_contaminant_count": len(truth_keys),
        "loo_call_count": len(loo_calls),
        "loo_true_positive_count": len(loo_calls & truth_keys),
        "loo_false_positive_count": len(loo_calls - truth_keys),
        "nonzero_voter_dispersion_count": sum(value > 0.0 for value in voter_uppers),
        "screen_candidate_count": len(screen_candidates),
        "screen_true_positive_count": len(screen_candidates & truth_keys),
        "screen_false_positive_count": len(screen_candidates - truth_keys),
        "cmtv_call_count": len(cmtv_calls),
        "cmtv_true_positive_count": len(cmtv_calls & truth_keys),
        "cmtv_false_positive_count": len(cmtv_calls - truth_keys),
        "case_pass": case_pass,
        "events": event_results,
    }


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> dict:
    if total == 0:
        return {"successes": 0, "total": 0, "rate": None, "low": None, "high": None}
    rate = successes / total
    denominator = 1.0 + (z * z / total)
    center = (rate + z * z / (2.0 * total)) / denominator
    half_width = (
        z
        * ((rate * (1.0 - rate) / total + z * z / (4.0 * total * total)) ** 0.5)
        / denominator
    )
    return {
        "successes": successes,
        "total": total,
        "rate": rate,
        "low": max(0.0, center - half_width),
        "high": min(1.0, center + half_width),
    }


def _case_rate(
    cases: list[dict],
    event_class: str,
    predicate,
    *,
    include_interval: bool = True,
) -> dict:
    selected = [case for case in cases if case["event_class"] == event_class]
    successes = sum(bool(predicate(case)) for case in selected)
    if include_interval:
        result = wilson_interval(successes, len(selected))
    else:
        result = {
            "successes": successes,
            "total": len(selected),
            "rate": successes / len(selected) if selected else None,
        }
    result["unit"] = "synthetic_fixture_case"
    return result


def run_benchmark() -> dict:
    started = time.monotonic()
    fixtures = build_fixtures()
    cases = [score_fixture(fixture) for fixture in fixtures]
    scale_fixtures = build_scale_fixtures()
    scale_cases = [score_fixture(fixture) for fixture in scale_fixtures]
    elapsed = time.monotonic() - started
    base_panel_successes = sum(
        all(
            case["case_pass"]
            for case in cases
            if case["shape"] == shape and case["seed"] == seed
        )
        for shape in SHAPES
        for seed in SEEDS
    )
    base_panel_interval = wilson_interval(
        base_panel_successes,
        len(SHAPES) * len(SEEDS),
    )
    base_panel_interval["unit"] = "paired_shape_seed_panel"
    serialized_trees = [
        tree.write(format=1)
        for fixture in fixtures
        for tree in fixture["trees"].values()
    ]
    screen_true_positives = sum(case["screen_true_positive_count"] for case in cases)
    screen_calls = sum(case["screen_candidate_count"] for case in cases)
    gene_rich_cases = [
        case for case in cases if case["event_class"] == "gene_rich_contaminant"
    ]
    gene_rich_true_positives = sum(
        case["screen_true_positive_count"] for case in gene_rich_cases
    )
    metrics = {
        "base_panel_expectations": base_panel_interval,
        "clean_screen_safety": _case_rate(
            cases,
            "clean",
            lambda case: case["screen_candidate_count"] == 0,
        ),
        "gene_rich_screen_detection": _case_rate(
            cases,
            "gene_rich_contaminant",
            lambda case: case["screen_true_positive_count"] == 1
            and case["screen_false_positive_count"] == 0,
        ),
        "solo_marker_screen_abstention": _case_rate(
            cases,
            "solo_marker_contaminant",
            lambda case: case["screen_candidate_count"] == 0,
        ),
        "safety_control_screen_safety": _case_rate(
            cases,
            "safety_controls",
            lambda case: case["screen_candidate_count"] == 0,
        ),
        "loo_gene_rich_detection": _case_rate(
            cases,
            "gene_rich_contaminant",
            lambda case: case["loo_true_positive_count"] == 1,
        ),
        "cmtv_gene_rich_detection": _case_rate(
            cases,
            "gene_rich_contaminant",
            lambda case: case["cmtv_true_positive_count"] == 1,
        ),
        "screen_marker_precision": {
            **wilson_interval(screen_true_positives, screen_calls),
            "unit": "hypothetical_marker_screen_candidate",
        },
        "screen_gene_rich_recall": {
            **wilson_interval(gene_rich_true_positives, len(gene_rich_cases)),
            "unit": "synthetic_gene_rich_fixture_case",
        },
    }
    scale_panel_successes = sum(
        all(
            case["case_pass"]
            for case in scale_cases
            if case["shape"] == shape and case["taxa_count"] == taxa_count
        )
        for taxa_count in SCALE_TAXA_COUNTS
        for shape in SHAPES
    )
    scale_panel_total = len(SCALE_TAXA_COUNTS) * len(SHAPES)
    scale_panel_rate = {
        "successes": scale_panel_successes,
        "total": scale_panel_total,
        "rate": scale_panel_successes / scale_panel_total,
        "unit": "paired_shape_taxa_panel",
    }
    scale_metrics = {
        "panel_expectations": scale_panel_rate,
        "clean_loo_safety": _case_rate(
            scale_cases,
            "scale_clean",
            lambda case: case["loo_call_count"] == 0,
            include_interval=False,
        ),
        "clean_screen_safety": _case_rate(
            scale_cases,
            "scale_clean",
            lambda case: case["screen_candidate_count"] == 0,
            include_interval=False,
        ),
        "clean_voter_dispersion": _case_rate(
            scale_cases,
            "scale_clean",
            lambda case: case["nonzero_voter_dispersion_count"] > 0,
            include_interval=False,
        ),
        "far_source_loo_detection": _case_rate(
            scale_cases,
            "scale_far_source",
            lambda case: case["loo_true_positive_count"] == 1
            and case["loo_false_positive_count"] == 0,
            include_interval=False,
        ),
        "far_source_screen_detection": _case_rate(
            scale_cases,
            "scale_far_source",
            lambda case: case["screen_true_positive_count"] == 1
            and case["screen_false_positive_count"] == 0,
            include_interval=False,
        ),
        "far_source_rf_improvement": _case_rate(
            scale_cases,
            "scale_far_source",
            lambda case: case["events"][0]["rf_guard_pass"]
            and case["events"][0]["rf_after"] < case["events"][0]["rf_before"],
            include_interval=False,
        ),
        "cmtv_far_source_detection": _case_rate(
            scale_cases,
            "scale_far_source",
            lambda case: case["cmtv_true_positive_count"] == 1,
            include_interval=False,
        ),
    }
    scale_serialized_trees = [
        tree.write(format=1)
        for fixture in scale_fixtures
        for tree in fixture["trees"].values()
    ]
    combined_serialized_trees = serialized_trees + scale_serialized_trees
    return {
        "benchmark_kind": "synthetic_tree_mechanism_screen",
        "production_pruning_enabled": False,
        "contig_vote_source": "preassigned_synthetic_fixture_clades",
        "comparison_scope": (
            "adversarial mechanism contrast conditioned on preserving each "
            "contaminant event's CMTV k=5 neighbor set"
        ),
        "limitations": [
            "contig votes are synthetic gate inputs, not sequence-derived assignments",
            "the CMTV contrast is not an unbiased model comparison",
            "the tree screen does not measure rebuilt species-tree RF or genome loss",
        ],
        "fixture_count": len(cases),
        "marker_tree_instance_count": len(serialized_trees),
        "distinct_marker_tree_count": len(set(serialized_trees)),
        "runtime_seconds": elapsed,
        "metrics": metrics,
        "cases": cases,
        "scale": {
            "comparison_scope": (
                "one-seed ideal-tree scale regression and safety check; "
                "CMTV detection is descriptive, not comparative superiority"
            ),
            "limitations": [
                "one fixed seed and idealized trees do not calibrate biological performance",
                "the scale tier uses synthetic contig votes and fixed far-source events",
                "near sources, missing taxa, low support, and inferred-tree error are not tested",
                "scale rates are descriptive counts; no confidence intervals are reported",
            ],
            "fixture_count": len(scale_cases),
            "marker_tree_instance_count": len(scale_serialized_trees),
            "distinct_marker_tree_count": len(set(scale_serialized_trees)),
            "taxa_counts": list(SCALE_TAXA_COUNTS),
            "metrics": scale_metrics,
            "cases": scale_cases,
        },
        "combined_fixture_count": len(cases) + len(scale_cases),
        "combined_marker_tree_instance_count": len(combined_serialized_trees),
        "combined_distinct_marker_tree_count": len(set(combined_serialized_trees)),
        "combined_panel_count": (
            len(SHAPES) * len(SEEDS)
            + len(SHAPES) * len(SCALE_TAXA_COUNTS)
        ),
        "combined_truth_event_count": sum(
            len(fixture["events"])
            for fixture in fixtures + scale_fixtures
        ),
    }


def check_benchmark(report: dict) -> list[str]:
    errors = []
    coordinates = {
        (case["shape"], case["event_class"], case["seed"])
        for case in report["cases"]
    }
    expected_coordinates = {
        (shape, event_class, seed)
        for shape in SHAPES
        for event_class in EVENT_CLASSES
        for seed in SEEDS
    }
    if coordinates != expected_coordinates or report["fixture_count"] != 24:
        errors.append("fixture matrix is not exactly 2 x 4 x 3")
    if any(case["marker_count"] != MARKER_COUNT for case in report["cases"]):
        errors.append("a fixture does not contain exactly eight marker trees")
    if {case["taxa_count"] for case in report["cases"]} != {10, 12}:
        errors.append("fixture taxa counts are not the fixed 10/12 design")
    if report["metrics"]["base_panel_expectations"]["successes"] != 6:
        errors.append("one or more paired shape-seed panels failed")
    for metric in (
        "clean_screen_safety",
        "gene_rich_screen_detection",
        "solo_marker_screen_abstention",
        "safety_control_screen_safety",
    ):
        if report["metrics"][metric]["rate"] != 1.0:
            errors.append(f"safety gate failed: {metric}")
    if (
        report["metrics"]["loo_gene_rich_detection"]["successes"]
        <= report["metrics"]["cmtv_gene_rich_detection"]["successes"]
    ):
        errors.append("LOO does not beat baseline CMTV in the fixed adversarial contrast")
    if report["metrics"]["screen_marker_precision"]["rate"] < 0.90:
        errors.append("hypothetical marker-screen precision is below 0.90")
    if report["metrics"]["screen_gene_rich_recall"]["rate"] < 0.75:
        errors.append("hypothetical gene-rich recall is below 0.75")
    if any(
        case["nonzero_voter_dispersion_count"] == 0
        for case in report["cases"]
        if case["event_class"] == "clean"
    ):
        errors.append("a clean fixture has zero voter dispersion throughout")
    scale = report["scale"]
    scale_coordinates = {
        (case["taxa_count"], case["shape"], case["event_class"], case["seed"])
        for case in scale["cases"]
    }
    expected_scale_coordinates = {
        (taxa_count, shape, event_class, SCALE_SEED)
        for taxa_count in SCALE_TAXA_COUNTS
        for shape in SHAPES
        for event_class in SCALE_EVENT_CLASSES
    }
    if scale_coordinates != expected_scale_coordinates or scale["fixture_count"] != 8:
        errors.append("scale fixture matrix is not exactly 2 x 2 x 2 x 1")
    if any(
        case["fixture_tier"] != "scale" or case["marker_count"] != MARKER_COUNT
        for case in scale["cases"]
    ):
        errors.append("a scale fixture has the wrong tier or marker count")
    if scale["taxa_counts"] != list(SCALE_TAXA_COUNTS):
        errors.append("scale fixture taxa counts are not the fixed 50/100 design")
    if scale["marker_tree_instance_count"] != 64:
        errors.append("scale fixture matrix does not contain exactly 64 marker trees")
    if scale["distinct_marker_tree_count"] != 36:
        errors.append("scale fixture matrix does not contain exactly 36 distinct trees")
    for metric in (
        "panel_expectations",
        "clean_loo_safety",
        "clean_screen_safety",
        "clean_voter_dispersion",
        "far_source_loo_detection",
        "far_source_screen_detection",
        "far_source_rf_improvement",
    ):
        if scale["metrics"][metric]["rate"] != 1.0:
            errors.append(f"scale gate failed: {metric}")
    if scale["metrics"]["cmtv_far_source_detection"]["total"] != 4:
        errors.append("scale CMTV detection result does not cover four far-source cases")
    combined_counts = (
        report["combined_fixture_count"],
        report["combined_marker_tree_instance_count"],
        report["combined_distinct_marker_tree_count"],
        report["combined_panel_count"],
        report["combined_truth_event_count"],
    )
    if combined_counts != (32, 256, 114, 10, 34):
        errors.append("combined scale and mechanism counts are not fixed at 32/256/114/10/34")
    if report["runtime_seconds"] > MAX_RUNTIME_SECONDS:
        errors.append(
            f"runtime {report['runtime_seconds']:.3f}s exceeds {MAX_RUNTIME_SECONDS:.0f}s"
        )
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail when a fixed screen gate fails")
    args = parser.parse_args(argv)
    report = run_benchmark()
    print(json.dumps(report, indent=2, sort_keys=True))
    if not args.check:
        return 0
    errors = check_benchmark(report)
    if errors:
        print("FAIL: " + "; ".join(errors))
        return 1
    print(
        f"PASS: {report['fixture_count']} mechanism fixtures + "
        f"{report['scale']['fixture_count']} scale fixtures, "
        f"{report['combined_marker_tree_instance_count']} marker-tree instances, "
        f"{report['combined_distinct_marker_tree_count']} distinct trees, "
        f"{report['runtime_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
