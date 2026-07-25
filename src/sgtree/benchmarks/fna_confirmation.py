"""Frozen empirical-FNA confirmation benchmark for US-010.

The benchmark keeps SGTree's LOO mode report-only.  It computes guarded
actions after inference so the benchmark can test the safety contract without
changing production pruning behavior.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tarfile
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree

from sgtree import __version__
from sgtree.benchmarks import (
    _extract_hmm_subset,
    _load_table,
    _read_normalized_proteomes,
)
from sgtree.id_schema import sanitize_token
from sgtree.benchmarks.loo_sequence_benchmark import _patristic_nrmse
from sgtree.benchmarks.gunc_low_memory import DIAMOND_BLOCK_SIZE
from sgtree.benchmarks.sequence_evidence import assign_contig_gene_split_votes
from sgtree.id_schema import parse_sequence_id
from sgtree.marker_selection import (
    _marker_name_from_tree_path,
    _rf_distance_between,
    build_singleton_output_tree,
    select_singleton_proposals,
)
from sgtree.marker_selection.contig_evidence import contig_gene_vote_gate
from sgtree.marker_selection.contig_evidence import (
    MIN_AGREEMENT_DENOMINATOR,
    MIN_AGREEMENT_NUMERATOR,
    MIN_CONFLICTING_GENES,
    MIN_INFORMATIVE_GENES,
)
from sgtree.marker_selection.loo_profile import (
    MIN_COORDINATES,
    MIN_MARKER_MARGIN,
    MIN_TARGET_DISCORDANCE,
    MIN_TARGET_SUPPORT,
    MIN_VOTERS,
)


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = PROJECT_ROOT / "runs/ml_contam_detection/us010_fna_confirmation"
DEFAULT_MODELS = PROJECT_ROOT / "resources/models/UNI56.hmm"
DEFAULT_GUNC_DB_DIR = PROJECT_ROOT / "resources/gunc_db"
SOURCE_BASE = PROJECT_ROOT / "runs/contig_variant_proof_of_concept/source_genomes"
STAGE_BASE = PROJECT_ROOT / "runs/contig_variant_proof_of_concept/benchmarks"
LINEAGES = ("alpha", "bacteroidota")
PANEL_SEEDS = (1103, 1301, 1601)
PANEL_GENOME_COUNT = 16
MARKER_COUNT = 8
GENE_RICH_NEIGHBOR_COUNT = 6
MAX_THREADS = 4
MAX_PEAK_RSS_KB = 8 * 1024 * 1024
MAX_MARKER_TREE_BUILDS = 120
EXPECTED_MARKER_TREE_BUILDS = 78
MANIFEST_SCHEMA = 1
BOOTSTRAP_REPLICATES = 10_000
CONTEXT_ORDER = ("clean", "near", "intermediate", "far", "solo")
EVENT_ORDER = ("near", "intermediate", "far", "solo", "sentinel")
DISTANCE_TARGET = {
    "near": 0.0,
    "intermediate": 0.5,
    "far": 1.0,
    "solo": 1.0,
    "sentinel": 1.0,
}


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_deterministic_tar_gz(
    destination: Path,
    members: list[tuple[Path, str]],
) -> None:
    def normalize(info: tarfile.TarInfo) -> tarfile.TarInfo:
        info.mtime = 0
        info.uid = 0
        info.gid = 0
        info.uname = ""
        info.gname = ""
        info.mode = 0o644
        info.pax_headers = {}
        return info

    with destination.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            mtime=0,
        ) as gzip_handle:
            with tarfile.open(
                fileobj=gzip_handle,
                mode="w",
                format=tarfile.PAX_FORMAT,
            ) as archive:
                for path, archive_name in members:
                    archive.add(
                        path,
                        arcname=archive_name,
                        recursive=False,
                        filter=normalize,
                    )


def _stable_key(*parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode()
    return hashlib.blake2s(payload, digest_size=16).hexdigest()


def _relative(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT))


def _project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _source_code_hash() -> str:
    paths = sorted((PROJECT_ROOT / "src/sgtree").rglob("*.py"))
    paths.extend(
        path
        for path in (PROJECT_ROOT / "pixi.toml", PROJECT_ROOT / "pixi.lock")
        if path.exists()
    )
    digest = hashlib.sha256()
    for path in paths:
        digest.update(_relative(path).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _git_metadata() -> dict[str, str | None]:
    def run(*args: str) -> str | None:
        result = subprocess.run(
            ["git", *args],
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        return result.stdout.strip() if result.returncode == 0 else None

    diff = run("diff", "--binary", "--", "src", "pixi.toml", "pixi.lock")
    return {
        "git_head": run("rev-parse", "HEAD"),
        "tracked_worktree_diff_sha256": (
            _sha256_bytes(diff.encode()) if diff is not None else None
        ),
        "source_code_sha256": _source_code_hash(),
    }


def _tool_versions() -> dict[str, str]:
    versions = {"sgtree": __version__}
    for module_name in ("Bio", "pandas", "pyhmmer", "pyrodigal"):
        module = __import__(module_name)
        versions[module_name] = str(getattr(module, "__version__", "unknown"))
    result = subprocess.run(
        ["VeryFastTree", "-version"],
        text=True,
        capture_output=True,
        check=False,
    )
    text = (result.stdout or result.stderr).strip().splitlines()
    versions["VeryFastTree"] = text[0] if text else "version unavailable"
    for executable, arguments in (
        ("gunc", ("--version",)),
        ("diamond", ("version",)),
        ("prodigal", ("-v",)),
    ):
        result = subprocess.run(
            [executable, *arguments],
            text=True,
            capture_output=True,
            check=False,
        )
        text = (result.stdout or result.stderr).strip().splitlines()
        versions[executable] = text[0] if text else "version unavailable"
    return versions


def _resolve_gunc_db(path: Path) -> Path:
    if path.is_file():
        return path.resolve()
    candidates = sorted(path.glob("*.dmnd")) if path.is_dir() else []
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected one GUNC DIAMOND database under {path}; found {len(candidates)}"
        )
    return candidates[0].resolve()


def _lineage_paths(lineage: str) -> tuple[Path, Path]:
    source_dir = SOURCE_BASE / lineage / "fna"
    stage_dir = STAGE_BASE / lineage / "_shared_staging"
    required = (
        stage_dir / "tree.nwk",
        stage_dir / "table_elim_dups",
        stage_dir / "gene_calls.tsv",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            f"Incomplete shared staging for {lineage}: {', '.join(missing)}"
        )
    if len(list(source_dir.glob("*.fna"))) < PANEL_GENOME_COUNT:
        raise FileNotFoundError(f"Too few source FNA files under {source_dir}")
    return source_dir, stage_dir


def _source_files(source_dir: Path) -> dict[str, Path]:
    files = {path.stem: path for path in sorted(source_dir.glob("*.fna"))}
    if len(files) != len(list(source_dir.glob("*.fna"))):
        raise ValueError(f"Duplicate FNA stems under {source_dir}")
    return files


def _load_lineage(lineage: str) -> dict:
    source_dir, stage_dir = _lineage_paths(lineage)
    table = _load_table(stage_dir / "table_elim_dups")
    gene_calls = pd.read_csv(stage_dir / "gene_calls.tsv", sep="\t")
    tree = Tree(str(stage_dir / "tree.nwk"), format=1)
    sources = _source_files(source_dir)
    tree_genomes = {str(leaf.name) for leaf in tree.iter_leaves()}
    shared = set(sources) & set(table["genome"].astype(str)) & tree_genomes
    if len(shared) < PANEL_GENOME_COUNT * len(PANEL_SEEDS):
        raise ValueError(
            f"{lineage} has only {len(shared)} genomes shared by FNA, hits, and source tree"
        )
    return {
        "lineage": lineage,
        "source_dir": source_dir,
        "stage_dir": stage_dir,
        "table": table[table["genome"].isin(shared)].copy(),
        "gene_calls": gene_calls[gene_calls["genome_id"].isin(shared)].copy(),
        "tree": tree,
        "source_files": {genome: sources[genome] for genome in sorted(shared)},
        "genomes": sorted(shared),
    }


def _select_markers(data: dict, genomes: list[str]) -> list[str]:
    table = data["table"][data["table"]["genome"].isin(genomes)]
    counts = table.groupby(["marker", "genome"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=genomes, fill_value=0)
    universal = counts.index[(counts == 1).all(axis=1)]
    if len(universal) < MARKER_COUNT:
        raise ValueError(
            f"{data['lineage']} panel has only {len(universal)} shared single-copy markers"
        )
    stats = (
        table[table["marker"].isin(universal)]
        .groupby("marker")["bitscore"]
        .agg(["median", "min"])
        .reset_index()
        .sort_values(
            ["median", "min", "marker"],
            ascending=[False, False, True],
        )
    )
    return stats.head(MARKER_COUNT)["marker"].astype(str).tolist()


def _select_panel_genomes(data: dict) -> list[tuple[int, list[str]]]:
    remaining = set(data["genomes"])
    counts = (
        data["table"]
        .groupby(["genome", "marker"])
        .size()
        .unstack(fill_value=0)
    )
    exact_single_copy = {
        genome: {
            str(marker)
            for marker, count in counts.loc[genome].items()
            if int(count) == 1
        }
        for genome in remaining
    }
    panels = []
    for seed in PANEL_SEEDS:
        selected = []
        shared_markers: set[str] | None = None
        while len(selected) < PANEL_GENOME_COUNT:
            candidates = []
            for genome in remaining - set(selected):
                candidate_shared = (
                    set(exact_single_copy[genome])
                    if shared_markers is None
                    else shared_markers & exact_single_copy[genome]
                )
                candidates.append(
                    (
                        -len(candidate_shared),
                        _stable_key("panel", data["lineage"], seed, genome),
                        genome,
                        candidate_shared,
                    )
                )
            if not candidates:
                break
            _negative_count, _tie_break, genome, candidate_shared = min(candidates)
            selected.append(genome)
            shared_markers = candidate_shared
        if len(selected) != PANEL_GENOME_COUNT:
            raise ValueError(f"Could not select 16 disjoint genomes for {data['lineage']}")
        if shared_markers is None or len(shared_markers) < MARKER_COUNT:
            raise ValueError(
                f"{data['lineage']} seed {seed} has only "
                f"{len(shared_markers or set())} shared single-copy markers"
            )
        panels.append((seed, sorted(selected)))
        remaining.difference_update(selected)
    return panels


def _marker_rows(data: dict, markers: list[str]) -> dict[tuple[str, str], dict]:
    selected = data["table"][data["table"]["marker"].isin(markers)]
    rows: dict[tuple[str, str], dict] = {}
    for (genome, marker), group in selected.groupby(["genome", "marker"]):
        if len(group) != 1:
            continue
        row = group.iloc[0].to_dict()
        row["normalized_header"] = str(row["savedname"]).replace("/", "|")
        rows[(str(genome), str(marker))] = row
    return rows


def _gene_rows(data: dict) -> dict[str, dict]:
    rows = {}
    for row in data["gene_calls"].to_dict("records"):
        header = str(row["normalized_header"])
        if header in rows:
            raise ValueError(f"Duplicate gene-call header: {header}")
        rows[header] = row
    return rows


def _stage_marker_proteins(data: dict) -> dict[str, str]:
    cached = data.get("stage_marker_proteins")
    if cached is not None:
        return cached
    marker_headers = {
        str(savedname).replace("/", "|")
        for savedname in data["table"]["savedname"]
    }
    proteins = {}
    with (data["stage_dir"] / "proteomes").open() as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in marker_headers:
                proteins[record.id] = str(record.seq)
    missing = marker_headers - set(proteins)
    if missing:
        raise ValueError(
            f"Shared staging is missing {len(missing)} marker protein records"
        )
    data["stage_marker_proteins"] = proteins
    return proteins


def _fragment_bounds(
    gene_calls: pd.DataFrame,
    target_header: str,
    selected_marker_headers: set[str],
) -> dict | None:
    matches = gene_calls[gene_calls["normalized_header"] == target_header]
    if len(matches) != 1:
        return None
    target = matches.iloc[0]
    contig = (
        gene_calls[
            (gene_calls["genome_id"] == target["genome_id"])
            & (gene_calls["contig_id"] == target["contig_id"])
        ]
        .sort_values(["begin", "end", "normalized_header"])
        .reset_index(drop=True)
    )
    target_positions = contig.index[contig["normalized_header"] == target_header].tolist()
    if len(target_positions) != 1:
        return None
    target_index = target_positions[0]

    left_barrier = -1
    for index in range(target_index - 1, -1, -1):
        if str(contig.iloc[index]["normalized_header"]) in selected_marker_headers:
            left_barrier = index
            break
    right_barrier = len(contig)
    for index in range(target_index + 1, len(contig)):
        if str(contig.iloc[index]["normalized_header"]) in selected_marker_headers:
            right_barrier = index
            break

    candidates = [
        index
        for index in range(left_barrier + 1, right_barrier)
        if index != target_index
    ]
    candidates.sort(
        key=lambda index: (
            abs(index - target_index),
            index,
        )
    )
    chosen_neighbors = sorted(candidates[:GENE_RICH_NEIGHBOR_COUNT])
    if len(chosen_neighbors) < GENE_RICH_NEIGHBOR_COUNT:
        return None
    selected_indices = sorted([target_index, *chosen_neighbors])
    first_index, last_index = selected_indices[0], selected_indices[-1]

    first = contig.iloc[first_index]
    last = contig.iloc[last_index]
    start0 = max(0, int(first["begin"]) - 1)
    end0 = int(last["end"])
    if first_index > 0:
        previous = contig.iloc[first_index - 1]
        start0 = max(start0, (int(previous["end"]) + int(first["begin"]) - 1) // 2)
    if last_index + 1 < len(contig):
        following = contig.iloc[last_index + 1]
        end0 = min(end0, (int(last["end"]) + int(following["begin"]) - 1) // 2)

    return {
        "source_contig_id": str(target["contig_id"]),
        "fragment_start0": start0,
        "fragment_end0": end0,
        "planned_non_marker_gene_count": len(chosen_neighbors),
        "planned_gene_headers": [
            str(contig.iloc[index]["normalized_header"])
            for index in selected_indices
        ],
    }


def _distance_ranks(
    tree: Tree,
    recipients: list[str],
    donors: list[str],
) -> dict[str, dict[str, tuple[float, int, float]]]:
    result = {}
    for recipient in recipients:
        ranked = sorted(
            (
                (float(tree.get_distance(recipient, donor)), donor)
                for donor in donors
                if donor != recipient
            ),
            key=lambda item: (item[0], item[1]),
        )
        denominator = max(1, len(ranked) - 1)
        result[recipient] = {
            donor: (distance, rank, rank / denominator)
            for rank, (distance, donor) in enumerate(ranked)
        }
    return result


def _read_assembly(path: Path) -> list[SeqRecord]:
    with path.open() as handle:
        records = [record for record in SeqIO.parse(handle, "fasta")]
    if not records:
        raise ValueError(f"Empty source assembly: {path}")
    return records


def _contig_by_normalized_id(records: list[SeqRecord], contig_id: str) -> SeqRecord:
    matches = [
        record
        for record in records
        if sanitize_token(record.id, record.id) == contig_id
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one source contig {contig_id}; found {len(matches)}")
    return matches[0]


def _coding_sequence(records: list[SeqRecord], gene_row: dict) -> Seq:
    contig = _contig_by_normalized_id(records, str(gene_row["contig_id"]))
    sequence = contig.seq[int(gene_row["begin"]) - 1 : int(gene_row["end"])]
    if int(gene_row["strand"]) < 0:
        sequence = sequence.reverse_complement()
    return sequence


def _source_assembly(data: dict, genome: str) -> list[SeqRecord]:
    cache = data.setdefault("source_assembly_cache", {})
    if genome not in cache:
        cache[genome] = _read_assembly(data["source_files"][genome])
    return cache[genome]


def _donor_gene_terminal_stop_codon(
    data: dict,
    gene_row: dict,
) -> str | None:
    cache = data.setdefault("donor_gene_terminal_stop_cache", {})
    header = str(gene_row["normalized_header"])
    if header not in cache:
        cds = _coding_sequence(
            _source_assembly(data, str(gene_row["genome_id"])),
            gene_row,
        )
        terminal_codon = str(cds[-3:]).upper()
        cache[header] = (
            terminal_codon
            if len(cds) >= 3
            and len(cds) % 3 == 0
            and terminal_codon in {"TAA", "TAG", "TGA"}
            else None
        )
    return cache[header]


def _donor_gene_has_terminal_stop(data: dict, gene_row: dict) -> bool:
    return _donor_gene_terminal_stop_codon(data, gene_row) is not None


def _donor_cds_matches_staged_protein(
    data: dict,
    gene_row: dict,
    staged_protein: str,
) -> bool:
    cache = data.setdefault("donor_cds_protein_cache", {})
    header = str(gene_row["normalized_header"])
    if header not in cache:
        cds = _coding_sequence(
            _source_assembly(data, str(gene_row["genome_id"])),
            gene_row,
        )
        cache[header] = _translated_protein(
            cds,
            int(gene_row["translation_table"]),
        )
    return str(cache[header]) == staged_protein


def _translated_protein(cds: Seq, translation_table: int) -> str:
    return str(cds.translate(table=translation_table)).rstrip("*")


def _event_candidate(
    *,
    data: dict,
    panel_id: str,
    role: str,
    panel_genomes: list[str],
    donor_pool: list[str],
    markers: list[str],
    marker_rows: dict[tuple[str, str], dict],
    gene_rows: dict[str, dict],
    distance_ranks: dict[str, dict[str, tuple[float, int, float]]],
    used_recipients: set[str],
    used_donors: set[str],
    used_markers: set[str],
    fragment_cache: dict[tuple[str, str], dict | None],
) -> dict:
    target = DISTANCE_TARGET[role]
    candidates = []
    marker_proteins = _stage_marker_proteins(data)
    for recipient in panel_genomes:
        if recipient in used_recipients:
            continue
        for marker in markers:
            if marker in used_markers:
                continue
            recipient_hit = marker_rows.get((recipient, marker))
            if recipient_hit is None:
                continue
            for donor in donor_pool:
                if donor in used_donors:
                    continue
                donor_hit = marker_rows.get((donor, marker))
                if donor_hit is None:
                    continue
                native_protein = marker_proteins[
                    str(recipient_hit["normalized_header"])
                ]
                donor_protein = marker_proteins[
                    str(donor_hit["normalized_header"])
                ]
                if native_protein == donor_protein:
                    continue
                donor_gene = gene_rows[str(donor_hit["normalized_header"])]
                if not _donor_gene_has_terminal_stop(data, donor_gene):
                    continue
                if (
                    role == "solo"
                    and _donor_gene_terminal_stop_codon(data, donor_gene)
                    not in {"TAA", "TAG"}
                ):
                    continue
                if role == "sentinel" and not _donor_cds_matches_staged_protein(
                    data,
                    donor_gene,
                    donor_protein,
                ):
                    continue
                fragment = None
                if role in {"near", "intermediate", "far"}:
                    cache_key = (donor, marker)
                    if cache_key not in fragment_cache:
                        donor_headers = {
                            str(marker_rows[(donor, selected_marker)]["normalized_header"])
                            for selected_marker in markers
                            if (donor, selected_marker) in marker_rows
                        }
                        fragment_cache[cache_key] = _fragment_bounds(
                            data["gene_calls"],
                            str(donor_hit["normalized_header"]),
                            donor_headers,
                        )
                    fragment = fragment_cache[cache_key]
                    if fragment is None:
                        continue
                distance, rank, quantile = distance_ranks[recipient][donor]
                candidates.append(
                    (
                        abs(quantile - target),
                        _stable_key(
                            "event",
                            panel_id,
                            role,
                            recipient,
                            donor,
                            marker,
                        ),
                        recipient,
                        donor,
                        marker,
                        distance,
                        rank,
                        quantile,
                        fragment,
                    )
                )
    if not candidates:
        raise ValueError(f"No feasible event candidate for {panel_id} {role}")
    (
        _distance_error,
        _tie_break,
        recipient,
        donor,
        marker,
        distance,
        rank,
        quantile,
        fragment,
    ) = min(candidates)

    recipient_hit = marker_rows[(recipient, marker)]
    donor_hit = marker_rows[(donor, marker)]
    recipient_gene = gene_rows[str(recipient_hit["normalized_header"])]
    donor_gene = gene_rows[str(donor_hit["normalized_header"])]
    native_protein = marker_proteins[
        str(recipient_hit["normalized_header"])
    ]
    donor_protein = marker_proteins[
        str(donor_hit["normalized_header"])
    ]
    if not donor_protein:
        raise ValueError(f"Empty translated donor target for {panel_id} {role}")

    event_id = f"{panel_id}_{role}"
    foreign_contig_id = sanitize_token(f"us010_{event_id}_foreign", event_id)
    expected_contig_id = (
        str(recipient_gene["contig_id"])
        if role == "sentinel"
        else foreign_contig_id
    )
    event = {
        "event_id": event_id,
        "panel_id": panel_id,
        "role": role,
        "event_class": (
            "gene_rich_replacement"
            if role in {"near", "intermediate", "far"}
            else "solo_marker_control"
            if role == "solo"
            else "native_contig_sentinel"
        ),
        "source_distance_stratum": (
            role if role in {"near", "intermediate", "far"} else "far_control"
        ),
        "recipient_genome": recipient,
        "donor_genome": donor,
        "marker": marker,
        "source_tree_distance": distance,
        "donor_distance_rank": rank,
        "donor_distance_quantile": quantile,
        "donor_pool_size": len(distance_ranks[recipient]),
        "distance_target_quantile": target,
        "distance_tie_break": "blake2s(panel,event,recipient,donor,marker)",
        "native_record_id": str(recipient_hit["normalized_header"]),
        "native_contig_id": str(recipient_gene["contig_id"]),
        "native_begin": int(recipient_gene["begin"]),
        "native_end": int(recipient_gene["end"]),
        "native_strand": int(recipient_gene["strand"]),
        "donor_record_id": str(donor_hit["normalized_header"]),
        "donor_contig_id": str(donor_gene["contig_id"]),
        "donor_begin": int(donor_gene["begin"]),
        "donor_end": int(donor_gene["end"]),
        "donor_strand": int(donor_gene["strand"]),
        "donor_translation_table": int(donor_gene["translation_table"]),
        "donor_cds_terminal_stop": True,
        "donor_terminal_stop_codon": _donor_gene_terminal_stop_codon(
            data,
            donor_gene,
        ),
        "native_protein_sha256": _sha256_bytes(native_protein.encode()),
        "donor_protein_sha256": _sha256_bytes(donor_protein.encode()),
        "donor_protein_differs_from_native": True,
        "expected_observed_contig_id": expected_contig_id,
        "foreign_contig_id": foreign_contig_id,
    }
    if fragment is not None:
        event.update(fragment)
    used_recipients.add(recipient)
    used_donors.add(donor)
    used_markers.add(marker)
    return event


def _panel_contexts(panel_id: str, events: list[dict]) -> list[dict]:
    event_by_role = {event["role"]: event["event_id"] for event in events}
    return [
        {
            "context_id": f"{panel_id}_{context}",
            "panel_id": panel_id,
            "context": context,
            "event_ids": (
                []
                if context == "clean"
                else [event_by_role[context], event_by_role["sentinel"]]
                if context in {"near", "intermediate", "far"}
                else [event_by_role["solo"]]
            ),
        }
        for context in CONTEXT_ORDER
    ]


def _freeze_panel(
    *,
    data: dict,
    panel_index: int,
    seed: int,
    panel_genomes: list[str],
    all_observed_genomes: set[str],
    markers: list[str],
    marker_rows: dict[tuple[str, str], dict],
    gene_rows: dict[str, dict],
) -> dict:
    panel_id = f"{data['lineage']}_p{panel_index}_seed{seed}"
    donor_pool = sorted(set(data["genomes"]) - all_observed_genomes)
    distance_ranks = _distance_ranks(
        data["tree"],
        panel_genomes,
        donor_pool,
    )
    used_recipients: set[str] = set()
    used_donors: set[str] = set()
    used_markers: set[str] = set()
    fragment_cache: dict[tuple[str, str], dict | None] = {}
    events = [
        _event_candidate(
            data=data,
            panel_id=panel_id,
            role=role,
            panel_genomes=panel_genomes,
            donor_pool=donor_pool,
            markers=markers,
            marker_rows=marker_rows,
            gene_rows=gene_rows,
            distance_ranks=distance_ranks,
            used_recipients=used_recipients,
            used_donors=used_donors,
            used_markers=used_markers,
            fragment_cache=fragment_cache,
        )
        for role in EVENT_ORDER
    ]
    return {
        "panel_id": panel_id,
        "lineage": data["lineage"],
        "panel_index": panel_index,
        "seed": seed,
        "genomes": panel_genomes,
        "reference_genomes": panel_genomes,
        "gunc_clean_control_genome": min(
            panel_genomes,
            key=lambda genome: _stable_key(
                "gunc_clean_control",
                panel_id,
                genome,
            ),
        ),
        "markers": markers,
        "donor_pool_excludes_all_observed_lineage_genomes": True,
        "donor_pool_size": len(donor_pool),
        "events": events,
        "contexts": _panel_contexts(panel_id, events),
    }


def _flatten_manifest(manifest: dict) -> tuple[list[dict], list[dict]]:
    events = []
    contexts = []
    for panel in manifest["panels"]:
        for event in panel["events"]:
            events.append(
                {
                    "lineage": panel["lineage"],
                    "seed": panel["seed"],
                    **event,
                }
            )
        for context in panel["contexts"]:
            contexts.append(
                {
                    "lineage": panel["lineage"],
                    "seed": panel["seed"],
                    **context,
                    "event_ids": ",".join(context["event_ids"]),
                }
            )
    return events, contexts


def _validate_frozen_manifest(manifest: dict) -> None:
    panels = manifest.get("panels", [])
    if len(panels) != 6:
        raise ValueError(f"Expected six panels; found {len(panels)}")
    contexts = [context for panel in panels for context in panel["contexts"]]
    events = [event for panel in panels for event in panel["events"]]
    if len(contexts) != 30:
        raise ValueError(f"Expected 30 contexts; found {len(contexts)}")
    counts = defaultdict(int)
    for event in events:
        counts[event["event_class"]] += 1
    expected = {
        "gene_rich_replacement": 18,
        "solo_marker_control": 6,
        "native_contig_sentinel": 6,
    }
    if dict(counts) != expected:
        raise ValueError(f"Unexpected event counts: {dict(counts)}")
    for panel in panels:
        if len(panel["genomes"]) != PANEL_GENOME_COUNT:
            raise ValueError(f"{panel['panel_id']} does not have 16 genomes")
        if len(panel["markers"]) != MARKER_COUNT:
            raise ValueError(f"{panel['panel_id']} does not have eight markers")
        if [row["context"] for row in panel["contexts"]] != list(CONTEXT_ORDER):
            raise ValueError(f"{panel['panel_id']} context order changed")
        if panel.get("gunc_clean_control_genome") not in panel["genomes"]:
            raise ValueError(f"{panel['panel_id']} has an invalid GUNC clean control")
        for field in ("recipient_genome", "donor_genome", "marker"):
            values = [event[field] for event in panel["events"]]
            if len(values) != len(set(values)):
                raise ValueError(f"{panel['panel_id']} reuses {field}")
        if set(event["donor_genome"] for event in panel["events"]) & set(panel["genomes"]):
            raise ValueError(f"{panel['panel_id']} uses an observed genome as donor")
        if not all(event.get("donor_cds_terminal_stop") is True for event in panel["events"]):
            raise ValueError(f"{panel['panel_id']} contains a partial donor marker")
        if not all(
            event.get("donor_terminal_stop_codon") in {"TAA", "TAG"}
            for event in panel["events"]
            if event["role"] == "solo"
        ):
            raise ValueError(
                f"{panel['panel_id']} contains a code-sensitive solo donor stop"
            )
        if not all(
            event.get("donor_protein_differs_from_native") is True
            and event.get("donor_protein_sha256")
            != event.get("native_protein_sha256")
            for event in panel["events"]
        ):
            raise ValueError(
                f"{panel['panel_id']} contains a sequence-identical replacement"
            )
    for lineage in LINEAGES:
        lineage_panels = [panel for panel in panels if panel["lineage"] == lineage]
        observed = [genome for panel in lineage_panels for genome in panel["genomes"]]
        if len(observed) != len(set(observed)):
            raise ValueError(f"{lineage} panels are not genome-disjoint")


def freeze_manifest(
    outdir: Path,
    models_path: Path,
    gunc_db_path: Path,
) -> Path:
    manifest_path = outdir / "frozen_manifest.json"
    if manifest_path.exists():
        raise FileExistsError(
            f"Frozen manifest already exists: {manifest_path}. "
            "Use the existing manifest rather than silently changing the confirmation set."
        )
    if not models_path.is_file():
        raise FileNotFoundError(models_path)
    gunc_db_path = _resolve_gunc_db(gunc_db_path)

    panels = []
    source_artifacts = []
    input_paths: set[Path] = set()
    for lineage in LINEAGES:
        data = _load_lineage(lineage)
        panel_specs = _select_panel_genomes(data)
        all_observed = {
            genome
            for _seed, genomes in panel_specs
            for genome in genomes
        }
        gene_rows = _gene_rows(data)
        lineage_panels = []
        for panel_index, (seed, genomes) in enumerate(panel_specs, start=1):
            markers = _select_markers(data, genomes)
            lineage_panels.append(
                _freeze_panel(
                    data=data,
                    panel_index=panel_index,
                    seed=seed,
                    panel_genomes=genomes,
                    all_observed_genomes=all_observed,
                    markers=markers,
                    marker_rows=_marker_rows(data, markers),
                    gene_rows=gene_rows,
                )
            )
        panels.extend(lineage_panels)
        selected_sources = {
            data["source_files"][genome]
            for panel in lineage_panels
            for genome in (
                list(panel["genomes"])
                + [event["donor_genome"] for event in panel["events"]]
            )
        }
        input_paths.update(selected_sources)
        for name in ("tree.nwk", "table_elim_dups", "gene_calls.tsv"):
            path = data["stage_dir"] / name
            source_artifacts.append(
                {
                    "lineage": lineage,
                    "role": name,
                    "path": _relative(path),
                    "sha256": _sha256_file(path),
                }
            )

    config = {
        "panel_seeds": list(PANEL_SEEDS),
        "panel_genome_count": PANEL_GENOME_COUNT,
        "marker_count": MARKER_COUNT,
        "gene_rich_neighbor_count": GENE_RICH_NEIGHBOR_COUNT,
        "solo_marker_terminal_stop_codons": ["TAA", "TAG"],
        "context_order": list(CONTEXT_ORDER),
        "event_order": list(EVENT_ORDER),
        "distance_targets": DISTANCE_TARGET,
        "event_selection_inputs": (
            "source-tree patristic distance, universal-marker hits, gene-call "
            "coordinates, and stable BLAKE2s tie-breaks only"
        ),
        "gunc_sampling": (
            "all 30 declared event assemblies plus one stable-hash clean "
            "control genome per base panel"
        ),
        "gunc_diamond_block_size": DIAMOND_BLOCK_SIZE,
        "scorer_output_used_for_event_selection": False,
        "loo_thresholds": {
            "min_voters": MIN_VOTERS,
            "min_coordinates": MIN_COORDINATES,
            "min_target_support": MIN_TARGET_SUPPORT,
            "min_target_discordance": MIN_TARGET_DISCORDANCE,
            "min_marker_margin": MIN_MARKER_MARGIN,
        },
        "contig_gate_thresholds": {
            "min_informative_genes": MIN_INFORMATIVE_GENES,
            "minimum_agreement": (
                f"{MIN_AGREEMENT_NUMERATOR}/{MIN_AGREEMENT_DENOMINATOR}"
            ),
            "min_conflicting_genes": MIN_CONFLICTING_GENES,
        },
        "resource_caps": {
            "threads": MAX_THREADS,
            "species_tree_threads": 1,
            "peak_rss_kb": MAX_PEAK_RSS_KB,
            "distinct_marker_tree_builds": MAX_MARKER_TREE_BUILDS,
            "target_distinct_marker_tree_builds": EXPECTED_MARKER_TREE_BUILDS,
        },
        "production_pruning": "disabled; loo_profile report-only",
    }
    source_fna_artifacts = [
        {
            "path": _relative(path),
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
        for path in sorted(input_paths)
    ]
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "story": "US-010",
        "wp1_selected_scorer": "loo",
        "wp1_decision": (
            "LOO and weighted CMTV tied under common gates; frozen simplicity "
            "order selected LOO and attributes the method claim to the shared contig gate."
        ),
        "prior_validation_seeds_excluded": [607, 809],
        "models_path": _relative(models_path),
        "models_sha256": _sha256_file(models_path),
        "models_bytes": models_path.stat().st_size,
        "gunc_database": {
            "name": "progenomes_2.1",
            "path": _relative(gunc_db_path),
            "sha256": _sha256_file(gunc_db_path),
            "bytes": gunc_db_path.stat().st_size,
        },
        "pixi_lock_sha256": _sha256_file(PROJECT_ROOT / "pixi.lock"),
        "code": _git_metadata(),
        "tool_versions": _tool_versions(),
        "config": config,
        "config_sha256": _sha256_bytes(
            json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
        ),
        "source_artifacts": source_artifacts,
        "source_fna_artifacts": source_fna_artifacts,
        "panels": panels,
    }
    _validate_frozen_manifest(manifest)
    manifest["manifest_content_sha256"] = _sha256_bytes(
        json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    )

    outdir.mkdir(parents=True, exist_ok=True)
    _write_json(manifest_path, manifest)
    events, contexts = _flatten_manifest(manifest)
    _write_tsv(
        outdir / "events.tsv",
        events,
        [
            "lineage",
            "panel_id",
            "seed",
            "event_id",
            "role",
            "event_class",
            "source_distance_stratum",
            "recipient_genome",
            "donor_genome",
            "marker",
            "source_tree_distance",
            "donor_distance_rank",
            "donor_distance_quantile",
            "donor_pool_size",
            "distance_target_quantile",
            "native_record_id",
            "native_contig_id",
            "donor_record_id",
            "donor_contig_id",
            "native_protein_sha256",
            "donor_protein_sha256",
            "donor_protein_differs_from_native",
            "donor_cds_terminal_stop",
            "donor_terminal_stop_codon",
            "expected_observed_contig_id",
            "planned_non_marker_gene_count",
            "fragment_start0",
            "fragment_end0",
        ],
    )
    _write_tsv(
        outdir / "contexts.tsv",
        contexts,
        [
            "lineage",
            "panel_id",
            "seed",
            "context_id",
            "context",
            "event_ids",
        ],
    )
    _write_tsv(
        outdir / "input_checksums.tsv",
        source_fna_artifacts,
        ["path", "sha256", "bytes"],
    )
    return manifest_path


def _load_manifest(outdir: Path) -> dict:
    path = outdir / "frozen_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(f"Freeze the confirmation manifest first: {path}")
    manifest = json.loads(path.read_text())
    _validate_frozen_manifest(manifest)
    expected = manifest.get("manifest_content_sha256")
    unsigned = dict(manifest)
    unsigned.pop("manifest_content_sha256", None)
    observed = _sha256_bytes(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":")).encode()
    )
    if expected != observed:
        raise ValueError(f"Frozen manifest content hash mismatch: {path}")
    frozen_source_hash = manifest.get("code", {}).get("source_code_sha256")
    current_source_hash = _source_code_hash()
    if frozen_source_hash != current_source_hash:
        raise RuntimeError(
            "Current source/environment hash differs from the frozen manifest: "
            f"frozen={frozen_source_hash}, current={current_source_hash}"
        )
    return manifest


def _event_by_id(manifest: dict) -> dict[str, dict]:
    return {
        event["event_id"]: event
        for panel in manifest["panels"]
        for event in panel["events"]
    }


def _verify_inference_inputs(manifest: dict) -> None:
    artifacts = [
        {
            "path": manifest["models_path"],
            "sha256": manifest["models_sha256"],
            "bytes": manifest["models_bytes"],
        },
        *manifest["source_fna_artifacts"],
    ]
    for artifact in artifacts:
        path = _project_path(str(artifact["path"]))
        if not path.is_file():
            raise FileNotFoundError(f"Frozen inference input is missing: {path}")
        if path.stat().st_size != int(artifact["bytes"]):
            raise ValueError(f"Frozen inference input size changed: {path}")
        if _sha256_file(path) != str(artifact["sha256"]):
            raise ValueError(f"Frozen inference input checksum changed: {path}")


def _replace_record_sequence(
    records: list[SeqRecord],
    contig_id: str,
    replacement: Seq,
) -> None:
    target = _contig_by_normalized_id(records, contig_id)
    target.seq = replacement


def _mask_native(records: list[SeqRecord], event: dict) -> None:
    contig = _contig_by_normalized_id(records, str(event["native_contig_id"]))
    start0 = int(event["native_begin"]) - 1
    end0 = int(event["native_end"])
    contig.seq = contig.seq[:start0] + Seq("N" * (end0 - start0)) + contig.seq[end0:]


def _donor_target(
    donor_records: list[SeqRecord],
    event: dict,
) -> tuple[Seq, str]:
    donor_gene = {
        "contig_id": event["donor_contig_id"],
        "begin": event["donor_begin"],
        "end": event["donor_end"],
        "strand": event["donor_strand"],
    }
    cds = _coding_sequence(donor_records, donor_gene)
    protein = _translated_protein(cds, int(event["donor_translation_table"]))
    if not protein:
        raise ValueError(f"Empty donor coding sequence for {event['event_id']}")
    return cds, protein


def _apply_event_to_records(
    recipient_records: list[SeqRecord],
    donor_records: list[SeqRecord],
    event: dict,
) -> None:
    donor_cds, _donor_protein = _donor_target(donor_records, event)
    role = event["role"]
    if role == "sentinel":
        contig = _contig_by_normalized_id(
            recipient_records,
            str(event["native_contig_id"]),
        )
        start0 = int(event["native_begin"]) - 1
        end0 = int(event["native_end"])
        genomic_cds = (
            donor_cds
            if int(event["native_strand"]) > 0
            else donor_cds.reverse_complement()
        )
        contig.seq = contig.seq[:start0] + genomic_cds + contig.seq[end0:]
        return

    _mask_native(recipient_records, event)
    donor_contig = _contig_by_normalized_id(
        donor_records,
        str(event["donor_contig_id"]),
    )
    if role in {"near", "intermediate", "far"}:
        fragment = donor_contig.seq[
            int(event["fragment_start0"]) : int(event["fragment_end0"])
        ]
    elif role == "solo":
        start0 = max(0, int(event["donor_begin"]) - 1 - 90)
        end0 = min(len(donor_contig.seq), int(event["donor_end"]) + 90)
        fragment = donor_contig.seq[start0:end0]
    else:
        raise ValueError(f"Unsupported event role: {role}")
    recipient_records.append(
        SeqRecord(
            fragment,
            id=str(event["foreign_contig_id"]),
            description="",
        )
    )


def _materialize_context(
    *,
    manifest: dict,
    panel: dict,
    context: dict,
    outdir: Path,
) -> Path:
    context_dir = outdir / "panels" / panel["panel_id"] / "contexts" / context["context"]
    inputs_dir = context_dir / "inputs"
    checksum_path = context_dir / "input_checksums.tsv"
    if checksum_path.is_file():
        rows = pd.read_csv(checksum_path, sep="\t")
        for row in rows.to_dict("records"):
            path = inputs_dir / str(row["file"])
            if not path.is_file() or _sha256_file(path) != str(row["sha256"]):
                raise ValueError(f"Existing materialized input failed checksum: {path}")
        return context_dir
    if inputs_dir.exists():
        raise FileExistsError(f"Partial context input directory requires inspection: {inputs_dir}")

    lineage = panel["lineage"]
    source_dir, _stage_dir = _lineage_paths(lineage)
    source_files = _source_files(source_dir)
    events = _event_by_id(manifest)
    context_events = [events[event_id] for event_id in context["event_ids"]]
    events_by_recipient: dict[str, list[dict]] = defaultdict(list)
    for event in context_events:
        events_by_recipient[event["recipient_genome"]].append(event)

    inputs_dir.mkdir(parents=True)
    for genome in panel["genomes"]:
        source_path = source_files[genome]
        output_path = inputs_dir / f"{genome}.fna"
        genome_events = events_by_recipient.get(genome, [])
        if not genome_events:
            shutil.copy2(source_path, output_path)
            continue
        records = _read_assembly(source_path)
        for event in genome_events:
            donor_records = _read_assembly(source_files[event["donor_genome"]])
            _apply_event_to_records(records, donor_records, event)
        with output_path.open("w") as handle:
            SeqIO.write(records, handle, "fasta")

    checksum_rows = [
        {
            "file": path.name,
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
        for path in sorted(inputs_dir.glob("*.fna"))
    ]
    _write_tsv(
        checksum_path,
        checksum_rows,
        ["file", "sha256", "bytes"],
    )
    _write_json(
        context_dir / "context_truth.json",
        {
            "panel_id": panel["panel_id"],
            "context": context["context"],
            "event_ids": context["event_ids"],
            "manifest_content_sha256": manifest["manifest_content_sha256"],
        },
    )
    return context_dir


def _resource_value(path: Path, label: str) -> float | None:
    if not path.is_file():
        return None
    prefix = f"{label}:"
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith(prefix):
            value = stripped.split(":", 1)[1].strip()
            try:
                return float(value)
            except ValueError:
                return None
    return None


def _run_context(
    *,
    manifest: dict,
    panel: dict,
    context: dict,
    outdir: Path,
    threads: int,
) -> dict:
    if not 1 <= threads <= MAX_THREADS:
        raise ValueError(f"threads must be between 1 and {MAX_THREADS}")
    context_dir = _materialize_context(
        manifest=manifest,
        panel=panel,
        context=context,
        outdir=outdir,
    )
    run_dir = context_dir / "sgtree"
    required = (
        run_dir / "tree_final.nwk",
        run_dir / "singleton_candidates.tsv",
        run_dir / "table_elim_dups",
        run_dir / "proteomes",
    )
    if all(path.is_file() for path in required):
        audit_path = context_dir / "truth_audit.json"
        if audit_path.is_file():
            audit = json.loads(audit_path.read_text())
            if audit.get("pass") is True:
                return audit
        audit = _truth_audit(manifest, panel, context, context_dir)
        _write_json(audit_path, audit)
        return audit
    if run_dir.exists():
        raise FileExistsError(f"Partial SGTree run requires inspection: {run_dir}")

    panel_dir = outdir / "panels" / panel["panel_id"]
    models_path = panel_dir / "markers.hmm"
    if not models_path.is_file():
        _extract_hmm_subset(
            _project_path(manifest["models_path"]),
            list(panel["markers"]),
            models_path,
        )
    cache_dir = outdir / "marker_tree_cache"
    resource_path = context_dir / "resource-usage.txt"
    log_path = context_dir / "sgtree.log"
    command = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(resource_path),
        sys.executable,
        "-m",
        "sgtree",
        str(context_dir / "inputs"),
        str(models_path),
        "--num_cpus",
        str(threads),
        "--percent_models",
        "0",
        "--save_dir",
        str(run_dir),
        "--marker_selection",
        "yes",
        "--singles",
        "yes",
        "--singles-mode",
        "loo_profile",
        "--selection_mode",
        "coordinate",
        "--selection_max_rounds",
        "5",
        "--selection_global_rounds",
        "1",
        "--keep_intermediates",
        "yes",
    ]
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONPATH": str(PROJECT_ROOT / "src"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "SGTREE_MARKER_TREE_CACHE_DIR": str(cache_dir),
            "SGTREE_SPECIES_TREE_THREADS": "1",
        }
    )
    started = time.monotonic()
    with log_path.open("w") as log_handle:
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            env=environment,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    wall_seconds = time.monotonic() - started
    if result.returncode != 0:
        raise RuntimeError(
            f"SGTree failed for {panel['panel_id']} {context['context']}; "
            f"see {log_path}"
        )
    cache_count = len(list(cache_dir.glob("*.nwk")))
    peak_rss_kb = _resource_value(
        resource_path,
        "Maximum resident set size (kbytes)",
    )
    if cache_count > MAX_MARKER_TREE_BUILDS:
        raise RuntimeError(
            f"Marker-tree cache contains {cache_count} builds; cap is "
            f"{MAX_MARKER_TREE_BUILDS}"
        )
    if peak_rss_kb is None or peak_rss_kb > MAX_PEAK_RSS_KB:
        raise RuntimeError(
            f"Peak RSS is {peak_rss_kb!r} kB; cap is {MAX_PEAK_RSS_KB} kB"
        )
    audit = _truth_audit(manifest, panel, context, context_dir)
    audit.update(
        {
            "wall_seconds": wall_seconds,
            "peak_rss_kb": int(peak_rss_kb),
            "marker_tree_cache_entries": cache_count,
            "command": command,
        }
    )
    _write_json(context_dir / "truth_audit.json", audit)
    return audit


def _truth_audit(
    manifest: dict,
    panel: dict,
    context: dict,
    context_dir: Path,
) -> dict:
    run_dir = context_dir / "sgtree"
    table = _load_table(run_dir / "table_elim_dups")
    selected = table[
        table["genome"].isin(panel["genomes"])
        & table["marker"].isin(panel["markers"])
    ]
    counts = selected.groupby(["genome", "marker"]).size()
    bad_cells = []
    for genome in panel["genomes"]:
        for marker in panel["markers"]:
            observed = int(counts.get((genome, marker), 0))
            if observed != 1:
                bad_cells.append(
                    {"genome": genome, "marker": marker, "observed_copy_count": observed}
                )

    proteomes = _read_normalized_proteomes(run_dir / "proteomes")
    selected_marker_record_ids = {
        str(row["savedname"]).replace("/", "|")
        for row in selected.to_dict("records")
    }
    event_lookup = _event_by_id(manifest)
    event_audits = []
    for event_id in context["event_ids"]:
        event = event_lookup[event_id]
        rows = selected[
            (selected["genome"] == event["recipient_genome"])
            & (selected["marker"] == event["marker"])
        ]
        observed_contigs = sorted(rows["contig_id"].astype(str).tolist())
        observed_hashes = []
        for row in rows.to_dict("records"):
            record_id = str(row["savedname"]).replace("/", "|")
            observed_hashes.append(
                _sha256_bytes(str(proteomes[event["recipient_genome"]][record_id].seq).encode())
            )
        same_contig_markers = sorted(
            selected[
                (selected["genome"] == event["recipient_genome"])
                & (
                    selected["contig_id"].astype(str)
                    == event["expected_observed_contig_id"]
                )
            ]["marker"]
            .astype(str)
            .tolist()
        )
        records_on_target_contig = {
            record_id
            for record_id in proteomes[event["recipient_genome"]]
            if parse_sequence_id(record_id)[1]
            == event["expected_observed_contig_id"]
        }
        non_marker_gene_count = len(
            records_on_target_contig - selected_marker_record_ids
        )
        contig_gene_layout_pass = (
            non_marker_gene_count >= MIN_INFORMATIVE_GENES
            if event["event_class"] == "gene_rich_replacement"
            else non_marker_gene_count < MIN_INFORMATIVE_GENES
            if event["event_class"] == "solo_marker_control"
            else True
        )
        target_contig_markers_pass = _target_contig_marker_layout_pass(
            event,
            same_contig_markers,
        )
        event_audits.append(
            {
                "event_id": event_id,
                "observed_marker_copy_count": len(rows),
                "observed_contigs": observed_contigs,
                "expected_contig": event["expected_observed_contig_id"],
                "observed_protein_sha256": observed_hashes,
                "expected_protein_sha256": event["donor_protein_sha256"],
                "selected_markers_on_target_contig": same_contig_markers,
                "observed_non_marker_gene_count": non_marker_gene_count,
                "contig_gene_layout_pass": contig_gene_layout_pass,
                "target_contig_marker_layout_pass": target_contig_markers_pass,
                "pass": (
                    len(rows) == 1
                    and observed_contigs == [event["expected_observed_contig_id"]]
                    and observed_hashes == [event["donor_protein_sha256"]]
                    and target_contig_markers_pass
                    and contig_gene_layout_pass
                ),
            }
        )

    clean_native_audits = []
    if context["context"] == "clean":
        for event in panel["events"]:
            rows = selected[
                (selected["genome"] == event["recipient_genome"])
                & (selected["marker"] == event["marker"])
            ]
            observed_contigs = sorted(rows["contig_id"].astype(str).tolist())
            observed_hashes = []
            for row in rows.to_dict("records"):
                record_id = str(row["savedname"]).replace("/", "|")
                observed_hashes.append(
                    _sha256_bytes(
                        str(
                            proteomes[event["recipient_genome"]][record_id].seq
                        ).encode()
                    )
                )
            clean_native_audits.append(
                {
                    "event_id": event["event_id"],
                    "observed_marker_copy_count": len(rows),
                    "observed_contigs": observed_contigs,
                    "expected_contig": event["native_contig_id"],
                    "observed_protein_sha256": observed_hashes,
                    "expected_protein_sha256": event[
                        "native_protein_sha256"
                    ],
                    "pass": (
                        len(rows) == 1
                        and observed_contigs == [event["native_contig_id"]]
                        and observed_hashes
                        == [event["native_protein_sha256"]]
                    ),
                }
            )

    final_tree = Tree(str(run_dir / "tree_final.nwk"), format=1)
    final_taxa = sorted(str(leaf.name) for leaf in final_tree.iter_leaves())
    report_only_tree_failures = []
    no_singles_paths = {
        _marker_name_from_tree_path(str(path)): path
        for path in (run_dir / "protTrees/no_singles").glob("*")
        if path.is_file()
    }
    for marker in panel["markers"]:
        path = no_singles_paths.get(marker)
        if path is None:
            report_only_tree_failures.append(
                {"marker": marker, "reason": "missing_no_singles_tree"}
            )
            continue
        tree = Tree(str(path), format=1)
        genomes = sorted(
            str(leaf.name).split("|", 1)[0]
            for leaf in tree.iter_leaves()
        )
        if genomes != sorted(panel["genomes"]):
            report_only_tree_failures.append(
                {
                    "marker": marker,
                    "reason": "report_only_tree_changed_genome_set",
                    "observed_genomes": genomes,
                }
            )
    audit = {
        "panel_id": panel["panel_id"],
        "context": context["context"],
        "selected_marker_cell_count": len(panel["genomes"]) * len(panel["markers"]),
        "bad_marker_cells": bad_cells,
        "event_audits": event_audits,
        "clean_native_audits": clean_native_audits,
        "final_taxa": final_taxa,
        "expected_final_taxa": sorted(panel["genomes"]),
        "report_only_marker_tree_failures": report_only_tree_failures,
        "pass": (
            not bad_cells
            and all(row["pass"] for row in event_audits)
            and all(row["pass"] for row in clean_native_audits)
            and final_taxa == sorted(panel["genomes"])
            and not report_only_tree_failures
        ),
    }
    if not audit["pass"]:
        _write_json(context_dir / "truth_audit.json", audit)
        raise RuntimeError(
            f"Exact truth audit failed for {panel['panel_id']} {context['context']}"
        )
    return audit


def _target_contig_marker_layout_pass(
    event: dict,
    observed_markers: list[str],
) -> bool:
    if event["role"] == "sentinel":
        return event["marker"] in observed_markers
    return observed_markers == [event["marker"]]


def _proteome_sequences(path: Path) -> dict[str, dict[str, str]]:
    records = _read_normalized_proteomes(path)
    return {
        genome: {
            record_id: str(record.seq)
            for record_id, record in genome_records.items()
        }
        for genome, genome_records in records.items()
    }


def _marker_tree_paths(run_dir: Path) -> dict[str, Path]:
    tree_dir = run_dir / "protTrees/no_duplicates/out"
    paths = {
        _marker_name_from_tree_path(str(path)): path
        for path in sorted(tree_dir.glob("*"))
        if path.is_file()
    }
    if not paths:
        raise FileNotFoundError(f"No no-duplicate marker trees under {tree_dir}")
    return paths


def _split_values(value: object) -> set[str]:
    if not isinstance(value, str):
        return set()
    return {item for item in value.split(",") if item}


def _evaluate_context(
    *,
    manifest: dict,
    panel: dict,
    context: dict,
    outdir: Path,
    clean_proteomes: dict[str, dict[str, str]],
    clean_tree: Tree,
) -> tuple[list[dict], dict]:
    context_dir = outdir / "panels" / panel["panel_id"] / "contexts" / context["context"]
    run_dir = context_dir / "sgtree"
    table = _load_table(run_dir / "table_elim_dups")
    selected = table[
        table["genome"].isin(panel["genomes"])
        & table["marker"].isin(panel["markers"])
    ].copy()
    candidates = pd.read_csv(
        run_dir / "singleton_candidates.tsv",
        sep="\t",
        keep_default_na=False,
    )
    candidate_lookup = {
        (str(row["marker_name"]), str(row["leaf_name"])): row
        for row in candidates.to_dict("records")
    }
    current_proteomes = _proteome_sequences(run_dir / "proteomes")
    marker_paths = _marker_tree_paths(run_dir)
    marker_trees = {
        marker: Tree(str(path), format=1)
        for marker, path in marker_paths.items()
    }
    marker_record_ids = {
        str(row["savedname"]).replace("/", "|")
        for row in selected.to_dict("records")
    }
    active_events = {
        (event["recipient_genome"], event["marker"]): event
        for event in (
            _event_by_id(manifest)[event_id]
            for event_id in context["event_ids"]
        )
    }

    raw_keys = {
        key
        for key, row in candidate_lookup.items()
        if str(row.get("loo_class", "")) == "discordant_marker"
    }
    gate_audits: dict[tuple[str, str], dict] = {}
    gated_proposals = []
    for marker, leaf_name in sorted(raw_keys):
        row = candidate_lookup[(marker, leaf_name)]
        genome, contig_id, _gene_id = parse_sequence_id(leaf_name)
        candidate_genes = {
            record_id: sequence
            for record_id, sequence in current_proteomes.get(genome, {}).items()
            if parse_sequence_id(record_id)[1] == contig_id
        }
        attachment = _split_values(row.get("loo_attachment_taxa"))
        marker_tree = marker_trees[marker]
        background = {
            str(leaf.name).split("|", 1)[0]
            for leaf in marker_tree.iter_leaves()
        } - {genome} - attachment
        evidence = assign_contig_gene_split_votes(
            candidate_genes,
            clean_proteomes,
            recipient_genome=genome,
            candidate_contig_id=contig_id,
            marker_record_ids=marker_record_ids,
            attachment_taxa=attachment,
            background_taxa=background,
        )
        gate = contig_gene_vote_gate(
            evidence["votes"],
            str(row.get("loo_attachment_clade", "")) or None,
        )
        gate_audits[(marker, leaf_name)] = {
            "evidence": evidence,
            "gate": gate,
        }
        if gate["contig_gate_pass"]:
            gated_proposals.append(
                {
                    "marker_name": marker,
                    "leaf_name": leaf_name,
                    "genome": genome,
                    "score": float(row.get("loo_score") or 0.0),
                    "delta_rf": 0.0,
                    "topoknn_score": 0.0,
                }
            )

    budgeted = select_singleton_proposals(
        gated_proposals,
        genome_marker_counts={genome: MARKER_COUNT for genome in panel["genomes"]},
        min_markers_per_genome=1,
        max_prunes_per_genome=1,
    )
    budget_keys = {
        (str(proposal["marker_name"]), str(proposal["leaf_name"]))
        for proposal in budgeted
    }
    proposals_by_marker: dict[str, list[str]] = defaultdict(list)
    for marker, leaf_name in budget_keys:
        proposals_by_marker[marker].append(leaf_name)

    action_keys: set[tuple[str, str]] = set()
    rf_decisions: dict[tuple[str, str], dict] = {}
    species_tree_path = run_dir / "tree_final.nwk"
    species_tree = Tree(str(species_tree_path), format=1)
    for marker, leaf_names in sorted(proposals_by_marker.items()):
        original = marker_trees[marker]
        chosen, decision = build_singleton_output_tree(
            marker_tree_path=str(marker_paths[marker]),
            species_tree_path=str(species_tree_path),
            accepted_leaf_names=sorted(leaf_names),
            mode="loo_profile",
        )
        chosen_leaves = {str(leaf.name) for leaf in chosen.iter_leaves()}
        remaining = [
            str(leaf.name)
            for leaf in original.iter_leaves()
            if str(leaf.name) not in set(leaf_names)
        ]
        candidate = original.copy(method="deepcopy")
        candidate.prune(remaining)
        rf_before = _rf_distance_between(species_tree, original)
        rf_after = _rf_distance_between(species_tree, candidate)
        for leaf_name in leaf_names:
            key = (marker, leaf_name)
            removed = leaf_name not in chosen_leaves
            if removed:
                action_keys.add(key)
            rf_decisions[key] = {
                "rf_guard_decision": decision,
                "rf_before": rf_before,
                "rf_after": rf_after,
                "guarded_action": removed,
            }

    rows = []
    for hit in selected.sort_values(["genome", "marker"]).to_dict("records"):
        genome = str(hit["genome"])
        marker = str(hit["marker"])
        leaf_name = str(hit["savedname"]).replace("/", "|")
        key = (marker, leaf_name)
        candidate = candidate_lookup.get(key, {})
        event = active_events.get((genome, marker))
        event_class = event["event_class"] if event else "unaffected_cell"
        truth_positive = event_class == "gene_rich_replacement"
        gate_audit = gate_audits.get(key)
        gate = gate_audit["gate"] if gate_audit else {}
        evidence = gate_audit["evidence"] if gate_audit else {}
        raw_call = key in raw_keys
        gate_pass = bool(gate.get("contig_gate_pass", False))
        budget_pass = key in budget_keys
        guarded_action = key in action_keys
        if not raw_call:
            final_reason = "no_raw_call"
        elif not gate_pass:
            final_reason = "contig_gate:" + str(
                gate.get("contig_abstention_reason", "not_evaluated")
            )
        elif not budget_pass:
            final_reason = "budget_rejected"
        elif not guarded_action:
            final_reason = "rf_guard_kept"
        else:
            final_reason = "guarded_action"
        rows.append(
            {
                "panel_id": panel["panel_id"],
                "lineage": panel["lineage"],
                "seed": panel["seed"],
                "context": context["context"],
                "event_id": event["event_id"] if event else "",
                "event_class": event_class,
                "source_distance_stratum": (
                    event["source_distance_stratum"] if event else ""
                ),
                "truth_positive": truth_positive,
                "recipient_genome": genome,
                "marker": marker,
                "leaf_name": leaf_name,
                "contig_id": str(hit["contig_id"]),
                "raw_call": raw_call,
                "loo_class": str(candidate.get("loo_class", "")),
                "loo_abstention_reason": str(
                    candidate.get("loo_abstention_reason", "")
                ),
                "loo_score": candidate.get("loo_score", ""),
                "loo_target_support": candidate.get("loo_target_support", ""),
                "loo_attachment_clade": candidate.get("loo_attachment_clade", ""),
                "contig_input_status": evidence.get("input_status", ""),
                "informative_gene_count": gate.get("informative_gene_count", ""),
                "contig_gate_pass": gate_pass,
                "contig_gate_reason": gate.get("contig_abstention_reason", ""),
                "budget_pass": budget_pass,
                "rf_guard_decision": rf_decisions.get(key, {}).get(
                    "rf_guard_decision", ""
                ),
                "rf_before": rf_decisions.get(key, {}).get("rf_before", ""),
                "rf_after": rf_decisions.get(key, {}).get("rf_after", ""),
                "guarded_action": guarded_action,
                "final_reason": final_reason,
                "false_action": guarded_action and not truth_positive,
            }
        )

    contaminated_tree = Tree(str(run_dir / "tree_final.nwk"), format=1)
    rf, maxrf, *_ = clean_tree.robinson_foulds(
        contaminated_tree,
        unrooted_trees=True,
    )
    summary = {
        "panel_id": panel["panel_id"],
        "lineage": panel["lineage"],
        "seed": panel["seed"],
        "context": context["context"],
        "truth_positive_count": sum(row["truth_positive"] for row in rows),
        "raw_true_positive_count": sum(
            row["truth_positive"] and row["raw_call"] for row in rows
        ),
        "guarded_true_positive_count": sum(
            row["truth_positive"] and row["guarded_action"] for row in rows
        ),
        "guarded_action_count": sum(row["guarded_action"] for row in rows),
        "false_action_count": sum(row["false_action"] for row in rows),
        "abstention_count": sum(
            bool(row["loo_abstention_reason"]) for row in rows
        ),
        "normalized_rf_vs_clean": rf / maxrf if maxrf else 0.0,
        "normalized_patristic_error_vs_clean": _patristic_nrmse(
            outdir
            / "panels"
            / panel["panel_id"]
            / "contexts/clean/sgtree/tree_final.nwk",
            run_dir / "tree_final.nwk",
        ),
    }
    return rows, summary


PER_CELL_FIELDS = [
    "panel_id",
    "lineage",
    "seed",
    "context",
    "event_id",
    "event_class",
    "source_distance_stratum",
    "truth_positive",
    "recipient_genome",
    "marker",
    "leaf_name",
    "contig_id",
    "raw_call",
    "loo_class",
    "loo_abstention_reason",
    "loo_score",
    "loo_target_support",
    "loo_attachment_clade",
    "contig_input_status",
    "informative_gene_count",
    "contig_gate_pass",
    "contig_gate_reason",
    "budget_pass",
    "rf_guard_decision",
    "rf_before",
    "rf_after",
    "guarded_action",
    "final_reason",
    "false_action",
]


def _evaluate_panel(manifest: dict, panel: dict, outdir: Path) -> dict:
    panel_dir = outdir / "panels" / panel["panel_id"]
    clean_run = panel_dir / "contexts/clean/sgtree"
    clean_proteomes = _proteome_sequences(clean_run / "proteomes")
    clean_tree = Tree(str(clean_run / "tree_final.nwk"), format=1)
    rows = []
    summaries = []
    for context in panel["contexts"]:
        context_rows, context_summary = _evaluate_context(
            manifest=manifest,
            panel=panel,
            context=context,
            outdir=outdir,
            clean_proteomes=clean_proteomes,
            clean_tree=clean_tree,
        )
        rows.extend(context_rows)
        summaries.append(context_summary)
    _write_tsv(panel_dir / "per_cell_results.tsv", rows, PER_CELL_FIELDS)
    _write_tsv(
        panel_dir / "context_summary.tsv",
        summaries,
        list(summaries[0]),
    )
    false_actions = sum(row["false_action"] for row in rows)
    audit = {
        "panel_id": panel["panel_id"],
        "context_count": len(summaries),
        "cell_count": len(rows),
        "gene_rich_truth_count": sum(row["truth_positive"] for row in rows),
        "gene_rich_raw_call_count": sum(
            row["truth_positive"] and row["raw_call"] for row in rows
        ),
        "gene_rich_guarded_action_count": sum(
            row["truth_positive"] and row["guarded_action"] for row in rows
        ),
        "false_action_count": false_actions,
        "control_guarded_action_count": sum(
            row["guarded_action"]
            and row["event_class"]
            in {"solo_marker_control", "native_contig_sentinel"}
            for row in rows
        ),
        "pass": false_actions == 0,
    }
    _write_json(panel_dir / "panel_audit.json", audit)
    return audit


def run_confirmation(
    outdir: Path,
    *,
    threads: int,
    pilot_only: bool,
) -> dict:
    manifest = _load_manifest(outdir)
    _verify_inference_inputs(manifest)
    panels = manifest["panels"]
    selected_panels = panels[:1] if pilot_only else panels
    if not pilot_only:
        pilot_path = outdir / "pilot_status.json"
        if not pilot_path.is_file():
            raise RuntimeError("Run and pass the one-panel pilot before the full benchmark")
        pilot_status = json.loads(pilot_path.read_text())
        if pilot_status.get("pass") is not True:
            raise RuntimeError("The saved one-panel pilot did not pass")

    audits = []
    for panel in selected_panels:
        context_audits = [
            _run_context(
                manifest=manifest,
                panel=panel,
                context=context,
                outdir=outdir,
                threads=threads,
            )
            for context in panel["contexts"]
        ]
        panel_audit = _evaluate_panel(manifest, panel, outdir)
        if not all(audit.get("pass") is True for audit in context_audits):
            raise RuntimeError(f"Truth audit failed in {panel['panel_id']}")
        if not panel_audit["pass"]:
            raise RuntimeError(f"Safety audit failed in {panel['panel_id']}")
        audits.append(
            {
                **panel_audit,
                "peak_rss_kb": max(
                    int(audit.get("peak_rss_kb", 0))
                    for audit in context_audits
                ),
                "wall_seconds": sum(
                    float(audit.get("wall_seconds", 0.0))
                    for audit in context_audits
                ),
                "marker_tree_cache_entries": max(
                    int(audit.get("marker_tree_cache_entries", 0))
                    for audit in context_audits
                ),
            }
        )
    status = {
        "mode": "pilot" if pilot_only else "full",
        "panels_completed": len(selected_panels),
        "contexts_completed": len(selected_panels) * len(CONTEXT_ORDER),
        "panel_audits": audits,
        "marker_tree_cache_entries": len(
            list((outdir / "marker_tree_cache").glob("*.nwk"))
        ),
        "peak_rss_kb": max(row["peak_rss_kb"] for row in audits),
        "wall_seconds": sum(row["wall_seconds"] for row in audits),
        "production_pruning": "disabled",
        "pass": all(row["pass"] for row in audits),
    }
    _write_json(
        outdir / ("pilot_status.json" if pilot_only else "run_status.json"),
        status,
    )
    return status


def _gunc_summary_paths(gunc_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in (gunc_dir / "output").rglob("*maxCSS_level.tsv")
        if path.is_file()
    )


def run_gunc(outdir: Path, *, threads: int) -> dict:
    if not 1 <= threads <= MAX_THREADS:
        raise ValueError(f"threads must be between 1 and {MAX_THREADS}")
    manifest = _load_manifest(outdir)
    full_status = check_confirmation(outdir, require_full=True)
    if full_status["pass"] is not True:
        raise RuntimeError("The full SGTree confirmation must pass before GUNC")
    database = _project_path(manifest["gunc_database"]["path"])
    if not database.is_file():
        raise FileNotFoundError(f"Frozen GUNC database is missing: {database}")
    if database.stat().st_size != int(manifest["gunc_database"]["bytes"]):
        raise ValueError("Frozen GUNC database size changed")
    database_sha256 = _sha256_file(database)
    if database_sha256 != manifest["gunc_database"]["sha256"]:
        raise ValueError("Frozen GUNC database checksum changed")

    gunc_dir = outdir / "gunc"
    unique_dir = gunc_dir / "unique_inputs"
    unique_dir.mkdir(parents=True, exist_ok=True)
    event_lookup = _event_by_id(manifest)
    mapping_rows = []
    unique_sources: dict[str, Path] = {}
    for panel in manifest["panels"]:
        for context in panel["contexts"]:
            context_dir = (
                outdir
                / "panels"
                / panel["panel_id"]
                / "contexts"
                / context["context"]
            )
            checksum_rows = pd.read_csv(
                context_dir / "input_checksums.tsv",
                sep="\t",
            )
            events = [event_lookup[event_id] for event_id in context["event_ids"]]
            selected_genomes = (
                {panel["gunc_clean_control_genome"]}
                if context["context"] == "clean"
                else {event["recipient_genome"] for event in events}
            )
            for row in checksum_rows.to_dict("records"):
                source = context_dir / "inputs" / str(row["file"])
                genome = Path(str(row["file"])).stem
                if genome not in selected_genomes:
                    continue
                checksum = str(row["sha256"])
                unique_sources.setdefault(checksum, source)
                mapping_rows.append(
                    {
                        "panel_id": panel["panel_id"],
                        "lineage": panel["lineage"],
                        "seed": panel["seed"],
                        "context": context["context"],
                        "genome": genome,
                        "assembly_sha256": checksum,
                        "gunc_input_id": checksum,
                        "event_ids": ",".join(
                            event["event_id"]
                            for event in events
                            if event["recipient_genome"] == genome
                        ),
                    }
                )
    for checksum, source in sorted(unique_sources.items()):
        target = unique_dir / f"{checksum}.fna"
        if target.exists():
            if _sha256_file(target) != checksum:
                raise ValueError(f"Corrupt cached GUNC input: {target}")
            continue
        try:
            os.link(source, target)
        except OSError:
            shutil.copy2(source, target)
    expected_mapping_count = sum(
        1
        if context["context"] == "clean"
        else len(context["event_ids"])
        for panel in manifest["panels"]
        for context in panel["contexts"]
    )
    expected_unique_count = sum(
        1 + len(panel["events"])
        for panel in manifest["panels"]
    )
    if len(mapping_rows) != expected_mapping_count:
        raise ValueError(
            f"Expected {expected_mapping_count} GUNC assembly mappings; "
            f"found {len(mapping_rows)}"
        )
    if len(unique_sources) != expected_unique_count:
        raise ValueError(
            f"Expected {expected_unique_count} unique GUNC assemblies; "
            f"found {len(unique_sources)}"
        )
    expected_input_names = {
        f"{checksum}.fna" for checksum in unique_sources
    }
    observed_input_names = {
        path.name for path in unique_dir.glob("*.fna") if path.is_file()
    }
    if observed_input_names != expected_input_names:
        raise ValueError("GUNC unique-input directory contains stale or missing files")
    _write_tsv(
        gunc_dir / "assembly_map.tsv",
        mapping_rows,
        [
            "panel_id",
            "lineage",
            "seed",
            "context",
            "genome",
            "assembly_sha256",
            "gunc_input_id",
            "event_ids",
        ],
    )
    _write_tsv(
        gunc_dir / "unique_inputs.tsv",
        [
            {
                "gunc_input_id": checksum,
                "source_path": _relative(source),
                "sha256": checksum,
                "bytes": source.stat().st_size,
            }
            for checksum, source in sorted(unique_sources.items())
        ],
        ["gunc_input_id", "source_path", "sha256", "bytes"],
    )

    output_dir = gunc_dir / "output"
    summary_paths = _gunc_summary_paths(gunc_dir)
    status_path = gunc_dir / "status.json"
    if summary_paths and status_path.is_file():
        status = json.loads(status_path.read_text())
        if (
            status.get("pass") is True
            and status.get("database_sha256") == database_sha256
        ):
            return status
    if output_dir.exists():
        raise FileExistsError(f"Partial GUNC output requires inspection: {output_dir}")
    output_dir.mkdir(parents=True)
    temp_dir = gunc_dir / "tmp"
    temp_dir.mkdir()
    resource_path = gunc_dir / "resource-usage.txt"
    log_path = gunc_dir / "gunc.log"
    command = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(resource_path),
        sys.executable,
        "-m",
        "sgtree.benchmarks.gunc_low_memory",
        "run",
        "--input_dir",
        str(unique_dir),
        "--file_suffix",
        ".fna",
        "--db_file",
        str(database),
        "--threads",
        str(threads),
        "--out_dir",
        str(output_dir),
        "--temp_dir",
        str(temp_dir),
        "--detailed_output",
        "--contig_taxonomy_output",
    ]
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    started = time.monotonic()
    with log_path.open("w") as log_handle:
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            env=environment,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    wall_seconds = time.monotonic() - started
    if result.returncode != 0:
        raise RuntimeError(f"GUNC failed; see {log_path}")
    summary_paths = _gunc_summary_paths(gunc_dir)
    if not summary_paths:
        raise RuntimeError("GUNC completed without a maxCSS summary")
    peak_rss_kb = _resource_value(
        resource_path,
        "Maximum resident set size (kbytes)",
    )
    if peak_rss_kb is None or peak_rss_kb > MAX_PEAK_RSS_KB:
        raise RuntimeError(
            f"GUNC peak RSS is {peak_rss_kb!r} kB; cap is {MAX_PEAK_RSS_KB} kB"
        )
    status = {
        "pass": True,
        "database_path": _relative(database),
        "database_sha256": database_sha256,
        "unique_assembly_count": len(unique_sources),
        "mapped_context_assembly_count": len(mapping_rows),
        "summary_paths": [_relative(path) for path in summary_paths],
        "wall_seconds": wall_seconds,
        "peak_rss_kb": int(peak_rss_kb),
        "threads": threads,
        "command": command,
    }
    _write_json(status_path, status)
    return status


def _boolean_series(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"true", "yes", "1"})


def _cluster_ratio(
    frame: pd.DataFrame,
    *,
    numerator_column: str,
    denominator_column: str,
) -> dict:
    grouped = (
        frame.groupby("panel_id")[[numerator_column, denominator_column]]
        .sum()
        .reindex(
            [f"{lineage}_p{index}_seed{seed}" for lineage in LINEAGES for index, seed in enumerate(PANEL_SEEDS, 1)],
            fill_value=0,
        )
    )
    numerator = int(grouped[numerator_column].sum())
    denominator = int(grouped[denominator_column].sum())
    estimate = numerator / denominator if denominator else None
    if not denominator:
        interval = (None, None)
    else:
        rng = random.Random(10_607)
        values = []
        rows = grouped.to_dict("records")
        for _ in range(BOOTSTRAP_REPLICATES):
            sample = [rows[rng.randrange(len(rows))] for _index in rows]
            sample_denominator = sum(
                int(row[denominator_column]) for row in sample
            )
            if sample_denominator:
                values.append(
                    sum(int(row[numerator_column]) for row in sample)
                    / sample_denominator
                )
        interval = (
            float(np.quantile(values, 0.025)),
            float(np.quantile(values, 0.975)),
        )
    return {
        "numerator": numerator,
        "denominator": denominator,
        "estimate": estimate,
        "cluster_bootstrap_ci_low": interval[0],
        "cluster_bootstrap_ci_high": interval[1],
        "independent_unit": "base_panel",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
    }


def _load_gunc_results(outdir: Path) -> pd.DataFrame:
    gunc_dir = outdir / "gunc"
    status_path = gunc_dir / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError("Run the frozen GUNC endpoint before analysis")
    status = json.loads(status_path.read_text())
    if status.get("pass") is not True:
        raise RuntimeError("Saved GUNC status is not passing")
    summaries = []
    for value in status["summary_paths"]:
        frame = pd.read_csv(_project_path(value), sep="\t")
        if "genome" not in frame.columns:
            raise ValueError(f"GUNC summary lacks genome column: {value}")
        frame = frame.copy()
        frame["gunc_input_id"] = (
            frame["genome"]
            .astype(str)
            .map(lambda value: Path(value).stem)
        )
        summaries.append(frame)
    scores = pd.concat(summaries, ignore_index=True)
    if scores["gunc_input_id"].duplicated().any():
        raise ValueError("GUNC returned duplicate rows for a unique assembly")
    mapping = pd.read_csv(gunc_dir / "assembly_map.tsv", sep="\t", keep_default_na=False)
    results = mapping.merge(
        scores.drop(columns=["genome"]),
        on="gunc_input_id",
        how="left",
        validate="many_to_one",
    )
    if "pass.GUNC" not in results.columns:
        raise ValueError("GUNC summary lacks pass.GUNC")
    if results["pass.GUNC"].isna().any():
        raise ValueError("Failed to map one or more GUNC result rows")
    results["gunc_contamination_call"] = ~_boolean_series(results["pass.GUNC"])
    _write_tsv(
        gunc_dir / "mapped_results.tsv",
        results.to_dict("records"),
        list(results.columns),
    )
    return results


def _metric_row(name: str, result: dict, unit: str = "fraction") -> dict:
    return {
        "metric": name,
        "numerator": result.get("numerator"),
        "denominator": result.get("denominator"),
        "estimate": result.get("estimate"),
        "cluster_bootstrap_ci_low": result.get("cluster_bootstrap_ci_low"),
        "cluster_bootstrap_ci_high": result.get("cluster_bootstrap_ci_high"),
        "unit": unit,
        "independent_unit": result.get("independent_unit", "base_panel"),
    }


def analyze_confirmation(outdir: Path) -> dict:
    manifest = _load_manifest(outdir)
    check_confirmation(outdir, require_full=True)
    frames = []
    context_frames = []
    panel_audits = []
    for panel in manifest["panels"]:
        panel_dir = outdir / "panels" / panel["panel_id"]
        frames.append(
            pd.read_csv(
                panel_dir / "per_cell_results.tsv",
                sep="\t",
                keep_default_na=False,
            )
        )
        context_frames.append(
            pd.read_csv(
                panel_dir / "context_summary.tsv",
                sep="\t",
                keep_default_na=False,
            )
        )
        panel_audits.append(json.loads((panel_dir / "panel_audit.json").read_text()))
    cells = pd.concat(frames, ignore_index=True)
    contexts = pd.concat(context_frames, ignore_index=True)
    for column in (
        "truth_positive",
        "raw_call",
        "contig_gate_pass",
        "budget_pass",
        "guarded_action",
        "false_action",
    ):
        cells[column] = _boolean_series(cells[column])

    event_rows = cells[cells["event_id"].astype(str).ne("")].copy()
    gunc = _load_gunc_results(outdir)
    gunc_columns = [
        column
        for column in (
            "panel_id",
            "context",
            "genome",
            "pass.GUNC",
            "gunc_contamination_call",
            "taxonomic_level",
            "clade_separation_score",
            "contamination_portion",
            "reference_representation_score",
            "n_genes_called",
            "n_genes_mapped",
            "n_contigs",
        )
        if column in gunc.columns
    ]
    event_rows = event_rows.merge(
        gunc[gunc_columns],
        left_on=["panel_id", "context", "recipient_genome"],
        right_on=["panel_id", "context", "genome"],
        how="left",
        validate="many_to_one",
    )
    if "genome" in event_rows.columns:
        event_rows = event_rows.drop(columns=["genome"])

    metrics = []
    truth = cells[cells["truth_positive"]].copy()
    truth["denominator"] = 1
    truth["raw_numerator"] = truth["raw_call"].astype(int)
    truth["guarded_numerator"] = truth["guarded_action"].astype(int)
    raw_recall = _cluster_ratio(
        truth,
        numerator_column="raw_numerator",
        denominator_column="denominator",
    )
    guarded_recall = _cluster_ratio(
        truth,
        numerator_column="guarded_numerator",
        denominator_column="denominator",
    )
    metrics.extend(
        [
            _metric_row("loo_raw_gene_rich_recall", raw_recall),
            _metric_row("guarded_gene_rich_action_recall", guarded_recall),
        ]
    )
    raw_calls = cells.copy()
    raw_calls["precision_denominator"] = raw_calls["raw_call"].astype(int)
    raw_calls["precision_numerator"] = (
        raw_calls["raw_call"] & raw_calls["truth_positive"]
    ).astype(int)
    metrics.append(
        _metric_row(
            "loo_raw_candidate_precision",
            _cluster_ratio(
                raw_calls,
                numerator_column="precision_numerator",
                denominator_column="precision_denominator",
            ),
        )
    )
    for stratum in ("near", "intermediate", "far"):
        subset = truth[truth["source_distance_stratum"] == stratum]
        metrics.append(
            _metric_row(
                f"loo_raw_{stratum}_recall",
                _cluster_ratio(
                    subset,
                    numerator_column="raw_numerator",
                    denominator_column="denominator",
                ),
            )
        )
        metrics.append(
            _metric_row(
                f"guarded_{stratum}_action_recall",
                _cluster_ratio(
                    subset,
                    numerator_column="guarded_numerator",
                    denominator_column="denominator",
                ),
            )
        )

    actions = cells.copy()
    actions["precision_denominator"] = actions["guarded_action"].astype(int)
    actions["precision_numerator"] = (
        actions["guarded_action"] & actions["truth_positive"]
    ).astype(int)
    metrics.append(
        _metric_row(
            "guarded_action_precision",
            _cluster_ratio(
                actions,
                numerator_column="precision_numerator",
                denominator_column="precision_denominator",
            ),
        )
    )
    gunc_clean = gunc[gunc["context"] == "clean"].copy()
    gunc_clean["denominator"] = 1
    gunc_clean["gunc_numerator"] = _boolean_series(
        gunc_clean["gunc_contamination_call"]
    ).astype(int)
    metrics.append(
        _metric_row(
            "gunc_clean_genome_call_rate",
            _cluster_ratio(
                gunc_clean,
                numerator_column="gunc_numerator",
                denominator_column="denominator",
            ),
        )
    )
    for metric_name, controls in (
        (
            "false_guarded_actions_per_1000_clean_cells",
            cells[cells["context"] == "clean"].copy(),
        ),
        (
            "false_guarded_actions_per_1000_unaffected_cells",
            cells[~cells["truth_positive"]].copy(),
        ),
    ):
        controls["control_denominator"] = 1
        controls["false_numerator"] = controls["guarded_action"].astype(int)
        false_rate = _cluster_ratio(
            controls,
            numerator_column="false_numerator",
            denominator_column="control_denominator",
        )
        false_rate["estimate"] = (
            false_rate["estimate"] * 1000
            if false_rate["estimate"] is not None
            else None
        )
        for key in ("cluster_bootstrap_ci_low", "cluster_bootstrap_ci_high"):
            if false_rate[key] is not None:
                false_rate[key] *= 1000
        metrics.append(
            _metric_row(
                metric_name,
                false_rate,
                unit="actions_per_1000_cells",
            )
        )
    truth["abstention_numerator"] = (
        truth["loo_abstention_reason"].astype(str).str.len() > 0
    ).astype(int)
    metrics.append(
        _metric_row(
            "gene_rich_abstention_coverage",
            _cluster_ratio(
                truth,
                numerator_column="abstention_numerator",
                denominator_column="denominator",
            ),
        )
    )

    event_retention = (
        event_rows.groupby(
            ["panel_id", "event_id", "event_class"],
            as_index=False,
        )["guarded_action"]
        .max()
    )
    for event_class, metric_name in (
        ("solo_marker_control", "solo_control_retention"),
        ("native_contig_sentinel", "native_contig_sentinel_retention"),
    ):
        subset = event_retention[event_retention["event_class"] == event_class].copy()
        subset["denominator"] = 1
        subset["retained"] = (~_boolean_series(subset["guarded_action"])).astype(int)
        metrics.append(
            _metric_row(
                metric_name,
                _cluster_ratio(
                    subset,
                    numerator_column="retained",
                    denominator_column="denominator",
                ),
            )
        )

    gunc_event_rows = event_rows.copy()
    gunc_event_rows["gunc_call"] = _boolean_series(
        gunc_event_rows["gunc_contamination_call"]
    )
    gunc_truth = gunc_event_rows[
        gunc_event_rows["event_class"] == "gene_rich_replacement"
    ].copy()
    gunc_truth["denominator"] = 1
    gunc_truth["gunc_numerator"] = gunc_truth["gunc_call"].astype(int)
    metrics.append(
        _metric_row(
            "gunc_gene_rich_genome_call_recall",
            _cluster_ratio(
                gunc_truth,
                numerator_column="gunc_numerator",
                denominator_column="denominator",
            ),
        )
    )
    for stratum in ("near", "intermediate", "far"):
        subset = gunc_truth[
            gunc_truth["source_distance_stratum"] == stratum
        ]
        metrics.append(
            _metric_row(
                f"gunc_{stratum}_genome_call_recall",
                _cluster_ratio(
                    subset,
                    numerator_column="gunc_numerator",
                    denominator_column="denominator",
                ),
            )
        )
    gunc_control_events = (
        gunc_event_rows.groupby(
            ["panel_id", "event_id", "event_class"],
            as_index=False,
        )["gunc_call"]
        .max()
    )
    for event_class, metric_name in (
        ("solo_marker_control", "gunc_solo_control_genome_call_rate"),
        (
            "native_contig_sentinel",
            "gunc_native_contig_sentinel_genome_call_rate",
        ),
    ):
        subset = gunc_control_events[
            gunc_control_events["event_class"] == event_class
        ].copy()
        subset["denominator"] = 1
        subset["gunc_numerator"] = _boolean_series(
            subset["gunc_call"]
        ).astype(int)
        metrics.append(
            _metric_row(
                metric_name,
                _cluster_ratio(
                    subset,
                    numerator_column="gunc_numerator",
                    denominator_column="denominator",
                ),
            )
        )

    panel_metrics = []
    for panel_id, group in cells.groupby("panel_id"):
        truth_group = group[group["truth_positive"]]
        panel_metrics.append(
            {
                "panel_id": panel_id,
                "lineage": str(group.iloc[0]["lineage"]),
                "seed": int(group.iloc[0]["seed"]),
                "gene_rich_truth_count": len(truth_group),
                "gene_rich_raw_call_count": int(truth_group["raw_call"].sum()),
                "gene_rich_guarded_action_count": int(
                    truth_group["guarded_action"].sum()
                ),
                "guarded_action_count": int(group["guarded_action"].sum()),
                "false_action_count": int(group["false_action"].sum()),
                "any_false_action": bool(group["false_action"].any()),
            }
        )

    results_dir = outdir / "results"
    results_dir.mkdir(exist_ok=True)
    _write_tsv(
        results_dir / "per_cell_results.tsv",
        cells.to_dict("records"),
        list(cells.columns),
    )
    _write_tsv(
        results_dir / "per_event_results.tsv",
        event_rows.to_dict("records"),
        list(event_rows.columns),
    )
    _write_tsv(
        results_dir / "context_summary.tsv",
        contexts.to_dict("records"),
        list(contexts.columns),
    )
    _write_tsv(
        results_dir / "panel_metrics.tsv",
        panel_metrics,
        list(panel_metrics[0]),
    )
    _write_tsv(
        results_dir / "summary.tsv",
        metrics,
        list(metrics[0]),
    )
    status = json.loads((outdir / "run_status.json").read_text())
    gunc_status = json.loads((outdir / "gunc/status.json").read_text())
    summary = {
        "story": "US-010",
        "manifest_content_sha256": manifest["manifest_content_sha256"],
        "independent_unit": "base_panel",
        "base_panel_count": 6,
        "context_count": 30,
        "gene_rich_event_count": 18,
        "solo_control_count": 6,
        "native_contig_sentinel_count": 6,
        "metrics": {row["metric"]: row for row in metrics},
        "panels_with_any_false_action": sum(
            row["any_false_action"] for row in panel_metrics
        ),
        "marker_tree_build_count": status["marker_tree_cache_entries"],
        "sgtree_peak_rss_kb": status["peak_rss_kb"],
        "sgtree_wall_seconds": status["wall_seconds"],
        "gunc_unique_assembly_count": gunc_status["unique_assembly_count"],
        "gunc_peak_rss_kb": gunc_status["peak_rss_kb"],
        "gunc_wall_seconds": gunc_status["wall_seconds"],
        "production_pruning": "disabled",
        "gunc_endpoint_note": (
            "GUNC genome/chimerism calls are reported separately from marker-level "
            "LOO and guarded actions; no cross-tool F1 is calculated."
        ),
        "interval_method": (
            "Deterministic percentile cluster bootstrap over six base panels "
            f"({BOOTSTRAP_REPLICATES} resamples); contexts and events remain nested."
        ),
        "tree_change_by_context": {
            context: {
                "panel_count": len(group),
                "normalized_rf_mean": float(
                    pd.to_numeric(group["normalized_rf_vs_clean"]).mean()
                ),
                "normalized_rf_median": float(
                    pd.to_numeric(group["normalized_rf_vs_clean"]).median()
                ),
                "normalized_patristic_error_mean": float(
                    pd.to_numeric(
                        group["normalized_patristic_error_vs_clean"]
                    ).mean()
                ),
                "normalized_patristic_error_median": float(
                    pd.to_numeric(
                        group["normalized_patristic_error_vs_clean"]
                    ).median()
                ),
            }
            for context, group in contexts.groupby("context")
        },
    }
    _write_json(results_dir / "summary.json", summary)
    return summary


def package_confirmation(outdir: Path) -> dict:
    manifest = _load_manifest(outdir)
    check_confirmation(outdir, require_full=True)
    summary_path = outdir / "results/summary.json"
    gunc_status_path = outdir / "gunc/status.json"
    if not summary_path.is_file() or not gunc_status_path.is_file():
        raise FileNotFoundError("Run GUNC and analysis before packaging")

    evidence_dir = outdir / "evidence"
    evidence_dir.mkdir(exist_ok=True)
    shutil.copy2(PROJECT_ROOT / "pixi.toml", evidence_dir / "pixi.toml")
    shutil.copy2(PROJECT_ROOT / "pixi.lock", evidence_dir / "pixi.lock")
    shutil.copy2(PROJECT_ROOT / "PRD.md", evidence_dir / "PRD.md")
    source_archive = evidence_dir / "source-and-tests.tar.gz"
    archive_members = [
        (path, _relative(path))
        for base in (PROJECT_ROOT / "src", PROJECT_ROOT / "tests")
        for path in sorted(base.rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    ]
    _write_deterministic_tar_gz(source_archive, archive_members)
    benchmark_command = (
        "env PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 "
        "pixi run python -m sgtree.benchmarks.fna_confirmation "
        f"--outdir {_relative(outdir)}"
    )
    commands = "\n".join(
        [
            "# Run from the repository root with the archived Pixi files restored.",
            "mkdir -p resources/gunc_db",
            "pixi run gunc download_db resources/gunc_db -db progenomes_2.1",
            f"{benchmark_command} freeze "
            "--models resources/models/UNI56.hmm --gunc-db resources/gunc_db",
            f"{benchmark_command} pilot --threads 4",
            f"{benchmark_command} run --threads 4",
            f"{benchmark_command} gunc --threads 4",
            f"{benchmark_command} analyze",
            f"{benchmark_command} check --full",
            f"{benchmark_command} package",
            "",
            "# To rebuild only the result tables from the archived run artifacts:",
            f"{benchmark_command} analyze",
            "",
        ]
    )
    (evidence_dir / "commands.txt").write_text(commands)

    checksum_rows = []
    checksum_path = evidence_dir / "checksums.tsv"
    for path in sorted(outdir.rglob("*")):
        if (
            not path.is_file()
            or path == checksum_path
            or path == evidence_dir / "package_status.json"
        ):
            continue
        checksum_rows.append(
            {
                "path": str(path.relative_to(outdir)),
                "sha256": _sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    _write_tsv(
        checksum_path,
        checksum_rows,
        ["path", "sha256", "bytes"],
    )
    status = {
        "pass": True,
        "manifest_content_sha256": manifest["manifest_content_sha256"],
        "source_code_sha256": manifest["code"]["source_code_sha256"],
        "pixi_lock_sha256": manifest["pixi_lock_sha256"],
        "file_count": len(checksum_rows),
        "checksums_sha256": _sha256_file(checksum_path),
        "per_event_results_sha256": _sha256_file(
            outdir / "results/per_event_results.tsv"
        ),
        "summary_json_sha256": _sha256_file(summary_path),
        "source_archive_sha256": _sha256_file(source_archive),
        "production_pruning": "disabled",
    }
    _write_json(evidence_dir / "package_status.json", status)
    return status


def check_confirmation(outdir: Path, *, require_full: bool) -> dict:
    manifest = _load_manifest(outdir)
    status_path = outdir / ("run_status.json" if require_full else "pilot_status.json")
    if not status_path.is_file():
        raise FileNotFoundError(status_path)
    status = json.loads(status_path.read_text())
    expected_panels = 6 if require_full else 1
    errors = []
    if status.get("pass") is not True:
        errors.append("saved status is not passing")
    if int(status.get("panels_completed", 0)) != expected_panels:
        errors.append(
            f"expected {expected_panels} completed panels; "
            f"found {status.get('panels_completed')}"
        )
    expected_contexts = expected_panels * len(CONTEXT_ORDER)
    if int(status.get("contexts_completed", 0)) != expected_contexts:
        errors.append(
            f"expected {expected_contexts} completed contexts; "
            f"found {status.get('contexts_completed')}"
        )
    cache_count = len(list((outdir / "marker_tree_cache").glob("*.nwk")))
    if cache_count > MAX_MARKER_TREE_BUILDS:
        errors.append(
            f"marker-tree build count {cache_count} exceeds {MAX_MARKER_TREE_BUILDS}"
        )
    target_cache_count = (
        EXPECTED_MARKER_TREE_BUILDS
        if require_full
        else EXPECTED_MARKER_TREE_BUILDS // len(manifest["panels"])
    )
    minimum_cache_count = expected_panels * MARKER_COUNT
    if not minimum_cache_count <= cache_count <= target_cache_count:
        errors.append(
            f"expected {minimum_cache_count} to {target_cache_count} distinct "
            f"marker-tree builds; found {cache_count}"
        )
    if int(status.get("peak_rss_kb", 0)) > MAX_PEAK_RSS_KB:
        errors.append("peak RSS exceeds the frozen 8 GB cap")
    checked_panels = manifest["panels"][:expected_panels]
    for panel in checked_panels:
        audit_path = outdir / "panels" / panel["panel_id"] / "panel_audit.json"
        if not audit_path.is_file() or json.loads(audit_path.read_text()).get("pass") is not True:
            errors.append(f"missing or failing panel audit: {panel['panel_id']}")
    result = {
        "pass": not errors,
        "require_full": require_full,
        "panels_checked": expected_panels,
        "marker_tree_cache_entries": cache_count,
        "errors": errors,
    }
    if errors:
        raise RuntimeError("; ".join(errors))
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Freeze and run the bounded US-010 empirical-FNA confirmation."
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    subparsers = parser.add_subparsers(dest="command", required=True)

    freeze = subparsers.add_parser("freeze")
    freeze.add_argument("--models", type=Path, default=DEFAULT_MODELS)
    freeze.add_argument("--gunc-db", type=Path, default=DEFAULT_GUNC_DB_DIR)

    for name in ("pilot", "run"):
        command = subparsers.add_parser(name)
        command.add_argument("--threads", type=int, default=MAX_THREADS)

    gunc = subparsers.add_parser("gunc")
    gunc.add_argument("--threads", type=int, default=MAX_THREADS)
    subparsers.add_parser("analyze")
    subparsers.add_parser("package")

    check = subparsers.add_parser("check")
    check.add_argument("--full", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    outdir = args.outdir.resolve()
    if args.command == "freeze":
        result = {
            "manifest": str(
                freeze_manifest(
                    outdir,
                    args.models.resolve(),
                    args.gunc_db.resolve(),
                )
            )
        }
    elif args.command == "pilot":
        result = run_confirmation(
            outdir,
            threads=args.threads,
            pilot_only=True,
        )
    elif args.command == "run":
        result = run_confirmation(
            outdir,
            threads=args.threads,
            pilot_only=False,
        )
    elif args.command == "gunc":
        result = run_gunc(outdir, threads=args.threads)
    elif args.command == "analyze":
        result = analyze_confirmation(outdir)
    elif args.command == "package":
        result = package_confirmation(outdir)
    else:
        result = check_confirmation(outdir, require_full=args.full)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
