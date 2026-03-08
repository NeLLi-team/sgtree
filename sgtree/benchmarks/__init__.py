"""Synthetic benchmark generation and evaluation for SGTree."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from itertools import combinations
from pathlib import Path
from random import Random

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from pyhmmer import plan7


DEFAULT_SCENARIOS = {
    "duplicate_only": {
        "pair_blocks": 2,
        "markers_per_block": 2,
        "replacement_events": 0,
        "native_degrade_fraction": 0.08,
    },
    "replacement_only": {
        "pair_blocks": 0,
        "markers_per_block": 0,
        "replacement_events": 4,
        "native_degrade_fraction": 0.12,
    },
    "combined": {
        "pair_blocks": 2,
        "markers_per_block": 2,
        "replacement_events": 4,
        "native_degrade_fraction": 0.12,
    },
}

DEFAULT_CLEANUP_PROFILES = {
    "duplicate_only": {
        "name": "duplicate_cleanup",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": False,
        "singles_mode": "neighbor",
    },
    "replacement_only": {
        "name": "singles_neighbor",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "neighbor",
    },
    "combined": {
        "name": "duplicate_plus_singles_neighbor",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "neighbor",
    },
}


def make_contaminant_record(
    recipient_genome: str,
    donor_record: SeqRecord,
    marker: str,
    donor_genome: str,
    event_index: int,
) -> SeqRecord:
    record = deepcopy(donor_record)
    record.id = (
        f"{recipient_genome}|contam__{marker}__{donor_genome}__e{event_index:03d}"
    )
    record.name = record.id
    record.description = record.id
    return record


def apply_replacement_event(
    recipient_records: dict[str, SeqRecord],
    native_record_id: str,
    contaminant_record: SeqRecord,
) -> dict[str, SeqRecord]:
    updated = dict(recipient_records)
    if native_record_id not in updated:
        raise KeyError(f"Native marker record not found: {native_record_id}")
    del updated[native_record_id]
    updated[contaminant_record.id] = contaminant_record
    return updated


def _degrade_record_in_place(
    recipient_records: dict[str, SeqRecord],
    native_record_id: str,
    fraction: float,
    rng: Random,
) -> None:
    if fraction <= 0:
        return
    record = deepcopy(recipient_records[native_record_id])
    sequence = list(str(record.seq))
    n_mutations = max(1, round(len(sequence) * fraction))
    n_mutations = min(n_mutations, len(sequence))
    for idx in rng.sample(range(len(sequence)), n_mutations):
        sequence[idx] = "X"
    record.seq = Seq("".join(sequence))
    recipient_records[native_record_id] = record


def _run_cmd(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def _run_sgtree_python(
    genomedir: Path,
    modeldir: Path,
    outdir: Path,
    *,
    num_cpus: int,
    percent_models: int,
    marker_selection: bool,
    singles: bool,
    singles_mode: str,
    selection_mode: str,
    selection_max_rounds: int,
    selection_global_rounds: int,
    keep_intermediates: bool,
) -> None:
    cmd = [
        sys.executable,
        "-m",
        "sgtree",
        str(genomedir),
        str(modeldir),
        "--num_cpus",
        str(num_cpus),
        "--percent_models",
        str(percent_models),
        "--save_dir",
        str(outdir),
        "--selection_mode",
        selection_mode,
        "--selection_max_rounds",
        str(selection_max_rounds),
        "--selection_global_rounds",
        str(selection_global_rounds),
        "--keep_intermediates",
        "yes" if keep_intermediates else "no",
    ]
    if marker_selection:
        cmd.extend(["--marker_selection", "yes"])
    if singles:
        cmd.extend(["--singles", "yes"])
    cmd.extend(["--singles_mode", singles_mode])
    _run_cmd(cmd)


def _read_normalized_proteomes(path: Path) -> dict[str, dict[str, SeqRecord]]:
    by_genome: dict[str, dict[str, SeqRecord]] = {}
    inputs = sorted(entry for entry in path.iterdir() if entry.is_file()) if path.is_dir() else [path]
    for input_path in inputs:
        with open(input_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                genome = record.id.split("|")[0]
                by_genome.setdefault(genome, {})[record.id] = deepcopy(record)
    return by_genome


def _write_proteome_dir(records_by_genome: dict[str, dict[str, SeqRecord]], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for genome, records in sorted(records_by_genome.items()):
        out_path = outdir / f"{genome}.faa"
        with open(out_path, "w") as handle:
            for record_id in sorted(records):
                SeqIO.write(records[record_id], handle, "fasta")


def _load_table(table_path: Path) -> pd.DataFrame:
    df = pd.read_csv(table_path)
    score_col = "score_bits" if "score_bits" in df.columns else "7"
    if score_col not in df.columns:
        raise ValueError(f"Missing score column in {table_path}")
    marker_col = "namemodel"
    if marker_col not in df.columns:
        raise ValueError(f"Missing marker column '{marker_col}' in {table_path}")
    df = df.copy()
    df["genome"] = df["savedname"].str.split("/").str[0]
    df["marker"] = df[marker_col].str.split("/").str[-1]
    df["bitscore"] = pd.to_numeric(df[score_col], errors="coerce")
    if df["bitscore"].isna().any():
        raise ValueError(f"Failed to parse one or more bitscores from {table_path}")
    return df


def _count_single_copy_markers(df: pd.DataFrame, genomes: list[str]) -> pd.DataFrame:
    subset = df[df["genome"].isin(genomes)]
    counts = subset.groupby(["marker", "genome"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=genomes, fill_value=0)
    return counts


def _choose_genome_subset(
    df: pd.DataFrame,
    genome_sizes: dict[str, int],
    n_genomes: int,
) -> tuple[list[str], pd.DataFrame]:
    all_genomes = sorted(df["genome"].unique())
    if len(all_genomes) < n_genomes:
        counts = _count_single_copy_markers(df, all_genomes)
        return all_genomes, counts

    best_genomes: list[str] | None = None
    best_counts: pd.DataFrame | None = None
    best_score: tuple[int, int] | None = None

    n_to_exclude = len(all_genomes) - n_genomes
    if len(all_genomes) <= 25 and n_to_exclude <= 5:
        candidates = (
            [g for g in all_genomes if g not in excluded]
            for excluded in combinations(all_genomes, n_to_exclude)
        )
    else:
        ranked = sorted(all_genomes, key=lambda genome: (genome_sizes.get(genome, 0), genome))
        candidates = [ranked[:n_genomes]]

    for genomes in candidates:
        counts = _count_single_copy_markers(df, genomes)
        universal = counts[(counts == 1).all(axis=1)]
        score = (len(universal), -sum(genome_sizes.get(genome, 0) for genome in genomes))
        if best_score is None or score > best_score:
            best_score = score
            best_genomes = list(genomes)
            best_counts = counts

    assert best_genomes is not None and best_counts is not None
    return sorted(best_genomes), best_counts


def _rank_markers(
    df: pd.DataFrame,
    normalized_records: dict[str, dict[str, SeqRecord]],
    genomes: list[str],
    counts: pd.DataFrame,
    n_markers: int,
    min_presence_fraction: float,
) -> tuple[list[str], pd.DataFrame]:
    min_presence = max(1, int(len(genomes) * min_presence_fraction + 0.999999))
    marker_stats = pd.DataFrame(
        {
            "present_genomes": (counts > 0).sum(axis=1),
            "duplicated_genomes": (counts > 1).sum(axis=1),
        }
    )
    eligible = marker_stats[
        (marker_stats["duplicated_genomes"] == 0)
        & (marker_stats["present_genomes"] >= min_presence)
    ]
    if len(eligible) < n_markers:
        raise ValueError(
            f"Only found {len(eligible)} single-copy markers with presence >= {min_presence}/{len(genomes)}"
        )

    rows = []
    subset = df[df["genome"].isin(genomes)]
    for marker in eligible.index:
        marker_rows = subset[subset["marker"] == marker]
        lengths = []
        for row in marker_rows.itertuples(index=False):
            record_id = row.savedname.replace("/", "|")
            lengths.append(len(normalized_records[row.genome][record_id].seq))
        rows.append(
            {
                "marker": marker,
                "present_genomes": int(eligible.loc[marker, "present_genomes"]),
                "median_bitscore": float(marker_rows["bitscore"].median()),
                "min_bitscore": float(marker_rows["bitscore"].min()),
                "length_cv": float(pd.Series(lengths).std(ddof=0) / max(1.0, pd.Series(lengths).mean())),
            }
        )

    ranking = pd.DataFrame(rows).sort_values(
        ["present_genomes", "median_bitscore", "min_bitscore", "length_cv", "marker"],
        ascending=[False, False, False, True, True],
    )
    return ranking.head(n_markers)["marker"].tolist(), ranking


def _extract_hmm_subset(models_path: Path, markers: list[str], out_path: Path) -> None:
    requested = set(markers)
    found = set()
    with plan7.HMMFile(models_path) as hmm_file, open(out_path, "wb") as handle:
        for hmm_profile in hmm_file:
            name = hmm_profile.name
            if isinstance(name, bytes):
                name = name.decode("utf-8", errors="replace")
            marker = str(name)
            if marker in requested:
                hmm_profile.write(handle)
                found.add(marker)
    missing = sorted(requested - found)
    if missing:
        raise ValueError(f"Missing markers in HMM set: {', '.join(missing)}")


def _native_marker_map(df: pd.DataFrame, genomes: list[str], markers: list[str]) -> dict[str, dict[str, str]]:
    subset = df[df["genome"].isin(genomes) & df["marker"].isin(markers)]
    mapping: dict[str, dict[str, str]] = {}
    for row in subset.itertuples(index=False):
        mapping.setdefault(row.genome, {})[row.marker] = row.savedname.replace("/", "|")
    return mapping


def _load_truth_tree(path: Path) -> Tree:
    return Tree(str(path))


def _distance_ranked_donors(truth_tree: Tree, recipient: str, candidates: list[str]) -> list[str]:
    ranked = sorted(
        candidates,
        key=lambda genome: (-truth_tree.get_distance(recipient, genome), genome),
    )
    return ranked


def _choose_markers_for_pair(
    used_pairs: set[tuple[str, str]],
    genomes: tuple[str, str],
    markers: list[str],
    n_needed: int,
    rng: Random,
) -> list[str]:
    available = [
        marker
        for marker in markers
        if all((genome, marker) not in used_pairs for genome in genomes)
    ]
    if len(available) <= n_needed:
        return available
    return sorted(rng.sample(available, n_needed))


def _pair_distant_genomes(
    truth_tree: Tree,
    genomes: list[str],
    n_pairs: int,
    rng: Random,
) -> list[tuple[str, str]]:
    remaining = set(genomes)
    pairs: list[tuple[str, str]] = []
    while len(pairs) < n_pairs and len(remaining) >= 2:
        seed = rng.choice(sorted(remaining))
        partner = max(
            (genome for genome in remaining if genome != seed),
            key=lambda genome: (truth_tree.get_distance(seed, genome), genome),
        )
        pairs.append((seed, partner))
        remaining.remove(seed)
        remaining.remove(partner)
    return pairs


def _write_manifest_json(path: Path, payload: dict) -> None:
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_events_tsv(path: Path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)


def _source_relation(recipient_group: str, donor_group: str) -> str:
    return "within_group" if recipient_group == donor_group else "cross_group"


def _write_genome_summary_tsv(path: Path, rows: list[dict]) -> pd.DataFrame:
    columns = [
        "recipient_genome",
        "recipient_group",
        "contaminant_markers_added",
        "duplicate_events",
        "replacement_events",
        "within_group_events",
        "cross_group_events",
    ]
    if not rows:
        summary = pd.DataFrame(columns=columns)
    else:
        df = pd.DataFrame(rows).copy()
        df["duplicate_events"] = (df["event_type"] == "duplicate").astype(int)
        df["replacement_events"] = (df["event_type"] == "replacement").astype(int)
        df["within_group_events"] = (df["source_relation"] == "within_group").astype(int)
        df["cross_group_events"] = (df["source_relation"] == "cross_group").astype(int)
        summary = (
            df.groupby(["recipient_genome", "recipient_group"], as_index=False)
            .agg(
                contaminant_markers_added=("event_index", "count"),
                duplicate_events=("duplicate_events", "sum"),
                replacement_events=("replacement_events", "sum"),
                within_group_events=("within_group_events", "sum"),
                cross_group_events=("cross_group_events", "sum"),
            )
            .sort_values(["recipient_group", "recipient_genome"])
        )
    summary.to_csv(path, sep="\t", index=False)
    return summary


def _scenario_summary(events: list[dict], genome_summary: pd.DataFrame) -> dict:
    event_df = pd.DataFrame(events)
    if event_df.empty:
        duplicate_events = 0
        replacement_events = 0
        within_group_events = 0
        cross_group_events = 0
    else:
        duplicate_events = int((event_df["event_type"] == "duplicate").sum())
        replacement_events = int((event_df["event_type"] == "replacement").sum())
        within_group_events = int((event_df["source_relation"] == "within_group").sum())
        cross_group_events = int((event_df["source_relation"] == "cross_group").sum())
    return {
        "event_count": int(len(events)),
        "contaminant_markers_added": int(len(events)),
        "duplicate_events": duplicate_events,
        "replacement_events": replacement_events,
        "affected_genomes": int(len(genome_summary)),
        "mean_contaminants_per_affected_genome": float(genome_summary["contaminant_markers_added"].mean()) if not genome_summary.empty else 0.0,
        "max_contaminants_per_genome": int(genome_summary["contaminant_markers_added"].max()) if not genome_summary.empty else 0,
        "within_group_events": within_group_events,
        "cross_group_events": cross_group_events,
    }


def _infer_group_label(source_dir: Path) -> str:
    label = source_dir.name.lower()
    joined = str(source_dir).lower()
    if "flavo" in joined:
        return "flavo"
    if "gamma" in joined:
        return "gamma"
    return label or "all"


def _stage_dir_complete(stage_dir: Path) -> bool:
    required = [
        stage_dir / "proteomes",
        stage_dir / "table_elim_dups",
        stage_dir / "tree.nwk",
    ]
    return all(path.exists() for path in required)


def _materialize_benchmark_from_truth(
    *,
    truth_records: dict[str, dict[str, SeqRecord]],
    truth_markers: list[str],
    models_path: Path,
    outdir: Path,
    seed: int,
    num_cpus: int,
    group_labels: dict[str, str] | None = None,
    cross_group_only: bool = False,
    manifest_extra: dict | None = None,
) -> None:
    rng = Random(seed)
    outdir.mkdir(parents=True, exist_ok=True)

    truth_inputs = outdir / "truth_inputs"
    truth_run = outdir / "truth_run"
    scenarios_dir = outdir / "scenarios"
    for path in [truth_inputs, truth_run, scenarios_dir]:
        if path.exists():
            shutil.rmtree(path)
    truth_models = outdir / "truth_markers.hmm"
    if truth_models.exists():
        truth_models.unlink()
    _write_proteome_dir(truth_records, truth_inputs)

    _extract_hmm_subset(models_path, truth_markers, truth_models)

    _run_sgtree_python(
        genomedir=truth_inputs,
        modeldir=truth_models,
        outdir=truth_run,
        num_cpus=num_cpus,
        percent_models=100,
        marker_selection=False,
        singles=False,
        singles_mode="neighbor",
        selection_mode="coordinate",
        selection_max_rounds=5,
        selection_global_rounds=1,
        keep_intermediates=True,
    )

    table = _load_table(truth_run / "table_elim_dups")
    native_map = _native_marker_map(table, sorted(truth_records), truth_markers)
    truth_tree = _load_truth_tree(truth_run / "tree.nwk")
    selected_genomes = sorted(leaf.name for leaf in truth_tree.iter_leaves())
    truth_records = {genome: truth_records[genome] for genome in selected_genomes}
    native_map = {genome: native_map[genome] for genome in selected_genomes}

    scenarios_dir.mkdir(exist_ok=True)
    if group_labels is None:
        selected_group_labels = {genome: "all" for genome in selected_genomes}
    else:
        selected_group_labels = {genome: group_labels.get(genome, "all") for genome in selected_genomes}

    manifest = {
        "seed": seed,
        "models_path": str(models_path),
        "selected_genomes": selected_genomes,
        "selected_markers": truth_markers,
        "groups": selected_group_labels,
        "cross_group_only": cross_group_only,
        "scenarios": [],
    }
    if manifest_extra:
        manifest.update(manifest_extra)

    def _cross_pairs(n_pairs: int) -> list[tuple[str, str]]:
        by_group: dict[str, list[str]] = {}
        for genome, label in selected_group_labels.items():
            if genome in selected_genomes:
                by_group.setdefault(label, []).append(genome)
        labels = sorted(by_group)
        if len(labels) < 2:
            raise ValueError("Cross-group benchmark requires at least two groups")
        left = sorted(by_group[labels[0]])
        right = sorted(by_group[labels[1]])
        pairs = []
        used_left = set()
        used_right = set()
        while len(pairs) < n_pairs and len(used_left) < len(left) and len(used_right) < len(right):
            candidates = [
                (truth_tree.get_distance(a, b), a, b)
                for a in left if a not in used_left
                for b in right if b not in used_right
            ]
            if not candidates:
                break
            _dist, a, b = max(candidates, key=lambda item: (item[0], item[1], item[2]))
            pairs.append((a, b))
            used_left.add(a)
            used_right.add(b)
        return pairs

    for scenario_name, spec in DEFAULT_SCENARIOS.items():
        scenario_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
        used_pairs: set[tuple[str, str]] = set()
        events: list[dict] = []
        event_index = 1

        if spec["pair_blocks"] > 0:
            reciprocal_pairs = _cross_pairs(spec["pair_blocks"]) if cross_group_only else _pair_distant_genomes(truth_tree, selected_genomes, spec["pair_blocks"], rng)
            for recipient_a, recipient_b in reciprocal_pairs:
                available_markers = _choose_markers_for_pair(
                    used_pairs,
                    (recipient_a, recipient_b),
                    truth_markers,
                    spec["markers_per_block"],
                    rng,
                )
                if not available_markers:
                    continue
                for marker in available_markers:
                    for recipient, donor in ((recipient_a, recipient_b), (recipient_b, recipient_a)):
                        donor_record = truth_records[donor][native_map[donor][marker]]
                        _degrade_record_in_place(
                            scenario_records[recipient],
                            native_map[recipient][marker],
                            spec["native_degrade_fraction"],
                            rng,
                        )
                        contaminant = make_contaminant_record(
                            recipient_genome=recipient,
                            donor_record=donor_record,
                            marker=marker,
                            donor_genome=donor,
                            event_index=event_index,
                        )
                        scenario_records[recipient][contaminant.id] = contaminant
                        used_pairs.add((recipient, marker))
                        events.append(
                            {
                                "event_index": event_index,
                                "scenario": scenario_name,
                                "event_type": "duplicate",
                                "recipient_genome": recipient,
                                "recipient_group": selected_group_labels[recipient],
                                "marker": marker,
                                "native_record_id": native_map[recipient][marker],
                                "donor_genome": donor,
                                "donor_group": selected_group_labels[donor],
                                "source_relation": _source_relation(selected_group_labels[recipient], selected_group_labels[donor]),
                                "donor_record_id": native_map[donor][marker],
                                "contaminant_record_id": contaminant.id,
                                "expected_duplicate_status": "Removed",
                                "expected_native_status": "Kept",
                                "native_degrade_fraction": spec["native_degrade_fraction"],
                            }
                        )
                        event_index += 1

        replacement_pairs = [
            (recipient, marker)
            for recipient in selected_genomes
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
        ]
        rng.shuffle(replacement_pairs)
        for recipient, marker in replacement_pairs[: spec["replacement_events"]]:
            if cross_group_only:
                donor_group = "gamma" if selected_group_labels[recipient] == "flavo" else "flavo"
                donor_candidates = [
                    genome
                    for genome in selected_genomes
                    if selected_group_labels[genome] == donor_group
                ]
            else:
                donor_candidates = _distance_ranked_donors(
                    truth_tree,
                    recipient,
                    [g for g in selected_genomes if g != recipient],
                )
            if not donor_candidates:
                continue
            donor = donor_candidates[0]
            donor_record = truth_records[donor][native_map[donor][marker]]
            contaminant = make_contaminant_record(
                recipient_genome=recipient,
                donor_record=donor_record,
                marker=marker,
                donor_genome=donor,
                event_index=event_index,
            )
            _degrade_record_in_place(
                scenario_records[recipient],
                native_map[recipient][marker],
                spec["native_degrade_fraction"],
                rng,
            )
            scenario_records[recipient] = apply_replacement_event(
                scenario_records[recipient],
                native_record_id=native_map[recipient][marker],
                contaminant_record=contaminant,
            )
            events.append(
                {
                    "event_index": event_index,
                    "scenario": scenario_name,
                    "event_type": "replacement",
                    "recipient_genome": recipient,
                    "recipient_group": selected_group_labels[recipient],
                    "marker": marker,
                    "native_record_id": native_map[recipient][marker],
                    "donor_genome": donor,
                    "donor_group": selected_group_labels[donor],
                    "source_relation": _source_relation(selected_group_labels[recipient], selected_group_labels[donor]),
                    "donor_record_id": native_map[donor][marker],
                    "contaminant_record_id": contaminant.id,
                    "expected_replacement_outcome": "DropMarkerOrRemoveContaminant",
                    "native_degrade_fraction": spec["native_degrade_fraction"],
                }
            )
            event_index += 1

        scenario_dir = scenarios_dir / scenario_name
        proteomes_dir = scenario_dir / "proteomes"
        _write_proteome_dir(scenario_records, proteomes_dir)
        _write_events_tsv(scenario_dir / "events.tsv", events)
        genome_summary = _write_genome_summary_tsv(scenario_dir / "genome_summary.tsv", events)
        summary = _scenario_summary(events, genome_summary)
        manifest["scenarios"].append(
            {
                "name": scenario_name,
                "proteomes_dir": str(proteomes_dir),
                "events_path": str(scenario_dir / "events.tsv"),
                "genome_summary_path": str(scenario_dir / "genome_summary.tsv"),
                **summary,
            }
        )

    _write_manifest_json(outdir / "benchmark_manifest.json", manifest)


def generate_cross_benchmark_dataset(
    benchmark_a_dir: Path,
    benchmark_b_dir: Path,
    outdir: Path,
    *,
    n_markers: int,
    seed: int,
    num_cpus: int,
) -> None:
    manifest_a = json.loads((benchmark_a_dir / "benchmark_manifest.json").read_text())
    manifest_b = json.loads((benchmark_b_dir / "benchmark_manifest.json").read_text())

    truth_records_a = _read_normalized_proteomes(benchmark_a_dir / "truth_inputs")
    truth_records_b = _read_normalized_proteomes(benchmark_b_dir / "truth_inputs")
    truth_records = {**truth_records_a, **truth_records_b}
    if set(truth_records_a) & set(truth_records_b):
        raise ValueError("Overlapping genome ids between benchmark families")

    ranks_a = {row["marker"]: idx for idx, row in enumerate(manifest_a["marker_ranking"])}
    ranks_b = {row["marker"]: idx for idx, row in enumerate(manifest_b["marker_ranking"])}
    shared = sorted(
        set(ranks_a) & set(ranks_b),
        key=lambda marker: (ranks_a[marker] + ranks_b[marker], marker),
    )
    truth_markers = shared[:n_markers]
    if not truth_markers:
        raise ValueError("No shared markers between benchmark families")

    group_labels = {
        **{genome: "flavo" for genome in truth_records_a},
        **{genome: "gamma" for genome in truth_records_b},
    }
    _materialize_benchmark_from_truth(
        truth_records=truth_records,
        truth_markers=truth_markers,
        models_path=Path(manifest_a["models_path"]),
        outdir=outdir,
        seed=seed,
        num_cpus=num_cpus,
        group_labels=group_labels,
        cross_group_only=True,
        manifest_extra={
            "benchmark_a_dir": str(benchmark_a_dir),
            "benchmark_b_dir": str(benchmark_b_dir),
            "source_relationship": "cross_group_only",
        },
    )


def generate_benchmark_dataset(
    source_dir: Path,
    models_path: Path,
    outdir: Path,
    *,
    n_genomes: int,
    n_markers: int,
    min_marker_presence_fraction: float,
    seed: int,
    num_cpus: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    stage_dir = outdir / "stage_full_model_clean"
    if not _stage_dir_complete(stage_dir):
        _run_sgtree_python(
            genomedir=source_dir,
            modeldir=models_path,
            outdir=stage_dir,
            num_cpus=num_cpus,
            percent_models=0,
            marker_selection=False,
            singles=False,
            singles_mode="neighbor",
            selection_mode="coordinate",
            selection_max_rounds=5,
            selection_global_rounds=1,
            keep_intermediates=True,
        )

    normalized_records = _read_normalized_proteomes(stage_dir / "proteomes")
    genome_sizes = {genome: len(records) for genome, records in normalized_records.items()}
    table = _load_table(stage_dir / "table_elim_dups")
    selected_genomes, counts = _choose_genome_subset(table, genome_sizes, n_genomes)
    truth_markers, ranking = _rank_markers(
        table,
        normalized_records,
        selected_genomes,
        counts,
        n_markers,
        min_marker_presence_fraction,
    )
    truth_records = {genome: deepcopy(normalized_records[genome]) for genome in selected_genomes}
    group_label = _infer_group_label(source_dir)
    group_labels = {genome: group_label for genome in selected_genomes}
    _materialize_benchmark_from_truth(
        truth_records=truth_records,
        truth_markers=truth_markers,
        models_path=models_path,
        outdir=outdir,
        seed=seed,
        num_cpus=num_cpus,
        group_labels=group_labels,
        cross_group_only=False,
        manifest_extra={
            "source_dir": str(source_dir),
            "min_marker_presence_fraction": min_marker_presence_fraction,
            "marker_ranking": ranking.to_dict("records"),
            "dataset_group_label": group_label,
        },
    )


def _load_events(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def _load_rf_status(path: Path) -> dict[tuple[str, str], str]:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    return {(row.ProteinID, row.MarkerGene): row.Status for row in df.itertuples(index=False)}


def _normalized_status_id(record_id: str) -> str:
    return record_id.replace("|", "/")


def _rf_norm(truth_tree_path: Path, observed_tree_path: Path) -> float:
    truth = Tree(str(truth_tree_path))
    observed = Tree(str(observed_tree_path))
    shared = sorted(set(leaf.name for leaf in truth) & set(leaf.name for leaf in observed))
    truth.prune(shared)
    observed.prune(shared)
    rf, maxrf, *_ = truth.robinson_foulds(observed, unrooted_trees=True)
    return rf / maxrf if maxrf else 0.0


def _replacement_outcome(run_dir: Path, recipient: str, marker: str) -> str:
    aligned_path = run_dir / "aligned_final" / f"{marker}.faa"
    if not aligned_path.exists():
        return "unknown"
    with open(aligned_path) as handle:
        ids = [record.id for record in SeqIO.parse(handle, "fasta") if record.id.startswith(f"{recipient}|")]
    if not ids:
        return "marker_dropped"
    if any("|contam__" in record_id for record_id in ids):
        return "contaminant_retained"
    return "native_retained"


def evaluate_benchmark_run(
    benchmark_dir: Path,
    scenario_name: str,
    run_dir: Path,
    runtime_seconds: float,
) -> dict:
    truth_tree = benchmark_dir / "truth_run" / "tree.nwk"
    events = _load_events(benchmark_dir / "scenarios" / scenario_name / "events.tsv")
    result_tree = run_dir / "tree_final.nwk"
    if not result_tree.exists():
        result_tree = run_dir / "tree.nwk"
    initial_tree = run_dir / "tree.nwk"
    rf_status = _load_rf_status(run_dir / "marker_selection_rf_values.txt")

    duplicate_events = events[events["event_type"].isin(["duplicate", "triplicate"])]
    replacement_events = events[events["event_type"] == "replacement"]

    contaminant_correct = 0
    native_correct = 0
    for row in duplicate_events.itertuples(index=False):
        contaminant_correct += int(
            rf_status.get((_normalized_status_id(row.contaminant_record_id), row.marker)) == "Removed"
        )
        native_correct += int(
            rf_status.get((_normalized_status_id(row.native_record_id), row.marker)) == "Kept"
        )

    replacement_outcomes = [
        _replacement_outcome(run_dir, row.recipient_genome, row.marker)
        for row in replacement_events.itertuples(index=False)
    ]
    replacement_contaminant_retained = int(sum(outcome == "contaminant_retained" for outcome in replacement_outcomes))
    replacement_contaminant_removed = int(len(replacement_events) - replacement_contaminant_retained)
    total_contaminants = int(len(duplicate_events) + len(replacement_events))
    total_removed = int(contaminant_correct + replacement_contaminant_removed)

    return {
        "scenario": scenario_name,
        "run_dir": str(run_dir),
        "truth_tree_path": str(truth_tree),
        "initial_tree_path": str(initial_tree),
        "final_tree_path": str(result_tree),
        "initial_tree_rf_norm": _rf_norm(truth_tree, initial_tree),
        "tree_rf_norm": _rf_norm(truth_tree, result_tree),
        "tree_rf_delta": _rf_norm(truth_tree, initial_tree) - _rf_norm(truth_tree, result_tree),
        "duplicate_events": int(len(duplicate_events)),
        "duplicate_contaminant_removed": int(contaminant_correct),
        "duplicate_native_retained": int(native_correct),
        "replacement_events": int(len(replacement_events)),
        "replacement_marker_dropped": int(sum(outcome == "marker_dropped" for outcome in replacement_outcomes)),
        "replacement_contaminant_retained": replacement_contaminant_retained,
        "replacement_contaminant_removed": replacement_contaminant_removed,
        "replacement_native_retained": int(sum(outcome == "native_retained" for outcome in replacement_outcomes)),
        "contaminant_markers_added": total_contaminants,
        "contaminant_markers_removed": total_removed,
        "contaminant_markers_removed_fraction": total_removed / total_contaminants if total_contaminants else 0.0,
        "runtime_seconds": runtime_seconds,
        "status": "ok",
    }


def _write_report(results: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# SGTree Systematic Benchmark Report",
        "",
        "| Scenario | Cleanup Profile | Initial RF | Final RF | RF Delta | Contaminants Removed | Duplicate Removed | Replacement Removed | Runtime (s) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in results.sort_values(["scenario"]).itertuples(index=False):
        lines.append(
            f"| {row.scenario} | {row.cleanup_profile} | {row.initial_tree_rf_norm:.3f} | {row.tree_rf_norm:.3f} | {row.tree_rf_delta:.3f} | "
            f"{row.contaminant_markers_removed}/{row.contaminant_markers_added} | "
            f"{row.duplicate_contaminant_removed}/{row.duplicate_events} | "
            f"{row.replacement_contaminant_removed}/{row.replacement_events} | "
            f"{row.runtime_seconds:.1f} |"
        )
    out_path.write_text("\n".join(lines) + "\n")


def run_benchmark_suite(
    benchmark_dir: Path,
    *,
    num_cpus: int,
    selection_max_rounds: int,
) -> None:
    manifest = json.loads((benchmark_dir / "benchmark_manifest.json").read_text())
    benchmark_models = benchmark_dir / "truth_markers.hmm"
    results_dir = benchmark_dir / "results"
    if results_dir.exists():
        shutil.rmtree(results_dir)
    results_dir.mkdir(exist_ok=True)

    summary_rows = []
    for scenario in manifest["scenarios"]:
        scenario_name = scenario["name"]
        if scenario_name not in DEFAULT_CLEANUP_PROFILES:
            raise ValueError(f"No cleanup profile defined for scenario: {scenario_name}")
        profile = DEFAULT_CLEANUP_PROFILES[scenario_name]
        proteomes_dir = Path(scenario["proteomes_dir"])
        run_dir = results_dir / f"{scenario_name}__{profile['name']}"
        start = time.time()
        try:
            _run_sgtree_python(
                genomedir=proteomes_dir,
                modeldir=benchmark_models,
                outdir=run_dir,
                num_cpus=num_cpus,
                percent_models=70,
                marker_selection=profile["marker_selection"],
                singles=profile["singles"],
                singles_mode=profile["singles_mode"],
                selection_mode=profile["selection_mode"],
                selection_max_rounds=selection_max_rounds,
                selection_global_rounds=profile["selection_global_rounds"],
                keep_intermediates=True,
            )
            runtime_seconds = time.time() - start
            row = evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name=scenario_name,
                run_dir=run_dir,
                runtime_seconds=runtime_seconds,
            )
            row["cleanup_profile"] = profile["name"]
            row["strategy"] = profile["name"]
        except subprocess.CalledProcessError as exc:
            row = {
                "cleanup_profile": profile["name"],
                "strategy": profile["name"],
                "scenario": scenario_name,
                "run_dir": str(run_dir),
                "status": "failed",
                "error": str(exc),
                "runtime_seconds": time.time() - start,
            }
        summary_rows.append(row)

    results = pd.DataFrame(summary_rows)
    results.to_csv(results_dir / "summary.tsv", sep="\t", index=False)
    ok_results = results[results["status"] == "ok"].copy()
    if not ok_results.empty:
        _write_report(ok_results, results_dir / "report.md")


def aggregate_benchmark_runs(benchmark_dirs: list[Path], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    frames = []
    for benchmark_dir in benchmark_dirs:
        summary_path = benchmark_dir / "results" / "summary.tsv"
        if not summary_path.exists():
            continue
        df = pd.read_csv(summary_path, sep="\t")
        df["benchmark_dir"] = str(benchmark_dir)
        frames.append(df)
    if not frames:
        raise FileNotFoundError("No benchmark summary.tsv files found for aggregation")
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(outdir / "combined_summary.tsv", sep="\t", index=False)
    ok = combined[combined["status"] == "ok"].copy()
    agg = (
        ok.groupby(["strategy", "scenario"], as_index=False)
        .agg(
            n_runs=("strategy", "size"),
            mean_initial_rf=("initial_tree_rf_norm", "mean"),
            mean_final_rf=("tree_rf_norm", "mean"),
            mean_rf_delta=("tree_rf_delta", "mean"),
            mean_dup_removed=("duplicate_contaminant_removed", "mean"),
            mean_dup_events=("duplicate_events", "mean"),
            mean_native_retained=("duplicate_native_retained", "mean"),
            mean_repl_dropped=("replacement_marker_dropped", "mean"),
            mean_repl_retained=("replacement_contaminant_retained", "mean"),
            mean_runtime=("runtime_seconds", "mean"),
        )
    )
    agg.to_csv(outdir / "aggregate_summary.tsv", sep="\t", index=False)

    lines = [
        "# Aggregate Benchmark Report",
        "",
        "| Strategy | Scenario | Runs | Initial RF | Final RF | RF Delta | Dup Removed | Native Retained | Replacement Dropped | Replacement Retained | Runtime (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in agg.sort_values(["strategy", "scenario"]).itertuples(index=False):
        lines.append(
            f"| {row.strategy} | {row.scenario} | {row.n_runs} | "
            f"{row.mean_initial_rf:.3f} | {row.mean_final_rf:.3f} | {row.mean_rf_delta:.3f} | "
            f"{row.mean_dup_removed:.1f}/{row.mean_dup_events:.1f} | "
            f"{row.mean_native_retained:.1f}/{row.mean_dup_events:.1f} | "
            f"{row.mean_repl_dropped:.1f} | {row.mean_repl_retained:.1f} | {row.mean_runtime:.1f} |"
        )
    (outdir / "aggregate_report.md").write_text("\n".join(lines) + "\n")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate and run SGTree contamination benchmarks")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gen = subparsers.add_parser("generate", help="Generate a benchmark dataset")
    gen.add_argument("--source-dir", default="testgenomes/Chloroflexi")
    gen.add_argument("--models", default="resources/models/UNI56.hmm")
    gen.add_argument("--outdir", default="runs/benchmarks/dev_chloroflexi")
    gen.add_argument("--n-genomes", type=int, default=20)
    gen.add_argument("--n-markers", type=int, default=10)
    gen.add_argument("--min-marker-presence-fraction", type=float, default=1.0)
    gen.add_argument("--seed", type=int, default=42)
    gen.add_argument("--num-cpus", type=int, default=8)

    prep = subparsers.add_parser("prepare-source-subset", help="Prepare a sampled source subset directory")
    prep.add_argument("--source-dir", required=True)
    prep.add_argument("--outdir", required=True)
    prep.add_argument("--list-file", default=None)
    prep.add_argument("--n-candidates", type=int, default=None)
    prep.add_argument("--seed", type=int, default=42)

    run = subparsers.add_parser("run", help="Run legacy vs coordinate benchmark comparisons")
    run.add_argument("--benchmark-dir", default="runs/benchmarks/dev_chloroflexi")
    run.add_argument("--num-cpus", type=int, default=8)
    run.add_argument("--selection-max-rounds", type=int, default=5)

    genrep = subparsers.add_parser("generate-replicates", help="Generate multiple benchmark replicates")
    genrep.add_argument("--source-dir", required=True)
    genrep.add_argument("--outbase", required=True)
    genrep.add_argument("--replicate-seeds", nargs="+", type=int, required=True)
    genrep.add_argument("--n-genomes", type=int, default=50)
    genrep.add_argument("--n-markers", type=int, default=15)
    genrep.add_argument("--min-marker-presence-fraction", type=float, default=1.0)
    genrep.add_argument("--num-cpus", type=int, default=8)
    genrep.add_argument("--list-file", default=None)
    genrep.add_argument("--n-candidates", type=int, default=90)

    runrep = subparsers.add_parser("run-replicates", help="Run benchmark suite for multiple replicate dirs")
    runrep.add_argument("--benchmark-dirs", nargs="+", required=True)
    runrep.add_argument("--num-cpus", type=int, default=8)
    runrep.add_argument("--selection-max-rounds", type=int, default=5)

    cross = subparsers.add_parser("generate-cross", help="Generate a cross-family benchmark from two benchmark dirs")
    cross.add_argument("--benchmark-a-dir", required=True)
    cross.add_argument("--benchmark-b-dir", required=True)
    cross.add_argument("--outdir", required=True)
    cross.add_argument("--n-markers", type=int, default=8)
    cross.add_argument("--seed", type=int, default=42)
    cross.add_argument("--num-cpus", type=int, default=8)

    agg = subparsers.add_parser("aggregate-replicates", help="Aggregate multiple benchmark result directories")
    agg.add_argument("--benchmark-dirs", nargs="+", required=True)
    agg.add_argument("--outdir", required=True)

    return parser.parse_args()


def prepare_source_subset(
    source_dir: Path,
    outdir: Path,
    *,
    list_file: Path | None,
    n_candidates: int | None,
    seed: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rng = Random(seed)
    candidates = []
    if list_file is not None:
        ids = [line.strip() for line in list_file.read_text().splitlines() if line.strip()]
        for genome_id in ids:
            prefixed = source_dir / f"FLAV__{genome_id}.faa"
            if prefixed.exists():
                candidates.append(prefixed)
    else:
        candidates = sorted(source_dir.glob("*.faa"))
    if n_candidates is not None and len(candidates) > n_candidates:
        candidates = sorted(rng.sample(candidates, n_candidates), key=lambda path: path.name)
    for file_path in candidates:
        target = outdir / file_path.name
        if target.exists():
            target.unlink()
        target.symlink_to(file_path.resolve())


def generate_replicates(
    source_dir: Path,
    outbase: Path,
    replicate_seeds: list[int],
    *,
    n_genomes: int,
    n_markers: int,
    min_marker_presence_fraction: float,
    num_cpus: int,
    list_file: Path | None,
    n_candidates: int | None,
) -> list[Path]:
    replicate_dirs = []
    for seed in replicate_seeds:
        source_subset = outbase / f"source_seed_{seed}"
        prepare_source_subset(
            source_dir=source_dir,
            outdir=source_subset,
            list_file=list_file,
            n_candidates=n_candidates,
            seed=seed,
        )
        replicate_dir = outbase / f"benchmark_seed_{seed}"
        generate_benchmark_dataset(
            source_dir=source_subset,
            models_path=Path("resources/models/UNI56.hmm"),
            outdir=replicate_dir,
            n_genomes=n_genomes,
            n_markers=n_markers,
            min_marker_presence_fraction=min_marker_presence_fraction,
            seed=seed,
            num_cpus=num_cpus,
        )
        replicate_dirs.append(replicate_dir)
    return replicate_dirs


def main() -> None:
    args = _parse_args()
    if args.command == "generate":
        generate_benchmark_dataset(
            source_dir=Path(args.source_dir),
            models_path=Path(args.models),
            outdir=Path(args.outdir),
            n_genomes=args.n_genomes,
            n_markers=args.n_markers,
            min_marker_presence_fraction=args.min_marker_presence_fraction,
            seed=args.seed,
            num_cpus=args.num_cpus,
        )
    elif args.command == "prepare-source-subset":
        prepare_source_subset(
            source_dir=Path(args.source_dir),
            outdir=Path(args.outdir),
            list_file=Path(args.list_file) if args.list_file else None,
            n_candidates=args.n_candidates,
            seed=args.seed,
        )
    elif args.command == "generate-replicates":
        generate_replicates(
            source_dir=Path(args.source_dir),
            outbase=Path(args.outbase),
            replicate_seeds=args.replicate_seeds,
            n_genomes=args.n_genomes,
            n_markers=args.n_markers,
            min_marker_presence_fraction=args.min_marker_presence_fraction,
            num_cpus=args.num_cpus,
            list_file=Path(args.list_file) if args.list_file else None,
            n_candidates=args.n_candidates,
        )
    elif args.command == "run":
        run_benchmark_suite(
            benchmark_dir=Path(args.benchmark_dir),
            num_cpus=args.num_cpus,
            selection_max_rounds=args.selection_max_rounds,
        )
    elif args.command == "generate-cross":
        generate_cross_benchmark_dataset(
            benchmark_a_dir=Path(args.benchmark_a_dir),
            benchmark_b_dir=Path(args.benchmark_b_dir),
            outdir=Path(args.outdir),
            n_markers=args.n_markers,
            seed=args.seed,
            num_cpus=args.num_cpus,
        )
    elif args.command == "run-replicates":
        for benchmark_dir in [Path(item) for item in args.benchmark_dirs]:
            run_benchmark_suite(
                benchmark_dir=benchmark_dir,
                num_cpus=args.num_cpus,
                selection_max_rounds=args.selection_max_rounds,
            )
    elif args.command == "aggregate-replicates":
        aggregate_benchmark_runs(
            benchmark_dirs=[Path(item) for item in args.benchmark_dirs],
            outdir=Path(args.outdir),
        )


if __name__ == "__main__":
    main()
