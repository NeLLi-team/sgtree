"""Synthetic benchmark generation and evaluation for SGTree."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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

from sgtree.benchmark_dataset import prepare_burkholderiaceae_benchmark_dataset
from sgtree.id_schema import build_sequence_id, parse_sequence_id, sanitize_token


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

SCENARIO_DESCRIPTION = {
    "duplicate_only": "Added-copy contamination only; no singleton replacement events.",
    "replacement_only": "Singleton replacement contamination only; no added duplicate copies.",
    "combined": "Added-copy and singleton-replacement contamination in the same panel.",
    "mixed_high_level": "Class/phylum-level mixed contamination with duplicate and replacement recipients plus explicit overlap.",
}

DEFAULT_CLEANUP_PROFILES = {
    "duplicate_only": {
        "name": "duplicate_cleanup",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": False,
        "singles_mode": "delta_rf",
    },
    "replacement_only": {
        "name": "singles_delta_rf",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "delta_rf",
    },
    "combined": {
        "name": "duplicate_plus_singles_delta_rf",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "delta_rf",
    },
    "mixed_high_level": {
        "name": "duplicate_plus_singles_delta_rf",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "delta_rf",
    },
}

DEFAULT_TAXONOMY_DB = Path("/home/fschulz/dev/nelli-genomes-db/resources/database/gtdb_genomes.duckdb")
PROJECT_ROOT = Path(__file__).resolve().parents[2]
ASSEMBLY_ACCESSION_RE = re.compile(r"(GC[AF])[-_]?(\d{9})[-._]?(\d+)")
TAXONOMY_SCOPE_RULES = {
    "genus": {
        "same_rank": "family",
        "different_rank": "genus",
        "scope_label": "same_family_different_genus",
    },
    "family": {
        "same_rank": "order_name",
        "different_rank": "family",
        "scope_label": "same_order_different_family",
    },
    "order": {
        "same_rank": "class",
        "different_rank": "order_name",
        "scope_label": "same_class_different_order",
    },
    "phylum": {
        "same_rank": "domain",
        "different_rank": "phylum",
        "scope_label": "same_domain_different_phylum",
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
    contaminant_gene = f"contam__{marker}__{donor_genome}__e{event_index:03d}"
    contaminant_contig = f"contig__contam__{marker}__{donor_genome}__e{event_index:03d}"
    record.id = build_sequence_id(
        recipient_genome,
        contaminant_contig,
        contaminant_gene,
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


def drop_native_marker(
    recipient_records: dict[str, SeqRecord],
    native_record_id: str,
) -> dict[str, SeqRecord]:
    updated = dict(recipient_records)
    if native_record_id not in updated:
        raise KeyError(f"Native marker record not found: {native_record_id}")
    del updated[native_record_id]
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


def _normalize_assembly_accession(token: str) -> str:
    match = ASSEMBLY_ACCESSION_RE.search(token)
    if not match:
        raise ValueError(f"Could not extract assembly accession from: {token}")
    return f"{match.group(1)}_{match.group(2)}.{match.group(3)}"


def _discover_input_files(source_dir: Path) -> list[Path]:
    exts = ("*.faa", "*.fna", "*.fa", "*.fasta")
    seen: list[Path] = []
    for pattern in exts:
        seen.extend(sorted(source_dir.glob(pattern)))
    return sorted({path for path in seen if path.is_file()}, key=lambda path: path.name)


def _source_genome_index(source_dir: Path) -> pd.DataFrame:
    rows = []
    for path in _discover_input_files(source_dir):
        rows.append(
            {
                "file_path": str(path.resolve()),
                "source_file": path.name,
                "genome_id": sanitize_token(path.stem, path.stem),
                "assembly_accession": _normalize_assembly_accession(path.stem),
            }
        )
    if not rows:
        raise FileNotFoundError(f"No input genome files found under {source_dir}")
    return pd.DataFrame(rows).sort_values(["genome_id"]).reset_index(drop=True)


def _load_taxonomy_rows(
    assembly_accessions: list[str],
    taxonomy_db_path: Path,
) -> pd.DataFrame:
    if not taxonomy_db_path.exists():
        raise FileNotFoundError(taxonomy_db_path)
    try:
        import duckdb
    except ImportError as exc:
        raise RuntimeError(
            "Taxonomic benchmark generation requires duckdb in the SGTree environment"
        ) from exc

    request_df = pd.DataFrame({"assembly_accession": sorted(set(assembly_accessions))})
    with duckdb.connect(str(taxonomy_db_path), read_only=True) as con:
        con.register("requested_accessions", request_df)
        rows = con.execute(
            """
            with combined as (
                select
                    assembly_accession,
                    organism_name,
                    phylum,
                    class,
                    order_name,
                    family,
                    genus,
                    species,
                    'gtdb' as taxonomy_source
                from gtdb_genomes
                union all
                select
                    assembly_accession,
                    organism_name,
                    phylum,
                    class,
                    order_name,
                    family,
                    genus,
                    species,
                    'non_gtdb' as taxonomy_source
                from non_gtdb_genomes
            ),
            ranked as (
                select
                    combined.*,
                    row_number() over (
                        partition by combined.assembly_accession
                        order by case when taxonomy_source = 'gtdb' then 0 else 1 end
                    ) as rn
                from combined
                inner join requested_accessions using (assembly_accession)
            )
            select
                assembly_accession,
                organism_name,
                phylum,
                class,
                order_name,
                family,
                genus,
                species,
                taxonomy_source
            from ranked
            where rn = 1
            """
        ).fetchdf()
    return rows


def _taxonomy_string(meta: dict[str, str]) -> str:
    return ";".join(
        [
            f"p__{meta.get('phylum', '')}",
            f"c__{meta.get('class', '')}",
            f"o__{meta.get('order_name', '')}",
            f"f__{meta.get('family', '')}",
            f"g__{meta.get('genus', '')}",
            f"s__{meta.get('species', '')}",
        ]
    )


def _load_source_taxonomy(
    source_dir: Path,
    taxonomy_db_path: Path,
) -> pd.DataFrame:
    index_df = _source_genome_index(source_dir)
    taxonomy_df = _load_taxonomy_rows(index_df["assembly_accession"].tolist(), taxonomy_db_path)
    merged = index_df.merge(taxonomy_df, on="assembly_accession", how="left", validate="one_to_one")
    missing = merged[merged["family"].isna() | merged["genus"].isna()]
    if not missing.empty:
        raise ValueError(
            "Missing taxonomy for one or more genomes: "
            + ", ".join(sorted(missing["source_file"].astype(str).tolist()))
        )
    merged["taxonomy"] = merged.apply(
        lambda row: _taxonomy_string(row.to_dict()),
        axis=1,
    )
    return merged.sort_values(["genome_id"]).reset_index(drop=True)


def _taxonomy_lookup(df: pd.DataFrame) -> dict[str, dict[str, str]]:
    return {
        str(row["genome_id"]): {
            "assembly_accession": str(row["assembly_accession"]),
            "organism_name": str(row["organism_name"]),
            "domain": "Bacteria",
            "phylum": str(row["phylum"]),
            "class": str(row["class"]),
            "order_name": str(row["order_name"]),
            "family": str(row["family"]),
            "genus": str(row["genus"]),
            "species": str(row["species"]),
            "taxonomy_source": str(row["taxonomy_source"]),
            "taxonomy": str(row["taxonomy"]),
            "source_file": str(row["source_file"]),
            "file_path": str(row["file_path"]),
        }
        for row in df.to_dict("records")
    }


def _taxonomy_field(meta: dict[str, str], key: str) -> str:
    return str(meta.get(key, "")).strip()


def _taxonomy_scope_matches(
    recipient_meta: dict[str, str],
    donor_meta: dict[str, str],
    scope: str,
) -> bool:
    if scope not in TAXONOMY_SCOPE_RULES:
        raise KeyError(f"Unsupported taxonomy scope: {scope}")
    rule = TAXONOMY_SCOPE_RULES[scope]
    same_value = _taxonomy_field(recipient_meta, rule["same_rank"])
    donor_same_value = _taxonomy_field(donor_meta, rule["same_rank"])
    different_value = _taxonomy_field(recipient_meta, rule["different_rank"])
    donor_different_value = _taxonomy_field(donor_meta, rule["different_rank"])
    if not same_value or not donor_same_value or not different_value or not donor_different_value:
        return False
    return same_value == donor_same_value and different_value != donor_different_value


def _taxonomic_donor_candidates(
    recipient_genome: str,
    marker: str,
    scope: str,
    recipient_taxonomy: dict[str, dict[str, str]],
    donor_taxonomy: dict[str, dict[str, str]],
    donor_native_map: dict[str, dict[str, str]],
    truth_tree: Tree | None = None,
) -> list[str]:
    recipient_meta = recipient_taxonomy[recipient_genome]
    candidates = [
        donor_genome
        for donor_genome, donor_meta in donor_taxonomy.items()
        if donor_genome != recipient_genome
        and marker in donor_native_map.get(donor_genome, {})
        and _taxonomy_scope_matches(recipient_meta, donor_meta, scope)
    ]
    if truth_tree is None:
        return sorted(candidates)
    truth_leaves = {leaf.name for leaf in truth_tree.iter_leaves()}
    if recipient_genome not in truth_leaves:
        return sorted(candidates)
    in_tree = [genome for genome in candidates if genome in truth_leaves]
    out_of_tree = sorted(genome for genome in candidates if genome not in set(in_tree))
    ranked_in_tree = _distance_ranked_donors(truth_tree, recipient_genome, in_tree) if in_tree else []
    return ranked_in_tree + out_of_tree


def _stage_source_dir(
    source_dir: Path,
    models_path: Path,
    stage_dir: Path,
    *,
    num_cpus: int,
) -> tuple[dict[str, dict[str, SeqRecord]], pd.DataFrame]:
    if not _stage_dir_complete(stage_dir):
        _run_sgtree_python(
            genomedir=source_dir,
            modeldir=models_path,
            outdir=stage_dir,
            num_cpus=num_cpus,
            percent_models=0,
            marker_selection=False,
            singles=False,
            singles_mode="delta_rf",
            selection_mode="coordinate",
            selection_max_rounds=5,
            selection_global_rounds=1,
            keep_intermediates=True,
        )
    return _read_normalized_proteomes(stage_dir / "proteomes"), _load_table(stage_dir / "table_elim_dups")


def _duplicate_event_target(spec: dict) -> int:
    return int(spec["pair_blocks"]) * int(spec["markers_per_block"]) * 2


def _event_taxonomy_fields(prefix: str, meta: dict[str, str]) -> dict[str, str]:
    return {
        f"{prefix}_assembly_accession": str(meta.get("assembly_accession", "")),
        f"{prefix}_organism_name": str(meta.get("organism_name", "")),
        f"{prefix}_phylum": str(meta.get("phylum", "")),
        f"{prefix}_class": str(meta.get("class", "")),
        f"{prefix}_order": str(meta.get("order_name", "")),
        f"{prefix}_family": str(meta.get("family", "")),
        f"{prefix}_genus": str(meta.get("genus", "")),
        f"{prefix}_species": str(meta.get("species", "")),
        f"{prefix}_taxonomy": str(meta.get("taxonomy", "")),
    }


def _choose_recipient_sets(
    genomes: list[str],
    *,
    duplicate_recipients: int,
    replacement_recipients: int,
    overlap_recipients: int,
    rng: Random,
) -> tuple[list[str], list[str]]:
    if overlap_recipients > duplicate_recipients or overlap_recipients > replacement_recipients:
        raise ValueError("overlap_recipients cannot exceed duplicate or replacement recipient counts")
    total_unique = duplicate_recipients + replacement_recipients - overlap_recipients
    if total_unique > len(genomes):
        raise ValueError("Not enough genomes to satisfy duplicate/replacement overlap layout")
    ordered = sorted(genomes)
    overlap = sorted(rng.sample(ordered, overlap_recipients))
    remaining = [genome for genome in ordered if genome not in set(overlap)]
    duplicate_only_count = duplicate_recipients - overlap_recipients
    replacement_only_count = replacement_recipients - overlap_recipients
    duplicate_only = sorted(rng.sample(remaining, duplicate_only_count))
    remaining = [genome for genome in remaining if genome not in set(duplicate_only)]
    replacement_only = sorted(rng.sample(remaining, replacement_only_count))
    return sorted(overlap + duplicate_only), sorted(overlap + replacement_only)


def _round_robin_labels(labels: list[str], n_events: int) -> list[str]:
    if not labels:
        raise ValueError("At least one donor label is required")
    return [labels[idx % len(labels)] for idx in range(n_events)]


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
    native_map: dict[str, dict[str, str]],
    n_needed: int,
    rng: Random,
) -> list[str]:
    available = [
        marker
        for marker in markers
        if all((genome, marker) not in used_pairs for genome in genomes)
        and all(marker in native_map.get(genome, {}) for genome in genomes)
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


def _expected_reference_taxa(
    manifest: dict,
    scenario_meta: dict,
    truth_tree_path: Path,
) -> list[str]:
    scenario_taxa = scenario_meta.get("reference_taxa")
    if scenario_taxa:
        return sorted(str(name) for name in scenario_taxa)
    selected = manifest.get("selected_genomes")
    if selected:
        return sorted(str(name) for name in selected)
    return _sorted_leaf_names(truth_tree_path)


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
    if "chlam" in joined:
        return "chlam"
    return label or "all"


def _stage_dir_complete(stage_dir: Path) -> bool:
    required = [
        stage_dir / "proteomes",
        stage_dir / "table_elim_dups",
        stage_dir / "tree.nwk",
    ]
    return all(path.exists() for path in required)


def _resolve_project_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def _link_existing_path(source: Path, dest: Path) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    if dest.is_symlink() or dest.exists():
        if dest.is_dir() and not dest.is_symlink():
            shutil.rmtree(dest)
        else:
            dest.unlink()
    dest.symlink_to(source.resolve())


def _taxonomy_map_from_manifest_rows(rows: list[dict]) -> dict[str, dict[str, str]]:
    taxonomy: dict[str, dict[str, str]] = {}
    for row in rows:
        genome_id = str(row.get("genome_id", "")).strip()
        if not genome_id:
            continue
        taxonomy[genome_id] = {
            str(key): "" if value is None else str(value)
            for key, value in row.items()
        }
    return taxonomy


def _load_staged_records_and_native_map(
    stage_dir: Path,
    truth_markers: list[str],
) -> tuple[dict[str, dict[str, SeqRecord]], dict[str, dict[str, str]]]:
    if not _stage_dir_complete(stage_dir):
        raise FileNotFoundError(f"Missing staged benchmark inputs: {stage_dir}")
    records = _read_normalized_proteomes(stage_dir / "proteomes")
    table = _load_table(stage_dir / "table_elim_dups")
    native_map = _native_marker_map(table, sorted(records), truth_markers)
    return records, native_map


def _prepare_replicate_output_dir(
    base_benchmark_dir: Path,
    outdir: Path,
) -> tuple[Path, Path, Path]:
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    truth_inputs_dir = outdir / "truth_inputs"
    truth_run_dir = outdir / "truth_run"
    truth_models = outdir / "truth_markers.hmm"
    _link_existing_path(base_benchmark_dir / "truth_inputs", truth_inputs_dir)
    _link_existing_path(base_benchmark_dir / "truth_run", truth_run_dir)
    _link_existing_path(base_benchmark_dir / "truth_markers.hmm", truth_models)
    scenarios_dir = outdir / "scenarios"
    scenarios_dir.mkdir(exist_ok=True)
    return truth_inputs_dir, truth_run_dir, scenarios_dir


def _write_reference_state(
    *,
    scenario_dir: Path,
    reference_records: dict[str, dict[str, SeqRecord]],
    truth_inputs_dir: Path,
    truth_run_dir: Path,
    truth_models: Path,
    num_cpus: int,
    reuse_truth_panel: bool,
) -> tuple[Path, Path]:
    reference_inputs_dir = scenario_dir / "reference_inputs"
    reference_run_dir = scenario_dir / "reference_run"
    if reuse_truth_panel:
        _link_existing_path(truth_inputs_dir, reference_inputs_dir)
        _link_existing_path(truth_run_dir, reference_run_dir)
        return reference_inputs_dir, reference_run_dir

    _write_proteome_dir(reference_records, reference_inputs_dir)
    _run_sgtree_python(
        genomedir=reference_inputs_dir,
        modeldir=truth_models,
        outdir=reference_run_dir,
        num_cpus=num_cpus,
        percent_models=0,
        marker_selection=False,
        singles=False,
        singles_mode="delta_rf",
        selection_mode="coordinate",
        selection_max_rounds=5,
        selection_global_rounds=1,
        keep_intermediates=True,
    )
    return reference_inputs_dir, reference_run_dir


def _parallel_benchmark_jobs(
    jobs: list[tuple[Path, callable]],
    *,
    max_parallel: int,
) -> None:
    if max_parallel <= 1:
        for _benchmark_dir, fn in jobs:
            fn()
        return

    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        future_to_dir = {
            executor.submit(fn): benchmark_dir
            for benchmark_dir, fn in jobs
        }
        for future in as_completed(future_to_dir):
            benchmark_dir = future_to_dir[future]
            try:
                future.result()
            except Exception as exc:
                raise RuntimeError(f"Benchmark job failed for {benchmark_dir}") from exc


def _donor_candidates_with_marker(
    candidates: list[str],
    marker: str,
    native_map: dict[str, dict[str, str]],
) -> list[str]:
    return [genome for genome in candidates if marker in native_map.get(genome, {})]


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
        # Keep the intended truth panel taxa even when a few chosen markers are
        # missing in some genomes; contamination events are filtered to loci that
        # exist for each recipient/donor pair.
        percent_models=0,
        marker_selection=False,
        singles=False,
        singles_mode="delta_rf",
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
        reference_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
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
                    native_map,
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
                                "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                                "donor_genome": donor,
                                "donor_group": selected_group_labels[donor],
                                "source_relation": _source_relation(selected_group_labels[recipient], selected_group_labels[donor]),
                                "donor_record_id": native_map[donor][marker],
                                "donor_contig_id": parse_sequence_id(native_map[donor][marker])[1],
                                "contaminant_record_id": contaminant.id,
                                "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
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
            and marker in native_map.get(recipient, {})
        ]
        rng.shuffle(replacement_pairs)
        for recipient, marker in replacement_pairs[: spec["replacement_events"]]:
            if cross_group_only:
                donor_group = "gamma" if selected_group_labels[recipient] == "flavo" else "flavo"
                donor_candidates = _donor_candidates_with_marker(
                    [
                        genome
                        for genome in selected_genomes
                        if selected_group_labels[genome] == donor_group
                    ],
                    marker,
                    native_map,
                )
            else:
                donor_candidates = _donor_candidates_with_marker(
                    _distance_ranked_donors(
                        truth_tree,
                        recipient,
                        [g for g in selected_genomes if g != recipient],
                    ),
                    marker,
                    native_map,
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
            reference_records[recipient] = drop_native_marker(
                reference_records[recipient],
                native_record_id=native_map[recipient][marker],
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
                    "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                    "donor_genome": donor,
                    "donor_group": selected_group_labels[donor],
                    "source_relation": _source_relation(selected_group_labels[recipient], selected_group_labels[donor]),
                    "donor_record_id": native_map[donor][marker],
                    "donor_contig_id": parse_sequence_id(native_map[donor][marker])[1],
                    "contaminant_record_id": contaminant.id,
                    "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                    "expected_replacement_outcome": "DropMarkerOrRemoveContaminant",
                    "native_degrade_fraction": spec["native_degrade_fraction"],
                }
            )
            event_index += 1

        scenario_dir = scenarios_dir / scenario_name
        proteomes_dir = scenario_dir / "proteomes"
        reference_inputs_dir = scenario_dir / "reference_inputs"
        reference_run_dir = scenario_dir / "reference_run"
        _write_proteome_dir(scenario_records, proteomes_dir)
        _write_proteome_dir(reference_records, reference_inputs_dir)
        _run_sgtree_python(
            genomedir=reference_inputs_dir,
            modeldir=truth_models,
            outdir=reference_run_dir,
            num_cpus=num_cpus,
            # The benchmark reference should preserve the original selected
            # taxa even when replacement events remove a few native markers.
            percent_models=0,
            marker_selection=False,
            singles=False,
            singles_mode="delta_rf",
            selection_mode="coordinate",
            selection_max_rounds=5,
            selection_global_rounds=1,
            keep_intermediates=True,
        )
        _write_events_tsv(scenario_dir / "events.tsv", events)
        genome_summary = _write_genome_summary_tsv(scenario_dir / "genome_summary.tsv", events)
        summary = _scenario_summary(events, genome_summary)
        manifest["scenarios"].append(
            {
                "name": scenario_name,
                "proteomes_dir": str(proteomes_dir),
                "reference_inputs_dir": str(reference_inputs_dir),
                "reference_tree_path": str(reference_run_dir / "tree.nwk"),
                "reference_taxa": selected_genomes,
                "events_path": str(scenario_dir / "events.tsv"),
                "genome_summary_path": str(scenario_dir / "genome_summary.tsv"),
                **summary,
            }
        )

    _write_manifest_json(outdir / "benchmark_manifest.json", manifest)


def _materialize_taxonomic_benchmark_from_truth(
    *,
    truth_records: dict[str, dict[str, SeqRecord]],
    truth_markers: list[str],
    truth_taxonomy: dict[str, dict[str, str]],
    donor_records: dict[str, dict[str, SeqRecord]],
    donor_native_map: dict[str, dict[str, str]],
    donor_taxonomy: dict[str, dict[str, str]],
    models_path: Path,
    outdir: Path,
    seed: int,
    num_cpus: int,
    lineage_label: str,
    donor_lineage_label: str,
    taxonomic_scope: str,
    manifest_extra: dict | None = None,
) -> None:
    if taxonomic_scope not in TAXONOMY_SCOPE_RULES:
        raise KeyError(f"Unsupported taxonomy scope: {taxonomic_scope}")

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
        percent_models=0,
        marker_selection=False,
        singles=False,
        singles_mode="delta_rf",
        selection_mode="coordinate",
        selection_max_rounds=5,
        selection_global_rounds=1,
        keep_intermediates=True,
    )

    truth_table = _load_table(truth_run / "table_elim_dups")
    native_map = _native_marker_map(truth_table, sorted(truth_records), truth_markers)
    truth_tree = _load_truth_tree(truth_run / "tree.nwk")
    selected_genomes = sorted(leaf.name for leaf in truth_tree.iter_leaves())
    truth_records = {genome: truth_records[genome] for genome in selected_genomes}
    native_map = {genome: native_map.get(genome, {}) for genome in selected_genomes}
    recipient_taxonomy = {genome: truth_taxonomy[genome] for genome in selected_genomes}

    manifest = {
        "seed": seed,
        "models_path": str(models_path),
        "selected_genomes": selected_genomes,
        "selected_markers": truth_markers,
        "taxonomic_scope": taxonomic_scope,
        "taxonomic_rule": TAXONOMY_SCOPE_RULES[taxonomic_scope],
        "lineage_label": lineage_label,
        "donor_lineage_label": donor_lineage_label,
        "selected_genome_taxonomy": [
            {
                "genome_id": genome,
                **recipient_taxonomy[genome],
            }
            for genome in selected_genomes
        ],
        "donor_pool_size": len(donor_records),
        "scenarios": [],
    }
    if manifest_extra:
        manifest.update(manifest_extra)

    scenarios_dir.mkdir(exist_ok=True)
    for scenario_name, spec in DEFAULT_SCENARIOS.items():
        scenario_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
        reference_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
        used_pairs: set[tuple[str, str]] = set()
        events: list[dict] = []
        event_index = 1

        duplicate_target = _duplicate_event_target(spec)
        duplicate_candidates = [
            (recipient, marker)
            for recipient in selected_genomes
            for marker in truth_markers
            if marker in native_map.get(recipient, {})
            and _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
        ]
        rng.shuffle(duplicate_candidates)
        for recipient, marker in duplicate_candidates:
            if len(events) >= duplicate_target:
                break
            if (recipient, marker) in used_pairs:
                continue
            donor_candidates = _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
            if not donor_candidates:
                continue
            donor = rng.choice(donor_candidates)
            donor_record = donor_records[donor][donor_native_map[donor][marker]]
            _degrade_record_in_place(
                scenario_records[recipient],
                native_record_id=native_map[recipient][marker],
                fraction=spec["native_degrade_fraction"],
                rng=rng,
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
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": TAXONOMY_SCOPE_RULES[taxonomic_scope]["scope_label"],
                    "recipient_genome": recipient,
                    "recipient_group": lineage_label,
                    "marker": marker,
                    "native_record_id": native_map[recipient][marker],
                    "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                    "donor_genome": donor,
                    "donor_group": donor_lineage_label,
                    "source_relation": _source_relation(lineage_label, donor_lineage_label),
                    "donor_record_id": donor_native_map[donor][marker],
                    "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                    "contaminant_record_id": contaminant.id,
                    "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                    "expected_duplicate_status": "Removed",
                    "expected_native_status": "Kept",
                    "native_degrade_fraction": spec["native_degrade_fraction"],
                    **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                    **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
                }
            )
            event_index += 1
        if len([row for row in events if row["event_type"] == "duplicate"]) != duplicate_target:
            raise ValueError(
                f"Could not schedule {duplicate_target} duplicate events for scope '{taxonomic_scope}'"
            )

        replacement_pairs = [
            (recipient, marker)
            for recipient in selected_genomes
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
            and marker in native_map.get(recipient, {})
            and _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
        ]
        rng.shuffle(replacement_pairs)
        replacement_events_written = 0
        for recipient, marker in replacement_pairs:
            if replacement_events_written >= spec["replacement_events"]:
                break
            donor_candidates = _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
            if not donor_candidates:
                continue
            donor = rng.choice(donor_candidates)
            donor_record = donor_records[donor][donor_native_map[donor][marker]]
            contaminant = make_contaminant_record(
                recipient_genome=recipient,
                donor_record=donor_record,
                marker=marker,
                donor_genome=donor,
                event_index=event_index,
            )
            reference_records[recipient] = drop_native_marker(
                reference_records[recipient],
                native_record_id=native_map[recipient][marker],
            )
            scenario_records[recipient] = apply_replacement_event(
                scenario_records[recipient],
                native_record_id=native_map[recipient][marker],
                contaminant_record=contaminant,
            )
            used_pairs.add((recipient, marker))
            events.append(
                {
                    "event_index": event_index,
                    "scenario": scenario_name,
                    "event_type": "replacement",
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": TAXONOMY_SCOPE_RULES[taxonomic_scope]["scope_label"],
                    "recipient_genome": recipient,
                    "recipient_group": lineage_label,
                    "marker": marker,
                    "native_record_id": native_map[recipient][marker],
                    "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                    "donor_genome": donor,
                    "donor_group": donor_lineage_label,
                    "source_relation": _source_relation(lineage_label, donor_lineage_label),
                    "donor_record_id": donor_native_map[donor][marker],
                    "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                    "contaminant_record_id": contaminant.id,
                    "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                    "expected_replacement_outcome": "DropMarkerOrRemoveContaminant",
                    "native_degrade_fraction": spec["native_degrade_fraction"],
                    **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                    **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
                }
            )
            replacement_events_written += 1
            event_index += 1
        if replacement_events_written != spec["replacement_events"]:
            raise ValueError(
                f"Could not schedule {spec['replacement_events']} replacement events for scope '{taxonomic_scope}'"
            )

        scenario_dir = scenarios_dir / scenario_name
        proteomes_dir = scenario_dir / "proteomes"
        reference_inputs_dir = scenario_dir / "reference_inputs"
        reference_run_dir = scenario_dir / "reference_run"
        _write_proteome_dir(scenario_records, proteomes_dir)
        _write_proteome_dir(reference_records, reference_inputs_dir)
        _run_sgtree_python(
            genomedir=reference_inputs_dir,
            modeldir=truth_models,
            outdir=reference_run_dir,
            num_cpus=num_cpus,
            percent_models=0,
            marker_selection=False,
            singles=False,
            singles_mode="delta_rf",
            selection_mode="coordinate",
            selection_max_rounds=5,
            selection_global_rounds=1,
            keep_intermediates=True,
        )
        _write_events_tsv(scenario_dir / "events.tsv", events)
        genome_summary = _write_genome_summary_tsv(scenario_dir / "genome_summary.tsv", events)
        summary = _scenario_summary(events, genome_summary)
        manifest["scenarios"].append(
            {
                "name": scenario_name,
                "taxonomic_scope": taxonomic_scope,
                "proteomes_dir": str(proteomes_dir),
                "reference_inputs_dir": str(reference_inputs_dir),
                "reference_tree_path": str(reference_run_dir / "tree.nwk"),
                "reference_taxa": selected_genomes,
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
            singles_mode="delta_rf",
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


def generate_taxonomic_benchmark_dataset(
    truth_source_dir: Path,
    donor_source_dir: Path | None,
    models_path: Path,
    outdir: Path,
    *,
    taxonomic_scope: str,
    taxonomy_db_path: Path,
    lineage_label: str | None,
    donor_lineage_label: str | None,
    n_genomes: int,
    n_markers: int,
    min_marker_presence_fraction: float,
    seed: int,
    num_cpus: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    truth_stage_dir = outdir / "stage_truth_full_model_clean"
    truth_normalized_records, truth_table = _stage_source_dir(
        truth_source_dir,
        models_path,
        truth_stage_dir,
        num_cpus=num_cpus,
    )
    truth_taxonomy_df = _load_source_taxonomy(truth_source_dir, taxonomy_db_path)
    truth_taxonomy_map = _taxonomy_lookup(truth_taxonomy_df)

    genome_sizes = {genome: len(records) for genome, records in truth_normalized_records.items()}
    selected_genomes, counts = _choose_genome_subset(truth_table, genome_sizes, n_genomes)
    truth_markers, ranking = _rank_markers(
        truth_table,
        truth_normalized_records,
        selected_genomes,
        counts,
        n_markers,
        min_marker_presence_fraction,
    )
    truth_records = {
        genome: deepcopy(truth_normalized_records[genome])
        for genome in selected_genomes
    }

    actual_donor_source_dir = donor_source_dir if donor_source_dir is not None else truth_source_dir
    if actual_donor_source_dir.resolve() == truth_source_dir.resolve():
        donor_normalized_records = truth_normalized_records
        donor_table = truth_table
        donor_taxonomy_map = truth_taxonomy_map
    else:
        donor_stage_dir = outdir / "stage_donor_full_model_clean"
        donor_normalized_records, donor_table = _stage_source_dir(
            actual_donor_source_dir,
            models_path,
            donor_stage_dir,
            num_cpus=num_cpus,
        )
        donor_taxonomy_df = _load_source_taxonomy(actual_donor_source_dir, taxonomy_db_path)
        donor_taxonomy_map = _taxonomy_lookup(donor_taxonomy_df)

    donor_native_map = _native_marker_map(
        donor_table,
        sorted(donor_normalized_records),
        truth_markers,
    )
    effective_lineage_label = lineage_label or _infer_group_label(truth_source_dir)
    effective_donor_lineage_label = donor_lineage_label or effective_lineage_label
    _materialize_taxonomic_benchmark_from_truth(
        truth_records=truth_records,
        truth_markers=truth_markers,
        truth_taxonomy={genome: truth_taxonomy_map[genome] for genome in selected_genomes},
        donor_records=donor_normalized_records,
        donor_native_map=donor_native_map,
        donor_taxonomy=donor_taxonomy_map,
        models_path=models_path,
        outdir=outdir,
        seed=seed,
        num_cpus=num_cpus,
        lineage_label=effective_lineage_label,
        donor_lineage_label=effective_donor_lineage_label,
        taxonomic_scope=taxonomic_scope,
        manifest_extra={
            "truth_source_dir": str(truth_source_dir),
            "donor_source_dir": str(actual_donor_source_dir),
            "taxonomy_db_path": str(taxonomy_db_path),
            "min_marker_presence_fraction": min_marker_presence_fraction,
            "marker_ranking": ranking.to_dict("records"),
        },
    )


def generate_taxonomic_benchmark_suite(
    truth_source_dir: Path,
    outbase: Path,
    *,
    models_path: Path,
    taxonomy_db_path: Path,
    lineage_label: str | None,
    scopes: list[str],
    donor_source_by_scope: dict[str, Path | None],
    donor_lineage_label_by_scope: dict[str, str | None],
    n_genomes: int,
    n_markers: int,
    min_marker_presence_fraction: float,
    seed: int,
    num_cpus: int,
) -> list[Path]:
    outbase.mkdir(parents=True, exist_ok=True)
    benchmark_dirs: list[Path] = []
    for scope in scopes:
        benchmark_dir = outbase / f"{scope}_scope"
        generate_taxonomic_benchmark_dataset(
            truth_source_dir=truth_source_dir,
            donor_source_dir=donor_source_by_scope.get(scope),
            models_path=models_path,
            outdir=benchmark_dir,
            taxonomic_scope=scope,
            taxonomy_db_path=taxonomy_db_path,
            lineage_label=lineage_label,
            donor_lineage_label=donor_lineage_label_by_scope.get(scope),
            n_genomes=n_genomes,
            n_markers=n_markers,
            min_marker_presence_fraction=min_marker_presence_fraction,
            seed=seed,
            num_cpus=num_cpus,
        )
        benchmark_dirs.append(benchmark_dir)
    _write_manifest_json(
        outbase / "suite_manifest.json",
        {
            "truth_source_dir": str(truth_source_dir),
            "models_path": str(models_path),
            "taxonomy_db_path": str(taxonomy_db_path),
            "lineage_label": lineage_label,
            "scopes": scopes,
            "benchmark_dirs": [str(path) for path in benchmark_dirs],
        },
    )
    return benchmark_dirs


def _materialize_mixed_lineage_benchmark_from_truth(
    *,
    truth_records: dict[str, dict[str, SeqRecord]],
    truth_markers: list[str],
    truth_taxonomy: dict[str, dict[str, str]],
    donor_bundles: dict[str, dict[str, object]],
    models_path: Path,
    outdir: Path,
    seed: int,
    num_cpus: int,
    lineage_label: str,
    duplicate_recipients: int,
    replacement_recipients: int,
    overlap_recipients: int,
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
        percent_models=0,
        marker_selection=False,
        singles=False,
        singles_mode="delta_rf",
        selection_mode="coordinate",
        selection_max_rounds=5,
        selection_global_rounds=1,
        keep_intermediates=True,
    )

    truth_table = _load_table(truth_run / "table_elim_dups")
    native_map = _native_marker_map(truth_table, sorted(truth_records), truth_markers)
    truth_tree = _load_truth_tree(truth_run / "tree.nwk")
    selected_genomes = sorted(leaf.name for leaf in truth_tree.iter_leaves())
    truth_records = {genome: truth_records[genome] for genome in selected_genomes}
    native_map = {genome: native_map.get(genome, {}) for genome in selected_genomes}
    recipient_taxonomy = {genome: truth_taxonomy[genome] for genome in selected_genomes}

    duplicate_recipients_list, replacement_recipients_list = _choose_recipient_sets(
        selected_genomes,
        duplicate_recipients=duplicate_recipients,
        replacement_recipients=replacement_recipients,
        overlap_recipients=overlap_recipients,
        rng=rng,
    )
    donor_labels = sorted(donor_bundles)
    duplicate_label_plan = _round_robin_labels(donor_labels, duplicate_recipients)
    replacement_label_plan = _round_robin_labels(list(reversed(donor_labels)), replacement_recipients)

    manifest = {
        "seed": seed,
        "models_path": str(models_path),
        "selected_genomes": selected_genomes,
        "selected_markers": truth_markers,
        "taxonomic_scope": "mixed_high_level",
        "lineage_label": lineage_label,
        "donor_labels": donor_labels,
        "duplicate_recipients": duplicate_recipients,
        "replacement_recipients": replacement_recipients,
        "overlap_recipients": overlap_recipients,
        "duplicate_recipient_genomes": duplicate_recipients_list,
        "replacement_recipient_genomes": replacement_recipients_list,
        "selected_genome_taxonomy": [
            {
                "genome_id": genome,
                **recipient_taxonomy[genome],
            }
            for genome in selected_genomes
        ],
        "scenarios": [],
    }
    if manifest_extra:
        manifest.update(manifest_extra)

    scenario_name = "mixed_high_level"
    scenario_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
    reference_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
    used_pairs: set[tuple[str, str]] = set()
    events: list[dict] = []
    event_index = 1

    for recipient, donor_label in zip(duplicate_recipients_list, duplicate_label_plan):
        donor_bundle = donor_bundles[donor_label]
        donor_records = donor_bundle["records"]
        donor_native_map = donor_bundle["native_map"]
        donor_taxonomy = donor_bundle["taxonomy"]
        marker_candidates = [
            marker
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
            and marker in native_map.get(recipient, {})
            and any(marker in donor_native_map.get(donor_genome, {}) for donor_genome in donor_records)
        ]
        if not marker_candidates:
            raise ValueError(f"No duplicate-event marker candidates for recipient {recipient} in mixed benchmark")
        rng.shuffle(marker_candidates)
        marker = marker_candidates[0]
        donor_candidates = [
            donor_genome
            for donor_genome in donor_records
            if marker in donor_native_map.get(donor_genome, {})
        ]
        donor = rng.choice(sorted(donor_candidates))
        donor_record = donor_records[donor][donor_native_map[donor][marker]]
        _degrade_record_in_place(
            scenario_records[recipient],
            native_record_id=native_map[recipient][marker],
            fraction=0.12,
            rng=rng,
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
                "taxonomic_scope": "mixed_high_level",
                "taxonomic_scope_label": f"{'+'.join(donor_labels)}_into_{lineage_label}",
                "recipient_genome": recipient,
                "recipient_group": lineage_label,
                "marker": marker,
                "native_record_id": native_map[recipient][marker],
                "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                "donor_genome": donor,
                "donor_group": donor_label,
                "source_relation": _source_relation(lineage_label, donor_label),
                "donor_record_id": donor_native_map[donor][marker],
                "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                "contaminant_record_id": contaminant.id,
                "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                "expected_duplicate_status": "Removed",
                "expected_native_status": "Kept",
                "native_degrade_fraction": 0.12,
                **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
            }
        )
        event_index += 1

    for recipient, donor_label in zip(replacement_recipients_list, replacement_label_plan):
        donor_bundle = donor_bundles[donor_label]
        donor_records = donor_bundle["records"]
        donor_native_map = donor_bundle["native_map"]
        donor_taxonomy = donor_bundle["taxonomy"]
        marker_candidates = [
            marker
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
            and marker in native_map.get(recipient, {})
            and any(marker in donor_native_map.get(donor_genome, {}) for donor_genome in donor_records)
        ]
        if not marker_candidates:
            raise ValueError(f"No replacement-event marker candidates for recipient {recipient} in mixed benchmark")
        rng.shuffle(marker_candidates)
        marker = marker_candidates[0]
        donor_candidates = [
            donor_genome
            for donor_genome in donor_records
            if marker in donor_native_map.get(donor_genome, {})
        ]
        donor = rng.choice(sorted(donor_candidates))
        donor_record = donor_records[donor][donor_native_map[donor][marker]]
        contaminant = make_contaminant_record(
            recipient_genome=recipient,
            donor_record=donor_record,
            marker=marker,
            donor_genome=donor,
            event_index=event_index,
        )
        reference_records[recipient] = drop_native_marker(
            reference_records[recipient],
            native_record_id=native_map[recipient][marker],
        )
        scenario_records[recipient] = apply_replacement_event(
            scenario_records[recipient],
            native_record_id=native_map[recipient][marker],
            contaminant_record=contaminant,
        )
        used_pairs.add((recipient, marker))
        events.append(
            {
                "event_index": event_index,
                "scenario": scenario_name,
                "event_type": "replacement",
                "taxonomic_scope": "mixed_high_level",
                "taxonomic_scope_label": f"{'+'.join(donor_labels)}_into_{lineage_label}",
                "recipient_genome": recipient,
                "recipient_group": lineage_label,
                "marker": marker,
                "native_record_id": native_map[recipient][marker],
                "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                "donor_genome": donor,
                "donor_group": donor_label,
                "source_relation": _source_relation(lineage_label, donor_label),
                "donor_record_id": donor_native_map[donor][marker],
                "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                "contaminant_record_id": contaminant.id,
                "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                "expected_replacement_outcome": "DropMarkerOrRemoveContaminant",
                "native_degrade_fraction": 0.12,
                **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
            }
        )
        event_index += 1

    scenario_dir = scenarios_dir / scenario_name
    proteomes_dir = scenario_dir / "proteomes"
    reference_inputs_dir = scenario_dir / "reference_inputs"
    reference_run_dir = scenario_dir / "reference_run"
    _write_proteome_dir(scenario_records, proteomes_dir)
    _write_proteome_dir(reference_records, reference_inputs_dir)
    _run_sgtree_python(
        genomedir=reference_inputs_dir,
        modeldir=truth_models,
        outdir=reference_run_dir,
        num_cpus=num_cpus,
        percent_models=0,
        marker_selection=False,
        singles=False,
        singles_mode="delta_rf",
        selection_mode="coordinate",
        selection_max_rounds=5,
        selection_global_rounds=1,
        keep_intermediates=True,
    )
    _write_events_tsv(scenario_dir / "events.tsv", events)
    genome_summary = _write_genome_summary_tsv(scenario_dir / "genome_summary.tsv", events)
    summary = _scenario_summary(events, genome_summary)
    manifest["scenarios"].append(
        {
            "name": scenario_name,
            "taxonomic_scope": "mixed_high_level",
            "proteomes_dir": str(proteomes_dir),
            "reference_inputs_dir": str(reference_inputs_dir),
            "reference_tree_path": str(reference_run_dir / "tree.nwk"),
            "reference_taxa": selected_genomes,
            "events_path": str(scenario_dir / "events.tsv"),
            "genome_summary_path": str(scenario_dir / "genome_summary.tsv"),
            **summary,
        }
    )
    _write_manifest_json(outdir / "benchmark_manifest.json", manifest)


def generate_mixed_lineage_benchmark_dataset(
    truth_source_dir: Path,
    donor_source_dirs: dict[str, Path],
    models_path: Path,
    outdir: Path,
    *,
    taxonomy_db_path: Path,
    lineage_label: str | None,
    n_genomes: int,
    n_markers: int,
    min_marker_presence_fraction: float,
    duplicate_recipients: int,
    replacement_recipients: int,
    overlap_recipients: int,
    seed: int,
    num_cpus: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    truth_stage_dir = outdir / "stage_truth_full_model_clean"
    truth_normalized_records, truth_table = _stage_source_dir(
        truth_source_dir,
        models_path,
        truth_stage_dir,
        num_cpus=num_cpus,
    )
    truth_taxonomy_df = _load_source_taxonomy(truth_source_dir, taxonomy_db_path)
    truth_taxonomy_map = _taxonomy_lookup(truth_taxonomy_df)
    genome_sizes = {genome: len(records) for genome, records in truth_normalized_records.items()}
    selected_genomes, counts = _choose_genome_subset(truth_table, genome_sizes, n_genomes)
    truth_markers, ranking = _rank_markers(
        truth_table,
        truth_normalized_records,
        selected_genomes,
        counts,
        n_markers,
        min_marker_presence_fraction,
    )
    truth_records = {
        genome: deepcopy(truth_normalized_records[genome])
        for genome in selected_genomes
    }

    donor_bundles: dict[str, dict[str, object]] = {}
    for donor_label, donor_source_dir in donor_source_dirs.items():
        donor_stage_dir = outdir / f"stage_donor_{sanitize_token(donor_label, donor_label)}"
        donor_records, donor_table = _stage_source_dir(
            donor_source_dir,
            models_path,
            donor_stage_dir,
            num_cpus=num_cpus,
        )
        donor_taxonomy_df = _load_source_taxonomy(donor_source_dir, taxonomy_db_path)
        donor_taxonomy = _taxonomy_lookup(donor_taxonomy_df)
        donor_bundles[donor_label] = {
            "records": donor_records,
            "native_map": _native_marker_map(donor_table, sorted(donor_records), truth_markers),
            "taxonomy": donor_taxonomy,
            "source_dir": str(donor_source_dir),
        }

    _materialize_mixed_lineage_benchmark_from_truth(
        truth_records=truth_records,
        truth_markers=truth_markers,
        truth_taxonomy={genome: truth_taxonomy_map[genome] for genome in selected_genomes},
        donor_bundles=donor_bundles,
        models_path=models_path,
        outdir=outdir,
        seed=seed,
        num_cpus=num_cpus,
        lineage_label=lineage_label or _infer_group_label(truth_source_dir),
        duplicate_recipients=duplicate_recipients,
        replacement_recipients=replacement_recipients,
        overlap_recipients=overlap_recipients,
        manifest_extra={
            "truth_source_dir": str(truth_source_dir),
            "donor_source_dirs": {label: str(path) for label, path in donor_source_dirs.items()},
            "taxonomy_db_path": str(taxonomy_db_path),
            "min_marker_presence_fraction": min_marker_presence_fraction,
            "marker_ranking": ranking.to_dict("records"),
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


def _resolve_manifest_path(path_value: str | Path, benchmark_dir: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    if path.exists():
        return path

    benchmark_relative = benchmark_dir / path
    if benchmark_relative.exists():
        return benchmark_relative

    resolved_benchmark_dir = benchmark_dir.resolve()
    if "runs" in resolved_benchmark_dir.parts:
        runs_index = resolved_benchmark_dir.parts.index("runs")
        root = Path(*resolved_benchmark_dir.parts[:runs_index])
        rooted_candidate = root / path
        if rooted_candidate.exists():
            return rooted_candidate

    return path


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


def _sorted_leaf_names(tree_path: Path) -> list[str]:
    tree = Tree(str(tree_path))
    return sorted(leaf.name for leaf in tree.iter_leaves())


def _format_taxa_list(names: list[str]) -> str:
    return ",".join(sorted(names))


def _marker_name_from_tree_path(path: Path) -> str:
    return path.stem.split("_")[-2]


def _singleton_pruning_pairs(run_dir: Path) -> list[tuple[str, str]]:
    no_duplicates_dir = run_dir / "protTrees" / "no_duplicates" / "out"
    no_singles_dir = run_dir / "protTrees" / "no_singles"
    if not no_duplicates_dir.exists() or not no_singles_dir.exists():
        return []

    removed_pairs: list[tuple[str, str]] = []
    for path in sorted(no_duplicates_dir.glob("*.nw")):
        candidate_path = no_singles_dir / path.name
        if not candidate_path.exists():
            continue
        before = {
            leaf.name.split("|")[0]
            for leaf in Tree(str(path)).iter_leaves()
        }
        after = {
            leaf.name.split("|")[0]
            for leaf in Tree(str(candidate_path)).iter_leaves()
        }
        marker = _marker_name_from_tree_path(path)
        for genome in sorted(before - after):
            removed_pairs.append((genome, marker))
    return removed_pairs


def _singleton_pruning_summary(run_dir: Path, replacement_events: pd.DataFrame) -> dict:
    replacement_pairs = {
        (str(row.recipient_genome), str(row.marker))
        for row in replacement_events.itertuples(index=False)
    }
    intended_removed: list[str] = []
    collateral_removed: list[str] = []
    collateral_genomes: set[str] = set()

    for genome, marker in _singleton_pruning_pairs(run_dir):
        label = f"{genome}:{marker}"
        if (genome, marker) in replacement_pairs:
            intended_removed.append(label)
        else:
            collateral_removed.append(label)
            collateral_genomes.add(genome)

    return {
        "singleton_pruned_total_count": len(intended_removed) + len(collateral_removed),
        "singleton_intended_removed_count": len(intended_removed),
        "singleton_intended_removed": _format_taxa_list(intended_removed),
        "singleton_collateral_removed_count": len(collateral_removed),
        "singleton_collateral_removed": _format_taxa_list(collateral_removed),
        "singleton_collateral_genome_count": len(collateral_genomes),
        "singleton_collateral_genomes": _format_taxa_list(sorted(collateral_genomes)),
    }


def evaluate_benchmark_run(
    benchmark_dir: Path,
    scenario_name: str,
    run_dir: Path,
    runtime_seconds: float,
) -> dict:
    manifest = json.loads((benchmark_dir / "benchmark_manifest.json").read_text())
    scenario_meta = next(
        (item for item in manifest["scenarios"] if item["name"] == scenario_name),
        None,
    )
    if scenario_meta is None:
        raise KeyError(f"Scenario not found in manifest: {scenario_name}")
    truth_tree = _resolve_manifest_path(
        scenario_meta.get("reference_tree_path", benchmark_dir / "truth_run" / "tree.nwk"),
        benchmark_dir=benchmark_dir,
    )
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
    reference_taxa = _expected_reference_taxa(manifest, scenario_meta, truth_tree)
    final_taxa = _sorted_leaf_names(result_tree)
    missing_taxa = sorted(set(reference_taxa) - set(final_taxa))
    extra_taxa = sorted(set(final_taxa) - set(reference_taxa))
    replacement_recipients = sorted(set(replacement_events["recipient_genome"])) if not replacement_events.empty else []
    replacement_recipient_losses = sorted(set(missing_taxa) & set(replacement_recipients))
    collateral_losses = sorted(set(missing_taxa) - set(replacement_recipients))
    singleton_summary = _singleton_pruning_summary(run_dir, replacement_events)

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
        "final_reference_taxa_count": len(reference_taxa),
        "final_observed_taxa_count": len(final_taxa),
        "final_taxa_match_reference": not missing_taxa and not extra_taxa,
        "final_missing_taxa_count": len(missing_taxa),
        "final_missing_taxa": _format_taxa_list(missing_taxa),
        "final_extra_taxa_count": len(extra_taxa),
        "final_extra_taxa": _format_taxa_list(extra_taxa),
        "replacement_recipient_genome_loss_count": len(replacement_recipient_losses),
        "replacement_recipient_genomes_lost": _format_taxa_list(replacement_recipient_losses),
        "collateral_genome_loss_count": len(collateral_losses),
        "collateral_genomes_lost": _format_taxa_list(collateral_losses),
        **singleton_summary,
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


def run_benchmark_suites(
    benchmark_dirs: list[Path],
    *,
    num_cpus: int,
    selection_max_rounds: int,
    max_parallel: int = 1,
) -> None:
    jobs = [
        (
            benchmark_dir,
            lambda benchmark_dir=benchmark_dir: run_benchmark_suite(
                benchmark_dir=benchmark_dir,
                num_cpus=num_cpus,
                selection_max_rounds=selection_max_rounds,
            ),
        )
        for benchmark_dir in benchmark_dirs
    ]
    _parallel_benchmark_jobs(jobs, max_parallel=max_parallel)


def aggregate_benchmark_runs(benchmark_dirs: list[Path], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    frames = []
    for benchmark_dir in benchmark_dirs:
        summary_path = benchmark_dir / "results" / "summary.tsv"
        if not summary_path.exists():
            continue
        df = pd.read_csv(summary_path, sep="\t")
        manifest_path = benchmark_dir / "benchmark_manifest.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            df["panel_id"] = _panel_id_from_manifest(manifest, benchmark_dir)
            df["panel_label"] = _panel_label_from_manifest(manifest, benchmark_dir)
        df["benchmark_dir"] = str(benchmark_dir)
        frames.append(df)
    if not frames:
        raise FileNotFoundError("No benchmark summary.tsv files found for aggregation")
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(outdir / "combined_summary.tsv", sep="\t", index=False)
    ok = combined[combined["status"] == "ok"].copy()
    group_cols = ["strategy", "scenario"]
    if "panel_id" in ok.columns:
        group_cols = ["panel_id", "panel_label"] + group_cols
    agg = (
        ok.groupby(group_cols, as_index=False)
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
        "| Panel | Strategy | Scenario | Runs | Initial RF | Final RF | RF Delta | Dup Removed | Native Retained | Replacement Dropped | Replacement Retained | Runtime (s) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in agg.sort_values(["strategy", "scenario"]).itertuples(index=False):
        panel_label = getattr(row, "panel_label", "all_panels")
        lines.append(
            f"| {panel_label} | {row.strategy} | {row.scenario} | {row.n_runs} | "
            f"{row.mean_initial_rf:.3f} | {row.mean_final_rf:.3f} | {row.mean_rf_delta:.3f} | "
            f"{row.mean_dup_removed:.1f}/{row.mean_dup_events:.1f} | "
            f"{row.mean_native_retained:.1f}/{row.mean_dup_events:.1f} | "
            f"{row.mean_repl_dropped:.1f} | {row.mean_repl_retained:.1f} | {row.mean_runtime:.1f} |"
        )
    (outdir / "aggregate_report.md").write_text("\n".join(lines) + "\n")


def _panel_id_from_manifest(manifest: dict, benchmark_dir: Path) -> str:
    lineage = str(manifest.get("lineage_label", "panel")).strip() or "panel"
    scope = str(manifest.get("taxonomic_scope", benchmark_dir.name)).strip() or benchmark_dir.name
    return f"{lineage}_{scope}"


def _panel_label_from_manifest(manifest: dict, benchmark_dir: Path) -> str:
    lineage = str(manifest.get("lineage_label", benchmark_dir.name)).strip() or benchmark_dir.name
    scope = str(manifest.get("taxonomic_scope", benchmark_dir.name)).strip() or benchmark_dir.name
    if scope == "mixed_high_level":
        return f"{lineage.capitalize()} Mixed High Level"
    return f"{lineage.capitalize()} {scope.capitalize()}"


def discover_taxonomic_benchmark_dirs(benchmarks_root: Path) -> list[Path]:
    best_by_panel: dict[str, tuple[bool, float, Path]] = {}
    for manifest_path in sorted(benchmarks_root.glob("**/benchmark_manifest.json")):
        benchmark_dir = manifest_path.parent
        try:
            manifest = json.loads(manifest_path.read_text())
        except json.JSONDecodeError:
            continue
        if manifest.get("replicate_source_benchmark_dir"):
            continue
        if "lineage_label" not in manifest or "taxonomic_scope" not in manifest:
            continue
        panel_id = _panel_id_from_manifest(manifest, benchmark_dir)
        summary_path = benchmark_dir / "results" / "summary.tsv"
        has_summary = summary_path.exists()
        mtime = summary_path.stat().st_mtime if has_summary else manifest_path.stat().st_mtime
        current = best_by_panel.get(panel_id)
        candidate = (has_summary, mtime, benchmark_dir)
        if current is None or candidate[:2] > current[:2]:
            best_by_panel[panel_id] = candidate
    return sorted(item[2] for item in best_by_panel.values())


def build_benchmark_exports(benchmark_dirs: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dataset_rows: list[dict] = []
    genome_frames: list[pd.DataFrame] = []
    result_frames: list[pd.DataFrame] = []
    for benchmark_dir in benchmark_dirs:
        manifest_path = benchmark_dir / "benchmark_manifest.json"
        if not manifest_path.exists():
            continue
        manifest = json.loads(manifest_path.read_text())
        panel_id = _panel_id_from_manifest(manifest, benchmark_dir)
        panel_label = _panel_label_from_manifest(manifest, benchmark_dir)
        lineage_label = str(manifest.get("lineage_label", ""))
        taxonomic_scope = str(manifest.get("taxonomic_scope", ""))
        taxonomic_scope_label = str(manifest.get("taxonomic_rule", {}).get("scope_label", taxonomic_scope))
        donor_labels = manifest.get("donor_labels")
        if donor_labels is None:
            donor_lineage_label = str(manifest.get("donor_lineage_label", lineage_label))
            donor_labels = [donor_lineage_label] if donor_lineage_label else []
        donor_labels_text = ",".join(str(label) for label in donor_labels)
        donor_source_dirs = manifest.get("donor_source_dirs", {})
        donor_source_text = (
            ",".join(f"{label}:{path}" for label, path in sorted(donor_source_dirs.items()))
            if donor_source_dirs
            else str(manifest.get("donor_source_dir", ""))
        )

        scenario_meta = {
            item["name"]: item
            for item in manifest.get("scenarios", [])
        }
        for scenario_name, scenario in scenario_meta.items():
            dataset_rows.append(
                {
                    "panel_id": panel_id,
                    "panel_label": panel_label,
                    "dataset_id": panel_id,
                    "dataset_label": panel_label,
                    "lineage_label": lineage_label,
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": taxonomic_scope_label,
                    "scenario": scenario_name,
                    "scenario_description": SCENARIO_DESCRIPTION.get(scenario_name, scenario_name),
                    "selected_genome_count": len(manifest.get("selected_genomes", [])),
                    "selected_marker_count": len(manifest.get("selected_markers", [])),
                    "contaminant_markers_added": scenario.get("contaminant_markers_added", 0),
                    "duplicate_events": scenario.get("duplicate_events", 0),
                    "replacement_events": scenario.get("replacement_events", 0),
                    "affected_genomes": scenario.get("affected_genomes", 0),
                    "mean_contaminants_per_affected_genome": scenario.get("mean_contaminants_per_affected_genome", 0.0),
                    "max_contaminants_per_genome": scenario.get("max_contaminants_per_genome", 0),
                    "within_group_events": scenario.get("within_group_events", 0),
                    "cross_group_events": scenario.get("cross_group_events", 0),
                    "truth_tree_path": scenario.get("reference_tree_path", ""),
                    "reference_tree_path": scenario.get("reference_tree_path", ""),
                    "events_path": scenario.get("events_path", ""),
                    "genome_summary_path": scenario.get("genome_summary_path", ""),
                    "donor_lineage_label": str(manifest.get("donor_lineage_label", lineage_label)),
                    "donor_labels": donor_labels_text,
                    "truth_source_dir": str(manifest.get("truth_source_dir", manifest.get("source_dir", ""))),
                    "donor_source_dirs": donor_source_text,
                    "duplicate_recipients": manifest.get("duplicate_recipients"),
                    "replacement_recipients": manifest.get("replacement_recipients"),
                    "overlap_recipients": manifest.get("overlap_recipients"),
                }
            )

            genome_summary_value = scenario.get("genome_summary_path")
            if genome_summary_value:
                genome_summary_path = Path(genome_summary_value)
                if not genome_summary_path.is_absolute():
                    genome_summary_path = PROJECT_ROOT / genome_summary_path
                if genome_summary_path.exists():
                    genome_df = pd.read_csv(genome_summary_path, sep="\t")
                    genome_df["panel_id"] = panel_id
                    genome_df["panel_label"] = panel_label
                    genome_df["dataset_id"] = panel_id
                    genome_df["dataset_label"] = panel_label
                    genome_df["lineage_label"] = lineage_label
                    genome_df["taxonomic_scope"] = taxonomic_scope
                    genome_df["taxonomic_scope_label"] = taxonomic_scope_label
                    genome_df["scenario"] = scenario_name
                    genome_df["donor_labels"] = donor_labels_text
                    genome_frames.append(genome_df)

        summary_path = benchmark_dir / "results" / "summary.tsv"
        if summary_path.exists():
            result_df = pd.read_csv(summary_path, sep="\t")
            scenario_df = pd.DataFrame(dataset_rows)
            scenario_df = scenario_df[scenario_df["panel_id"] == panel_id][
                [
                    "panel_id",
                    "panel_label",
                    "dataset_id",
                    "dataset_label",
                    "lineage_label",
                    "taxonomic_scope",
                    "taxonomic_scope_label",
                    "scenario",
                    "donor_lineage_label",
                    "donor_labels",
                    "duplicate_recipients",
                    "replacement_recipients",
                    "overlap_recipients",
                ]
            ].drop_duplicates()
            result_df = result_df.merge(scenario_df, on="scenario", how="left", validate="many_to_one")
            result_frames.append(result_df)

    dataset_overview = pd.DataFrame(dataset_rows).sort_values(["lineage_label", "taxonomic_scope", "scenario"])
    genome_overview = pd.concat(genome_frames, ignore_index=True) if genome_frames else pd.DataFrame()
    benchmark_results = pd.concat(result_frames, ignore_index=True) if result_frames else pd.DataFrame()
    return dataset_overview, genome_overview, benchmark_results


def _latest_results_dir(benchmark_dir: Path, prefix: str) -> Path | None:
    candidates: list[tuple[float, str, Path]] = []
    for summary_path in benchmark_dir.glob(f"{prefix}*/summary.tsv"):
        if summary_path.is_file():
            results_dir = summary_path.parent
            candidates.append((summary_path.stat().st_mtime, results_dir.name, results_dir))
    if not candidates:
        return None
    return max(candidates)[2]


def _frame_series(frame: pd.DataFrame, column: str, default: object = pd.NA) -> pd.Series:
    if column in frame.columns:
        return frame[column]
    return pd.Series([default] * len(frame), index=frame.index)


def build_alignment_comparison_exports(benchmark_dirs: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    comparison_frames: list[pd.DataFrame] = []
    for benchmark_dir in benchmark_dirs:
        manifest_path = benchmark_dir / "benchmark_manifest.json"
        baseline_summary_path = benchmark_dir / "results" / "summary.tsv"
        comparison_results_dir = _latest_results_dir(benchmark_dir, "results_mafft_")
        if not manifest_path.exists() or not baseline_summary_path.exists() or comparison_results_dir is None:
            continue

        comparison_summary_path = comparison_results_dir / "summary.tsv"
        if not comparison_summary_path.exists():
            continue

        manifest = json.loads(manifest_path.read_text())
        panel_id = _panel_id_from_manifest(manifest, benchmark_dir)
        panel_label = _panel_label_from_manifest(manifest, benchmark_dir)
        lineage_label = str(manifest.get("lineage_label", ""))
        taxonomic_scope = str(manifest.get("taxonomic_scope", ""))
        taxonomic_scope_label = str(manifest.get("taxonomic_rule", {}).get("scope_label", taxonomic_scope))

        baseline_df = pd.read_csv(baseline_summary_path, sep="\t")
        comparison_df = pd.read_csv(comparison_summary_path, sep="\t")
        merged = baseline_df.merge(
            comparison_df,
            on="scenario",
            how="inner",
            suffixes=("_baseline", "_comparison"),
            validate="one_to_one",
        )
        if merged.empty:
            continue

        baseline_tree_rf = pd.to_numeric(_frame_series(merged, "tree_rf_norm_baseline"), errors="coerce")
        comparison_tree_rf = pd.to_numeric(_frame_series(merged, "tree_rf_norm_comparison"), errors="coerce")
        baseline_removed_fraction = pd.to_numeric(
            _frame_series(merged, "contaminant_markers_removed_fraction_baseline"),
            errors="coerce",
        )
        comparison_removed_fraction = pd.to_numeric(
            _frame_series(merged, "contaminant_markers_removed_fraction_comparison"),
            errors="coerce",
        )
        baseline_runtime = pd.to_numeric(_frame_series(merged, "runtime_seconds_baseline"), errors="coerce")
        comparison_runtime = pd.to_numeric(_frame_series(merged, "runtime_seconds_comparison"), errors="coerce")
        runtime_ratio = comparison_runtime / baseline_runtime.replace(0, pd.NA)
        rf_delta = comparison_tree_rf - baseline_tree_rf
        removed_delta = comparison_removed_fraction - baseline_removed_fraction
        rf_change_direction = pd.Series("unchanged", index=merged.index)
        rf_change_direction.loc[rf_delta > 0] = "worse"
        rf_change_direction.loc[rf_delta < 0] = "improved"

        alignment_method = _frame_series(merged, "alignment_method_comparison", "mafft").fillna("mafft")
        comparison_frames.append(
            pd.DataFrame(
                {
                    "panel_id": panel_id,
                    "panel_label": panel_label,
                    "benchmark_dir": str(benchmark_dir),
                    "lineage_label": lineage_label,
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": taxonomic_scope_label,
                    "scenario": merged["scenario"],
                    "baseline_results_tag": "results",
                    "comparison_results_tag": comparison_results_dir.name,
                    "comparison_alignment_method": alignment_method,
                    "baseline_run_dir": _frame_series(merged, "run_dir_baseline", ""),
                    "comparison_run_dir": _frame_series(merged, "run_dir_comparison", ""),
                    "baseline_status": _frame_series(merged, "status_baseline", ""),
                    "comparison_status": _frame_series(merged, "status_comparison", ""),
                    "baseline_tree_rf_norm": baseline_tree_rf,
                    "comparison_tree_rf_norm": comparison_tree_rf,
                    "rf_delta_comparison_minus_baseline": rf_delta,
                    "rf_change_direction": rf_change_direction,
                    "baseline_contaminant_markers_removed_fraction": baseline_removed_fraction,
                    "comparison_contaminant_markers_removed_fraction": comparison_removed_fraction,
                    "removed_fraction_delta_comparison_minus_baseline": removed_delta,
                    "baseline_runtime_seconds": baseline_runtime,
                    "comparison_runtime_seconds": comparison_runtime,
                    "runtime_ratio_comparison_over_baseline": runtime_ratio,
                    "baseline_missing_taxa_count": pd.to_numeric(
                        _frame_series(merged, "final_missing_taxa_count_baseline"),
                        errors="coerce",
                    ),
                    "comparison_missing_taxa_count": pd.to_numeric(
                        _frame_series(merged, "final_missing_taxa_count_comparison"),
                        errors="coerce",
                    ),
                }
            )
        )

    comparison_results = pd.concat(comparison_frames, ignore_index=True) if comparison_frames else pd.DataFrame()
    if comparison_results.empty:
        return comparison_results, pd.DataFrame()

    comparison_summary = (
        comparison_results.groupby(
            [
                "panel_id",
                "panel_label",
                "benchmark_dir",
                "lineage_label",
                "taxonomic_scope",
                "taxonomic_scope_label",
                "comparison_results_tag",
                "comparison_alignment_method",
            ],
            dropna=False,
        )
        .agg(
            scenario_count=("scenario", "count"),
            improved_count=("rf_change_direction", lambda s: int((s == "improved").sum())),
            unchanged_count=("rf_change_direction", lambda s: int((s == "unchanged").sum())),
            worsened_count=("rf_change_direction", lambda s: int((s == "worse").sum())),
            baseline_tree_rf_norm_mean=("baseline_tree_rf_norm", "mean"),
            comparison_tree_rf_norm_mean=("comparison_tree_rf_norm", "mean"),
            rf_delta_mean=("rf_delta_comparison_minus_baseline", "mean"),
            rf_delta_min=("rf_delta_comparison_minus_baseline", "min"),
            rf_delta_max=("rf_delta_comparison_minus_baseline", "max"),
            baseline_removed_fraction_mean=("baseline_contaminant_markers_removed_fraction", "mean"),
            comparison_removed_fraction_mean=("comparison_contaminant_markers_removed_fraction", "mean"),
            removed_fraction_delta_mean=("removed_fraction_delta_comparison_minus_baseline", "mean"),
            baseline_runtime_seconds_mean=("baseline_runtime_seconds", "mean"),
            comparison_runtime_seconds_mean=("comparison_runtime_seconds", "mean"),
            runtime_ratio_mean=("runtime_ratio_comparison_over_baseline", "mean"),
        )
        .reset_index()
    )
    return comparison_results, comparison_summary


def _parse_singleton_result_file(path: Path) -> list[dict]:
    rows: list[dict] = []
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = raw_line.strip()
        if not line.startswith("sin"):
            continue
        leaf_field, *fields = line.split("\t")
        leaf_name = leaf_field[3:]
        genome_id, contig_id, gene_id = parse_sequence_id(leaf_name)
        row = {
            "marker_name": path.name,
            "leaf_name": leaf_name,
            "genome": genome_id,
            "contig_id": contig_id,
            "gene_id": gene_id,
            "is_contaminant_leaf": "contam__" in leaf_name,
        }
        for field in fields:
            if "=" not in field:
                continue
            key, value = field.split("=", 1)
            row[key] = value
        rows.append(row)
    return rows


def build_singleton_classification_exports(
    benchmark_dirs: list[Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    call_rows: list[dict] = []
    summary_rows: list[dict] = []
    for benchmark_dir in benchmark_dirs:
        manifest_path = benchmark_dir / "benchmark_manifest.json"
        summary_path = benchmark_dir / "results" / "summary.tsv"
        if not manifest_path.exists() or not summary_path.exists():
            continue

        manifest = json.loads(manifest_path.read_text())
        panel_id = _panel_id_from_manifest(manifest, benchmark_dir)
        panel_label = _panel_label_from_manifest(manifest, benchmark_dir)
        lineage_label = str(manifest.get("lineage_label", ""))
        taxonomic_scope = str(manifest.get("taxonomic_scope", ""))
        taxonomic_scope_label = str(manifest.get("taxonomic_rule", {}).get("scope_label", taxonomic_scope))
        donor_labels = manifest.get("donor_labels")
        if donor_labels is None:
            donor_lineage_label = str(manifest.get("donor_lineage_label", lineage_label))
            donor_labels = [donor_lineage_label] if donor_lineage_label else []
        donor_labels_text = ",".join(str(label) for label in donor_labels)

        result_df = pd.read_csv(summary_path, sep="\t")
        for row in result_df.to_dict("records"):
            scenario_name = str(row.get("scenario", ""))
            replacement_events = int(row.get("replacement_events", 0) or 0)
            if replacement_events <= 0:
                continue

            run_dir_value = row.get("run_dir", "")
            if not run_dir_value:
                continue
            run_dir = Path(str(run_dir_value))
            if not run_dir.is_absolute():
                run_dir = PROJECT_ROOT / run_dir
            removed_dir = run_dir / "removed"
            scenario_calls: list[dict] = []
            if removed_dir.exists():
                for removed_path in sorted(path for path in removed_dir.iterdir() if path.is_file()):
                    for call in _parse_singleton_result_file(removed_path):
                        enriched = {
                            **call,
                            "panel_id": panel_id,
                            "panel_label": panel_label,
                            "benchmark_dir": str(benchmark_dir),
                            "lineage_label": lineage_label,
                            "taxonomic_scope": taxonomic_scope,
                            "taxonomic_scope_label": taxonomic_scope_label,
                            "scenario": scenario_name,
                            "cleanup_profile": str(row.get("cleanup_profile", "")),
                            "run_dir": str(run_dir),
                            "donor_labels": donor_labels_text,
                        }
                        scenario_calls.append(enriched)
                        call_rows.append(enriched)

            scenario_call_df = pd.DataFrame(scenario_calls)
            contamination_candidate_calls = (
                scenario_call_df["singleton_class"].eq("contamination_candidate").sum()
                if not scenario_call_df.empty and "singleton_class" in scenario_call_df.columns
                else 0
            )
            true_contamination_candidate_calls = (
                (
                    scenario_call_df["singleton_class"].eq("contamination_candidate")
                    & scenario_call_df["is_contaminant_leaf"]
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            false_contamination_candidate_calls = (
                (
                    scenario_call_df["singleton_class"].eq("contamination_candidate")
                    & ~scenario_call_df["is_contaminant_leaf"]
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            hgt_candidate_kept_calls = (
                (
                    scenario_call_df["singleton_class"].eq("hgt_candidate")
                    & scenario_call_df["decision"].eq("kept_hgt_candidate")
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            non_contaminant_proposals = (
                (~scenario_call_df["is_contaminant_leaf"]).sum()
                if not scenario_call_df.empty
                else 0
            )
            replacement_proposed_as_contamination = (
                (
                    scenario_call_df["is_contaminant_leaf"]
                    & scenario_call_df["singleton_class"].eq("contamination_candidate")
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            replacement_proposed_as_hgt = (
                (
                    scenario_call_df["is_contaminant_leaf"]
                    & scenario_call_df["singleton_class"].eq("hgt_candidate")
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            native_pruned_count = (
                (
                    ~scenario_call_df["is_contaminant_leaf"]
                    & scenario_call_df["decision"].eq("pruned")
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            native_kept_hgt_candidate_count = (
                (
                    ~scenario_call_df["is_contaminant_leaf"]
                    & scenario_call_df["decision"].eq("kept_hgt_candidate")
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            native_blocked_by_genome_budget_count = (
                (
                    ~scenario_call_df["is_contaminant_leaf"]
                    & scenario_call_df["decision"].eq("blocked_by_genome_budget")
                ).sum()
                if not scenario_call_df.empty
                else 0
            )
            replacement_not_proposed_count = max(
                0,
                replacement_events - int(scenario_call_df["is_contaminant_leaf"].sum() if not scenario_call_df.empty else 0),
            )
            precision = (
                true_contamination_candidate_calls / contamination_candidate_calls
                if contamination_candidate_calls
                else pd.NA
            )
            specificity = (
                hgt_candidate_kept_calls / non_contaminant_proposals
                if non_contaminant_proposals
                else pd.NA
            )
            recall = (
                float(row.get("replacement_contaminant_removed", 0)) / replacement_events
                if replacement_events
                else pd.NA
            )

            summary_rows.append(
                {
                    "panel_id": panel_id,
                    "panel_label": panel_label,
                    "benchmark_dir": str(benchmark_dir),
                    "lineage_label": lineage_label,
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": taxonomic_scope_label,
                    "scenario": scenario_name,
                    "cleanup_profile": str(row.get("cleanup_profile", "")),
                    "run_dir": str(run_dir),
                    "donor_labels": donor_labels_text,
                    "replacement_events": replacement_events,
                    "replacement_contaminant_removed": int(row.get("replacement_contaminant_removed", 0) or 0),
                    "replacement_contaminant_retained": int(row.get("replacement_contaminant_retained", 0) or 0),
                    "singleton_pruned_total_count": int(row.get("singleton_pruned_total_count", 0) or 0),
                    "singleton_intended_removed_count": int(row.get("singleton_intended_removed_count", 0) or 0),
                    "singleton_collateral_removed_count": int(row.get("singleton_collateral_removed_count", 0) or 0),
                    "proposal_total_count": int(len(scenario_call_df)),
                    "non_contaminant_proposal_count": int(non_contaminant_proposals),
                    "contamination_candidate_call_count": int(contamination_candidate_calls),
                    "true_contamination_candidate_call_count": int(true_contamination_candidate_calls),
                    "false_contamination_candidate_call_count": int(false_contamination_candidate_calls),
                    "hgt_candidate_kept_count": int(hgt_candidate_kept_calls),
                    "replacement_proposed_as_contamination_count": int(replacement_proposed_as_contamination),
                    "replacement_proposed_as_hgt_count": int(replacement_proposed_as_hgt),
                    "replacement_not_proposed_count": int(replacement_not_proposed_count),
                    "native_pruned_count": int(native_pruned_count),
                    "native_kept_hgt_candidate_count": int(native_kept_hgt_candidate_count),
                    "native_blocked_by_genome_budget_count": int(native_blocked_by_genome_budget_count),
                    "replacement_contaminant_recall": recall,
                    "contamination_candidate_precision": precision,
                    "non_contaminant_hgt_keep_fraction": specificity,
                }
            )

    call_df = pd.DataFrame(call_rows)
    summary_df = pd.DataFrame(summary_rows)
    return call_df, summary_df


def export_benchmark_tables(benchmark_dirs: list[Path], outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    dataset_overview, genome_overview, benchmark_results = build_benchmark_exports(benchmark_dirs)
    alignment_comparison, alignment_comparison_summary = build_alignment_comparison_exports(benchmark_dirs)
    singleton_calls, singleton_summary = build_singleton_classification_exports(benchmark_dirs)
    dataset_overview.to_csv(outdir / "benchmark_dataset_overview.tsv", sep="\t", index=False)
    if not genome_overview.empty:
        genome_overview.to_csv(outdir / "benchmark_genome_contamination.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(outdir / "benchmark_genome_contamination.tsv", sep="\t", index=False)
    if not benchmark_results.empty:
        benchmark_results.sort_values(["panel_label", "scenario"]).to_csv(outdir / "benchmark_summary_all.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(outdir / "benchmark_summary_all.tsv", sep="\t", index=False)
    if not alignment_comparison.empty:
        alignment_comparison.sort_values(["panel_label", "scenario"]).to_csv(
            outdir / "benchmark_alignment_comparison.tsv",
            sep="\t",
            index=False,
        )
    else:
        pd.DataFrame().to_csv(outdir / "benchmark_alignment_comparison.tsv", sep="\t", index=False)
    if not alignment_comparison_summary.empty:
        alignment_comparison_summary.sort_values(["panel_label"]).to_csv(
            outdir / "benchmark_alignment_comparison_summary.tsv",
            sep="\t",
            index=False,
        )
    else:
        pd.DataFrame().to_csv(outdir / "benchmark_alignment_comparison_summary.tsv", sep="\t", index=False)
    if not singleton_calls.empty:
        singleton_calls.sort_values(["panel_label", "scenario", "marker_name", "leaf_name"]).to_csv(
            outdir / "benchmark_singleton_classification.tsv",
            sep="\t",
            index=False,
        )
    else:
        pd.DataFrame().to_csv(outdir / "benchmark_singleton_classification.tsv", sep="\t", index=False)
    if not singleton_summary.empty:
        singleton_summary.sort_values(["panel_label", "scenario"]).to_csv(
            outdir / "benchmark_singleton_classification_summary.tsv",
            sep="\t",
            index=False,
        )
    else:
        pd.DataFrame().to_csv(outdir / "benchmark_singleton_classification_summary.tsv", sep="\t", index=False)
    return dataset_overview, genome_overview, benchmark_results


def replicate_existing_taxonomic_benchmark_panel(
    benchmark_dir: Path,
    outdir: Path,
    *,
    seed: int,
    num_cpus: int,
) -> None:
    base_manifest = json.loads((benchmark_dir / "benchmark_manifest.json").read_text())
    truth_markers = [str(marker) for marker in base_manifest.get("selected_markers", [])]
    if not truth_markers:
        raise ValueError(f"Benchmark manifest is missing selected markers: {benchmark_dir}")

    truth_inputs_dir, truth_run_dir, scenarios_dir = _prepare_replicate_output_dir(benchmark_dir, outdir)
    truth_models = outdir / "truth_markers.hmm"
    truth_records = _read_normalized_proteomes(truth_inputs_dir)
    truth_table = _load_table(truth_run_dir / "table_elim_dups")
    truth_tree = _load_truth_tree(truth_run_dir / "tree.nwk")
    selected_genomes = sorted(leaf.name for leaf in truth_tree.iter_leaves())
    truth_records = {genome: truth_records[genome] for genome in selected_genomes}
    native_map = _native_marker_map(truth_table, selected_genomes, truth_markers)
    native_map = {genome: native_map.get(genome, {}) for genome in selected_genomes}

    recipient_taxonomy = _taxonomy_map_from_manifest_rows(
        list(base_manifest.get("selected_genome_taxonomy", []))
    )
    if set(selected_genomes) - set(recipient_taxonomy):
        raise ValueError(f"Missing selected-genome taxonomy rows in manifest: {benchmark_dir}")

    truth_source_dir = _resolve_project_path(base_manifest.get("truth_source_dir"))
    donor_source_dir = _resolve_project_path(base_manifest.get("donor_source_dir")) or truth_source_dir
    taxonomy_db_path = _resolve_project_path(base_manifest.get("taxonomy_db_path")) or DEFAULT_TAXONOMY_DB
    donor_stage_dir = (
        benchmark_dir / "stage_truth_full_model_clean"
        if truth_source_dir is not None and donor_source_dir is not None and donor_source_dir.resolve() == truth_source_dir.resolve()
        else benchmark_dir / "stage_donor_full_model_clean"
    )
    donor_records, donor_native_map = _load_staged_records_and_native_map(donor_stage_dir, truth_markers)
    donor_taxonomy = _taxonomy_lookup(_load_source_taxonomy(donor_source_dir, taxonomy_db_path))

    rng = Random(seed)
    manifest = deepcopy(base_manifest)
    manifest["seed"] = seed
    manifest["replicate_seed"] = seed
    manifest["replicate_source_benchmark_dir"] = str(benchmark_dir)
    manifest["replicate_parent_panel_id"] = _panel_id_from_manifest(base_manifest, benchmark_dir)
    manifest["replicate_parent_panel_label"] = _panel_label_from_manifest(base_manifest, benchmark_dir)
    manifest["scenarios"] = []

    taxonomic_scope = str(base_manifest["taxonomic_scope"])
    lineage_label = str(base_manifest.get("lineage_label", ""))
    donor_lineage_label = str(base_manifest.get("donor_lineage_label", lineage_label))

    for scenario_name, spec in DEFAULT_SCENARIOS.items():
        scenario_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
        reference_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
        used_pairs: set[tuple[str, str]] = set()
        events: list[dict] = []
        event_index = 1

        duplicate_target = _duplicate_event_target(spec)
        duplicate_candidates = [
            (recipient, marker)
            for recipient in selected_genomes
            for marker in truth_markers
            if marker in native_map.get(recipient, {})
            and _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
        ]
        rng.shuffle(duplicate_candidates)
        for recipient, marker in duplicate_candidates:
            if len(events) >= duplicate_target:
                break
            if (recipient, marker) in used_pairs:
                continue
            donor_candidates = _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
            if not donor_candidates:
                continue
            donor = rng.choice(donor_candidates)
            donor_record = donor_records[donor][donor_native_map[donor][marker]]
            _degrade_record_in_place(
                scenario_records[recipient],
                native_record_id=native_map[recipient][marker],
                fraction=spec["native_degrade_fraction"],
                rng=rng,
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
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": TAXONOMY_SCOPE_RULES[taxonomic_scope]["scope_label"],
                    "recipient_genome": recipient,
                    "recipient_group": lineage_label,
                    "marker": marker,
                    "native_record_id": native_map[recipient][marker],
                    "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                    "donor_genome": donor,
                    "donor_group": donor_lineage_label,
                    "source_relation": _source_relation(lineage_label, donor_lineage_label),
                    "donor_record_id": donor_native_map[donor][marker],
                    "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                    "contaminant_record_id": contaminant.id,
                    "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                    "expected_duplicate_status": "Removed",
                    "expected_native_status": "Kept",
                    "native_degrade_fraction": spec["native_degrade_fraction"],
                    **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                    **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
                }
            )
            event_index += 1
        if len([row for row in events if row["event_type"] == "duplicate"]) != duplicate_target:
            raise ValueError(
                f"Could not schedule {duplicate_target} duplicate events for scope '{taxonomic_scope}'"
            )

        replacement_pairs = [
            (recipient, marker)
            for recipient in selected_genomes
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
            and marker in native_map.get(recipient, {})
            and _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
        ]
        rng.shuffle(replacement_pairs)
        replacement_events_written = 0
        for recipient, marker in replacement_pairs:
            if replacement_events_written >= spec["replacement_events"]:
                break
            donor_candidates = _taxonomic_donor_candidates(
                recipient,
                marker,
                taxonomic_scope,
                recipient_taxonomy,
                donor_taxonomy,
                donor_native_map,
                truth_tree,
            )
            if not donor_candidates:
                continue
            donor = rng.choice(donor_candidates)
            donor_record = donor_records[donor][donor_native_map[donor][marker]]
            contaminant = make_contaminant_record(
                recipient_genome=recipient,
                donor_record=donor_record,
                marker=marker,
                donor_genome=donor,
                event_index=event_index,
            )
            reference_records[recipient] = drop_native_marker(
                reference_records[recipient],
                native_record_id=native_map[recipient][marker],
            )
            scenario_records[recipient] = apply_replacement_event(
                scenario_records[recipient],
                native_record_id=native_map[recipient][marker],
                contaminant_record=contaminant,
            )
            used_pairs.add((recipient, marker))
            events.append(
                {
                    "event_index": event_index,
                    "scenario": scenario_name,
                    "event_type": "replacement",
                    "taxonomic_scope": taxonomic_scope,
                    "taxonomic_scope_label": TAXONOMY_SCOPE_RULES[taxonomic_scope]["scope_label"],
                    "recipient_genome": recipient,
                    "recipient_group": lineage_label,
                    "marker": marker,
                    "native_record_id": native_map[recipient][marker],
                    "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                    "donor_genome": donor,
                    "donor_group": donor_lineage_label,
                    "source_relation": _source_relation(lineage_label, donor_lineage_label),
                    "donor_record_id": donor_native_map[donor][marker],
                    "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                    "contaminant_record_id": contaminant.id,
                    "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                    "expected_replacement_outcome": "DropMarkerOrRemoveContaminant",
                    "native_degrade_fraction": spec["native_degrade_fraction"],
                    **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                    **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
                }
            )
            replacement_events_written += 1
            event_index += 1
        if replacement_events_written != spec["replacement_events"]:
            raise ValueError(
                f"Could not schedule {spec['replacement_events']} replacement events for scope '{taxonomic_scope}'"
            )

        scenario_dir = scenarios_dir / scenario_name
        scenario_dir.mkdir(exist_ok=True)
        proteomes_dir = scenario_dir / "proteomes"
        _write_proteome_dir(scenario_records, proteomes_dir)
        _reference_inputs_dir, reference_run_dir = _write_reference_state(
            scenario_dir=scenario_dir,
            reference_records=reference_records,
            truth_inputs_dir=truth_inputs_dir,
            truth_run_dir=truth_run_dir,
            truth_models=truth_models,
            num_cpus=num_cpus,
            reuse_truth_panel=spec["replacement_events"] == 0,
        )
        _write_events_tsv(scenario_dir / "events.tsv", events)
        genome_summary = _write_genome_summary_tsv(scenario_dir / "genome_summary.tsv", events)
        summary = _scenario_summary(events, genome_summary)
        manifest["scenarios"].append(
            {
                "name": scenario_name,
                "taxonomic_scope": taxonomic_scope,
                "proteomes_dir": str(proteomes_dir),
                "reference_inputs_dir": str(scenario_dir / "reference_inputs"),
                "reference_tree_path": str(reference_run_dir / "tree.nwk"),
                "reference_taxa": selected_genomes,
                "events_path": str(scenario_dir / "events.tsv"),
                "genome_summary_path": str(scenario_dir / "genome_summary.tsv"),
                **summary,
            }
        )

    _write_manifest_json(outdir / "benchmark_manifest.json", manifest)


def replicate_existing_mixed_benchmark_panel(
    benchmark_dir: Path,
    outdir: Path,
    *,
    seed: int,
    num_cpus: int,
) -> None:
    base_manifest = json.loads((benchmark_dir / "benchmark_manifest.json").read_text())
    truth_markers = [str(marker) for marker in base_manifest.get("selected_markers", [])]
    if not truth_markers:
        raise ValueError(f"Benchmark manifest is missing selected markers: {benchmark_dir}")

    truth_inputs_dir, truth_run_dir, scenarios_dir = _prepare_replicate_output_dir(benchmark_dir, outdir)
    truth_models = outdir / "truth_markers.hmm"
    truth_records = _read_normalized_proteomes(truth_inputs_dir)
    truth_table = _load_table(truth_run_dir / "table_elim_dups")
    truth_tree = _load_truth_tree(truth_run_dir / "tree.nwk")
    selected_genomes = sorted(leaf.name for leaf in truth_tree.iter_leaves())
    truth_records = {genome: truth_records[genome] for genome in selected_genomes}
    native_map = _native_marker_map(truth_table, selected_genomes, truth_markers)
    native_map = {genome: native_map.get(genome, {}) for genome in selected_genomes}

    recipient_taxonomy = _taxonomy_map_from_manifest_rows(
        list(base_manifest.get("selected_genome_taxonomy", []))
    )
    if set(selected_genomes) - set(recipient_taxonomy):
        raise ValueError(f"Missing selected-genome taxonomy rows in manifest: {benchmark_dir}")

    taxonomy_db_path = _resolve_project_path(base_manifest.get("taxonomy_db_path")) or DEFAULT_TAXONOMY_DB
    donor_source_dirs = {
        str(label): _resolve_project_path(path)
        for label, path in dict(base_manifest.get("donor_source_dirs", {})).items()
    }
    donor_bundles: dict[str, dict[str, object]] = {}
    for donor_label, donor_source_dir in donor_source_dirs.items():
        stage_dir = benchmark_dir / f"stage_donor_{sanitize_token(donor_label, donor_label)}"
        donor_records, donor_native_map = _load_staged_records_and_native_map(stage_dir, truth_markers)
        donor_bundles[donor_label] = {
            "records": donor_records,
            "native_map": donor_native_map,
            "taxonomy": _taxonomy_lookup(_load_source_taxonomy(donor_source_dir, taxonomy_db_path)),
            "source_dir": str(donor_source_dir),
        }

    rng = Random(seed)
    manifest = deepcopy(base_manifest)
    manifest["seed"] = seed
    manifest["replicate_seed"] = seed
    manifest["replicate_source_benchmark_dir"] = str(benchmark_dir)
    manifest["replicate_parent_panel_id"] = _panel_id_from_manifest(base_manifest, benchmark_dir)
    manifest["replicate_parent_panel_label"] = _panel_label_from_manifest(base_manifest, benchmark_dir)
    manifest["scenarios"] = []

    duplicate_recipients_list, replacement_recipients_list = _choose_recipient_sets(
        selected_genomes,
        duplicate_recipients=int(base_manifest.get("duplicate_recipients", 5)),
        replacement_recipients=int(base_manifest.get("replacement_recipients", 5)),
        overlap_recipients=int(base_manifest.get("overlap_recipients", 2)),
        rng=rng,
    )
    donor_labels = sorted(donor_bundles)
    duplicate_label_plan = _round_robin_labels(
        donor_labels,
        int(base_manifest.get("duplicate_recipients", 5)),
    )
    replacement_label_plan = _round_robin_labels(
        list(reversed(donor_labels)),
        int(base_manifest.get("replacement_recipients", 5)),
    )

    scenario_name = "mixed_high_level"
    scenario_dir = scenarios_dir / scenario_name
    scenario_dir.mkdir(exist_ok=True)
    scenario_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
    reference_records = {genome: deepcopy(records) for genome, records in truth_records.items()}
    used_pairs: set[tuple[str, str]] = set()
    events: list[dict] = []
    event_index = 1
    lineage_label = str(base_manifest.get("lineage_label", ""))

    for recipient, donor_label in zip(duplicate_recipients_list, duplicate_label_plan):
        donor_bundle = donor_bundles[donor_label]
        donor_records = donor_bundle["records"]
        donor_native_map = donor_bundle["native_map"]
        donor_taxonomy = donor_bundle["taxonomy"]
        marker_candidates = [
            marker
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
            and marker in native_map.get(recipient, {})
            and any(marker in donor_native_map.get(donor_genome, {}) for donor_genome in donor_records)
        ]
        if not marker_candidates:
            raise ValueError(f"No duplicate-event marker candidates for recipient {recipient} in mixed benchmark")
        rng.shuffle(marker_candidates)
        marker = marker_candidates[0]
        donor_candidates = [
            donor_genome
            for donor_genome in donor_records
            if marker in donor_native_map.get(donor_genome, {})
        ]
        donor = rng.choice(sorted(donor_candidates))
        donor_record = donor_records[donor][donor_native_map[donor][marker]]
        _degrade_record_in_place(
            scenario_records[recipient],
            native_record_id=native_map[recipient][marker],
            fraction=0.12,
            rng=rng,
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
                "taxonomic_scope": "mixed_high_level",
                "taxonomic_scope_label": f"{'+'.join(donor_labels)}_into_{lineage_label}",
                "recipient_genome": recipient,
                "recipient_group": lineage_label,
                "marker": marker,
                "native_record_id": native_map[recipient][marker],
                "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                "donor_genome": donor,
                "donor_group": donor_label,
                "source_relation": _source_relation(lineage_label, donor_label),
                "donor_record_id": donor_native_map[donor][marker],
                "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                "contaminant_record_id": contaminant.id,
                "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                "expected_duplicate_status": "Removed",
                "expected_native_status": "Kept",
                "native_degrade_fraction": 0.12,
                **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
            }
        )
        event_index += 1

    for recipient, donor_label in zip(replacement_recipients_list, replacement_label_plan):
        donor_bundle = donor_bundles[donor_label]
        donor_records = donor_bundle["records"]
        donor_native_map = donor_bundle["native_map"]
        donor_taxonomy = donor_bundle["taxonomy"]
        marker_candidates = [
            marker
            for marker in truth_markers
            if (recipient, marker) not in used_pairs
            and marker in native_map.get(recipient, {})
            and any(marker in donor_native_map.get(donor_genome, {}) for donor_genome in donor_records)
        ]
        if not marker_candidates:
            raise ValueError(f"No replacement-event marker candidates for recipient {recipient} in mixed benchmark")
        rng.shuffle(marker_candidates)
        marker = marker_candidates[0]
        donor_candidates = [
            donor_genome
            for donor_genome in donor_records
            if marker in donor_native_map.get(donor_genome, {})
        ]
        donor = rng.choice(sorted(donor_candidates))
        donor_record = donor_records[donor][donor_native_map[donor][marker]]
        contaminant = make_contaminant_record(
            recipient_genome=recipient,
            donor_record=donor_record,
            marker=marker,
            donor_genome=donor,
            event_index=event_index,
        )
        reference_records[recipient] = drop_native_marker(
            reference_records[recipient],
            native_record_id=native_map[recipient][marker],
        )
        scenario_records[recipient] = apply_replacement_event(
            scenario_records[recipient],
            native_record_id=native_map[recipient][marker],
            contaminant_record=contaminant,
        )
        used_pairs.add((recipient, marker))
        events.append(
            {
                "event_index": event_index,
                "scenario": scenario_name,
                "event_type": "replacement",
                "taxonomic_scope": "mixed_high_level",
                "taxonomic_scope_label": f"{'+'.join(donor_labels)}_into_{lineage_label}",
                "recipient_genome": recipient,
                "recipient_group": lineage_label,
                "marker": marker,
                "native_record_id": native_map[recipient][marker],
                "native_contig_id": parse_sequence_id(native_map[recipient][marker])[1],
                "donor_genome": donor,
                "donor_group": donor_label,
                "source_relation": _source_relation(lineage_label, donor_label),
                "donor_record_id": donor_native_map[donor][marker],
                "donor_contig_id": parse_sequence_id(donor_native_map[donor][marker])[1],
                "contaminant_record_id": contaminant.id,
                "contaminant_contig_id": parse_sequence_id(contaminant.id)[1],
                "expected_replacement_outcome": "DropMarkerOrRemoveContaminant",
                "native_degrade_fraction": 0.12,
                **_event_taxonomy_fields("recipient", recipient_taxonomy[recipient]),
                **_event_taxonomy_fields("donor", donor_taxonomy[donor]),
            }
        )
        event_index += 1

    proteomes_dir = scenario_dir / "proteomes"
    _write_proteome_dir(scenario_records, proteomes_dir)
    _reference_inputs_dir, reference_run_dir = _write_reference_state(
        scenario_dir=scenario_dir,
        reference_records=reference_records,
        truth_inputs_dir=truth_inputs_dir,
        truth_run_dir=truth_run_dir,
        truth_models=truth_models,
        num_cpus=num_cpus,
        reuse_truth_panel=False,
    )
    _write_events_tsv(scenario_dir / "events.tsv", events)
    genome_summary = _write_genome_summary_tsv(scenario_dir / "genome_summary.tsv", events)
    summary = _scenario_summary(events, genome_summary)
    manifest["scenarios"].append(
        {
            "name": scenario_name,
            "taxonomic_scope": "mixed_high_level",
            "proteomes_dir": str(proteomes_dir),
            "reference_inputs_dir": str(scenario_dir / "reference_inputs"),
            "reference_tree_path": str(reference_run_dir / "tree.nwk"),
            "reference_taxa": selected_genomes,
            "events_path": str(scenario_dir / "events.tsv"),
            "genome_summary_path": str(scenario_dir / "genome_summary.tsv"),
            **summary,
        }
    )
    _write_manifest_json(outdir / "benchmark_manifest.json", manifest)


def replicate_existing_benchmark_panel(
    benchmark_dir: Path,
    outdir: Path,
    *,
    seed: int,
    num_cpus: int,
) -> None:
    manifest = json.loads((benchmark_dir / "benchmark_manifest.json").read_text())
    scope = str(manifest.get("taxonomic_scope", ""))
    if scope == "mixed_high_level":
        replicate_existing_mixed_benchmark_panel(
            benchmark_dir=benchmark_dir,
            outdir=outdir,
            seed=seed,
            num_cpus=num_cpus,
        )
        return
    if scope in TAXONOMY_SCOPE_RULES:
        replicate_existing_taxonomic_benchmark_panel(
            benchmark_dir=benchmark_dir,
            outdir=outdir,
            seed=seed,
            num_cpus=num_cpus,
        )
        return
    raise ValueError(f"Unsupported benchmark replicate source: {benchmark_dir}")


def generate_panel_replicates(
    benchmark_dirs: list[Path],
    outbase: Path,
    replicate_seeds: list[int],
    *,
    num_cpus: int,
    max_parallel: int = 1,
) -> list[Path]:
    outbase.mkdir(parents=True, exist_ok=True)
    replicate_dirs: list[Path] = []
    replicate_rows: list[dict] = []
    jobs: list[tuple[Path, callable]] = []

    for benchmark_dir in benchmark_dirs:
        manifest = json.loads((benchmark_dir / "benchmark_manifest.json").read_text())
        panel_id = _panel_id_from_manifest(manifest, benchmark_dir)
        panel_label = _panel_label_from_manifest(manifest, benchmark_dir)
        for seed in replicate_seeds:
            replicate_dir = outbase / f"{panel_id}__seed_{seed}"
            replicate_dirs.append(replicate_dir)
            replicate_rows.append(
                {
                    "panel_id": panel_id,
                    "panel_label": panel_label,
                    "source_benchmark_dir": str(benchmark_dir),
                    "replicate_seed": seed,
                    "replicate_dir": str(replicate_dir),
                }
            )
            jobs.append(
                (
                    benchmark_dir,
                    lambda benchmark_dir=benchmark_dir, replicate_dir=replicate_dir, seed=seed: replicate_existing_benchmark_panel(
                        benchmark_dir=benchmark_dir,
                        outdir=replicate_dir,
                        seed=seed,
                        num_cpus=num_cpus,
                    ),
                )
            )

    _parallel_benchmark_jobs(jobs, max_parallel=max_parallel)
    _write_manifest_json(
        outbase / "replicates_manifest.json",
        {
            "replicate_seeds": replicate_seeds,
            "num_cpus": num_cpus,
            "max_parallel_generation": max_parallel,
            "source_benchmark_dirs": [str(path) for path in benchmark_dirs],
            "replicates": replicate_rows,
        },
    )
    return replicate_dirs


def _benchmark_dirs_from_replicates_manifest(manifest_path: Path) -> list[Path]:
    payload = json.loads(manifest_path.read_text())
    return [Path(str(item["replicate_dir"])) for item in payload.get("replicates", []) if item.get("replicate_dir")]


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

    gentax = subparsers.add_parser("generate-taxonomic", help="Generate a taxonomy-aware benchmark dataset")
    gentax.add_argument("--truth-source-dir", required=True)
    gentax.add_argument("--donor-source-dir", default=None)
    gentax.add_argument("--models", default="resources/models/UNI56.hmm")
    gentax.add_argument("--outdir", required=True)
    gentax.add_argument("--taxonomy-scope", required=True, choices=sorted(TAXONOMY_SCOPE_RULES))
    gentax.add_argument("--taxonomy-db", default=str(DEFAULT_TAXONOMY_DB))
    gentax.add_argument("--lineage-label", default=None)
    gentax.add_argument("--donor-lineage-label", default=None)
    gentax.add_argument("--n-genomes", type=int, default=50)
    gentax.add_argument("--n-markers", type=int, default=10)
    gentax.add_argument("--min-marker-presence-fraction", type=float, default=0.8)
    gentax.add_argument("--seed", type=int, default=42)
    gentax.add_argument("--num-cpus", type=int, default=8)

    gensuite = subparsers.add_parser("generate-taxonomic-suite", help="Generate genus/family/order taxonomy-aware benchmark datasets")
    gensuite.add_argument("--truth-source-dir", required=True)
    gensuite.add_argument("--outbase", required=True)
    gensuite.add_argument("--models", default="resources/models/UNI56.hmm")
    gensuite.add_argument("--taxonomy-db", default=str(DEFAULT_TAXONOMY_DB))
    gensuite.add_argument("--lineage-label", default=None)
    gensuite.add_argument("--scopes", nargs="+", default=sorted(TAXONOMY_SCOPE_RULES))
    gensuite.add_argument("--genus-donor-source-dir", default=None)
    gensuite.add_argument("--family-donor-source-dir", default=None)
    gensuite.add_argument("--order-donor-source-dir", default=None)
    gensuite.add_argument("--genus-donor-lineage-label", default=None)
    gensuite.add_argument("--family-donor-lineage-label", default=None)
    gensuite.add_argument("--order-donor-lineage-label", default=None)
    gensuite.add_argument("--n-genomes", type=int, default=50)
    gensuite.add_argument("--n-markers", type=int, default=10)
    gensuite.add_argument("--min-marker-presence-fraction", type=float, default=0.8)
    gensuite.add_argument("--seed", type=int, default=42)
    gensuite.add_argument("--num-cpus", type=int, default=8)

    genmixed = subparsers.add_parser("generate-mixed", help="Generate a mixed high-level cross-lineage benchmark dataset")
    genmixed.add_argument("--truth-source-dir", required=True)
    genmixed.add_argument("--donor-a-source-dir", required=True)
    genmixed.add_argument("--donor-b-source-dir", required=True)
    genmixed.add_argument("--donor-a-label", required=True)
    genmixed.add_argument("--donor-b-label", required=True)
    genmixed.add_argument("--models", default="resources/models/UNI56.hmm")
    genmixed.add_argument("--outdir", required=True)
    genmixed.add_argument("--taxonomy-db", default=str(DEFAULT_TAXONOMY_DB))
    genmixed.add_argument("--lineage-label", default=None)
    genmixed.add_argument("--n-genomes", type=int, default=50)
    genmixed.add_argument("--n-markers", type=int, default=10)
    genmixed.add_argument("--min-marker-presence-fraction", type=float, default=0.8)
    genmixed.add_argument("--duplicate-recipients", type=int, default=5)
    genmixed.add_argument("--replacement-recipients", type=int, default=5)
    genmixed.add_argument("--overlap-recipients", type=int, default=2)
    genmixed.add_argument("--seed", type=int, default=42)
    genmixed.add_argument("--num-cpus", type=int, default=8)

    prep = subparsers.add_parser("prepare-source-subset", help="Prepare a sampled source subset directory")
    prep.add_argument("--source-dir", required=True)
    prep.add_argument("--outdir", required=True)
    prep.add_argument("--list-file", default=None)
    prep.add_argument("--n-candidates", type=int, default=None)
    prep.add_argument("--seed", type=int, default=42)

    prepburk = subparsers.add_parser("prepare-burkholderiaceae", help="Materialize the 50-genome Burkholderiaceae ANI benchmark panel")
    prepburk.add_argument("--outdir", default="testgenomes/Burkholderiaceae50")
    prepburk.add_argument("--lookup", default="testgenomes/burkholderiaceae50.lookup")
    prepburk.add_argument("--taxonomy-tsv", default="testgenomes/burkholderiaceae50_taxonomy.tsv")
    prepburk.add_argument("--selection-tsv", default="testgenomes/burkholderiaceae50_selection.tsv")
    prepburk.add_argument("--taxonomy-db", default=str(DEFAULT_TAXONOMY_DB))
    prepburk.add_argument("--prefix", default="BURK__")
    prepburk.add_argument("--overwrite", action="store_true")

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

    genpanelrep = subparsers.add_parser("generate-panel-replicates", help="Generate replicate copies from existing benchmark panels")
    genpanelrep.add_argument("--benchmark-dirs", nargs="+", default=None)
    genpanelrep.add_argument("--benchmarks-root", default="runs/benchmarks")
    genpanelrep.add_argument("--outbase", required=True)
    genpanelrep.add_argument("--replicate-seeds", nargs="+", type=int, required=True)
    genpanelrep.add_argument("--num-cpus", type=int, default=8)
    genpanelrep.add_argument("--max-parallel", type=int, default=1)

    runrep = subparsers.add_parser("run-replicates", help="Run benchmark suite for multiple replicate dirs")
    runrep.add_argument("--benchmark-dirs", nargs="+", default=None)
    runrep.add_argument("--replicates-manifest", default=None)
    runrep.add_argument("--num-cpus", type=int, default=8)
    runrep.add_argument("--selection-max-rounds", type=int, default=5)
    runrep.add_argument("--max-parallel", type=int, default=1)

    cross = subparsers.add_parser("generate-cross", help="Generate a cross-family benchmark from two benchmark dirs")
    cross.add_argument("--benchmark-a-dir", required=True)
    cross.add_argument("--benchmark-b-dir", required=True)
    cross.add_argument("--outdir", required=True)
    cross.add_argument("--n-markers", type=int, default=8)
    cross.add_argument("--seed", type=int, default=42)
    cross.add_argument("--num-cpus", type=int, default=8)

    agg = subparsers.add_parser("aggregate-replicates", help="Aggregate multiple benchmark result directories")
    agg.add_argument("--benchmark-dirs", nargs="+", default=None)
    agg.add_argument("--replicates-manifest", default=None)
    agg.add_argument("--outdir", required=True)

    exportdocs = subparsers.add_parser("export-docs", help="Export benchmark overview tables for docs/data")
    exportdocs.add_argument("--benchmarks-root", default="runs/benchmarks")
    exportdocs.add_argument("--outdir", default="docs/data")

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
    def _discover_candidates() -> list[Path]:
        exts = ("*.faa", "*.fna", "*.fa", "*.fasta")
        seen = []
        for pattern in exts:
            seen.extend(sorted(source_dir.glob(pattern)))
        return sorted({path for path in seen if path.is_file()}, key=lambda path: path.name)

    if list_file is not None:
        ids = [line.strip() for line in list_file.read_text().splitlines() if line.strip()]
        for genome_id in ids:
            prefixed = source_dir / f"FLAV__{genome_id}.faa"
            if prefixed.exists():
                candidates.append(prefixed)
                continue
            for ext in (".fna", ".fa", ".fasta", ".faa"):
                candidate = source_dir / f"{genome_id}{ext}"
                if candidate.exists():
                    candidates.append(candidate)
                    break
    else:
        candidates = _discover_candidates()
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
    elif args.command == "generate-taxonomic":
        generate_taxonomic_benchmark_dataset(
            truth_source_dir=Path(args.truth_source_dir),
            donor_source_dir=Path(args.donor_source_dir) if args.donor_source_dir else None,
            models_path=Path(args.models),
            outdir=Path(args.outdir),
            taxonomic_scope=args.taxonomy_scope,
            taxonomy_db_path=Path(args.taxonomy_db),
            lineage_label=args.lineage_label,
            donor_lineage_label=args.donor_lineage_label,
            n_genomes=args.n_genomes,
            n_markers=args.n_markers,
            min_marker_presence_fraction=args.min_marker_presence_fraction,
            seed=args.seed,
            num_cpus=args.num_cpus,
        )
    elif args.command == "generate-taxonomic-suite":
        scopes = list(dict.fromkeys(args.scopes))
        generate_taxonomic_benchmark_suite(
            truth_source_dir=Path(args.truth_source_dir),
            outbase=Path(args.outbase),
            models_path=Path(args.models),
            taxonomy_db_path=Path(args.taxonomy_db),
            lineage_label=args.lineage_label,
            scopes=scopes,
            donor_source_by_scope={
                "genus": Path(args.genus_donor_source_dir) if args.genus_donor_source_dir else None,
                "family": Path(args.family_donor_source_dir) if args.family_donor_source_dir else None,
                "order": Path(args.order_donor_source_dir) if args.order_donor_source_dir else None,
            },
            donor_lineage_label_by_scope={
                "genus": args.genus_donor_lineage_label,
                "family": args.family_donor_lineage_label,
                "order": args.order_donor_lineage_label,
            },
            n_genomes=args.n_genomes,
            n_markers=args.n_markers,
            min_marker_presence_fraction=args.min_marker_presence_fraction,
            seed=args.seed,
            num_cpus=args.num_cpus,
        )
    elif args.command == "generate-mixed":
        generate_mixed_lineage_benchmark_dataset(
            truth_source_dir=Path(args.truth_source_dir),
            donor_source_dirs={
                args.donor_a_label: Path(args.donor_a_source_dir),
                args.donor_b_label: Path(args.donor_b_source_dir),
            },
            models_path=Path(args.models),
            outdir=Path(args.outdir),
            taxonomy_db_path=Path(args.taxonomy_db),
            lineage_label=args.lineage_label,
            n_genomes=args.n_genomes,
            n_markers=args.n_markers,
            min_marker_presence_fraction=args.min_marker_presence_fraction,
            duplicate_recipients=args.duplicate_recipients,
            replacement_recipients=args.replacement_recipients,
            overlap_recipients=args.overlap_recipients,
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
    elif args.command == "prepare-burkholderiaceae":
        selection = prepare_burkholderiaceae_benchmark_dataset(
            taxonomy_db_path=Path(args.taxonomy_db),
            output_dir=Path(args.outdir),
            lookup_path=Path(args.lookup),
            taxonomy_tsv_path=Path(args.taxonomy_tsv),
            selection_tsv_path=Path(args.selection_tsv),
            prefix=args.prefix,
            overwrite=args.overwrite,
        )
        print(
            f"Prepared Burkholderiaceae benchmark panel: {len(selection)} genomes, "
            f"{selection['genus'].nunique()} genera, {selection['species'].nunique()} species"
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
    elif args.command == "generate-panel-replicates":
        benchmark_dirs = (
            [Path(item) for item in args.benchmark_dirs]
            if args.benchmark_dirs
            else discover_taxonomic_benchmark_dirs(Path(args.benchmarks_root))
        )
        generate_panel_replicates(
            benchmark_dirs=benchmark_dirs,
            outbase=Path(args.outbase),
            replicate_seeds=args.replicate_seeds,
            num_cpus=args.num_cpus,
            max_parallel=args.max_parallel,
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
        benchmark_dirs = (
            _benchmark_dirs_from_replicates_manifest(Path(args.replicates_manifest))
            if args.replicates_manifest
            else [Path(item) for item in args.benchmark_dirs or []]
        )
        run_benchmark_suites(
            benchmark_dirs=benchmark_dirs,
            num_cpus=args.num_cpus,
            selection_max_rounds=args.selection_max_rounds,
            max_parallel=args.max_parallel,
        )
    elif args.command == "aggregate-replicates":
        benchmark_dirs = (
            _benchmark_dirs_from_replicates_manifest(Path(args.replicates_manifest))
            if args.replicates_manifest
            else [Path(item) for item in args.benchmark_dirs or []]
        )
        aggregate_benchmark_runs(
            benchmark_dirs=benchmark_dirs,
            outdir=Path(args.outdir),
        )
    elif args.command == "export-docs":
        benchmark_dirs = discover_taxonomic_benchmark_dirs(Path(args.benchmarks_root))
        export_benchmark_tables(
            benchmark_dirs=benchmark_dirs,
            outdir=Path(args.outdir),
        )


if __name__ == "__main__":
    main()
