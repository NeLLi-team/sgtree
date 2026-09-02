"""Synthetic benchmark generation and evaluation for SGTree."""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
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
        "name": "singles_gcp",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "gcp",
    },
    "combined": {
        "name": "duplicate_plus_singles_gcp",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "gcp",
    },
    "mixed_high_level": {
        "name": "duplicate_plus_singles_gcp",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
        "marker_selection": True,
        "singles": True,
        "singles_mode": "gcp",
    },
}

# Path to the local GTDB genome DuckDB used by the taxonomy-aware generators.
# There is no portable default, so it comes from the environment and the CLI
# requires --taxonomy-db when the variable is unset.
DEFAULT_TAXONOMY_DB = os.environ.get("SGTREE_TAXONOMY_DB")
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
    "class": {
        "same_rank": "phylum",
        "different_rank": "class",
        "scope_label": "same_phylum_different_class",
    },
    "phylum": {
        "same_rank": "domain",
        "different_rank": "phylum",
        "scope_label": "same_domain_different_phylum",
    },
}
TAXONOMY_DISTANCE_RANKS = (
    ("phylum", "phylum"),
    ("class", "class"),
    ("order", "order_name"),
    ("family", "family"),
    ("genus", "genus"),
    ("species", "species"),
)
TAXONOMY_DISTANCE_LABELS = {
    "phylum": "different_phylum",
    "class": "different_class_same_phylum",
    "order": "different_order_same_class",
    "family": "different_family_same_order",
    "genus": "different_genus_same_family",
    "species": "different_species_same_genus",
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
    cmd.extend(["--singles-mode", singles_mode])
    subprocess.run(cmd, check=True, env=_child_env())


def _child_env() -> dict[str, str]:
    """Environment for a child ``python -m sgtree`` call.

    The package is run from the source tree rather than installed, so a child
    process needs the src directory on PYTHONPATH. Without it the child exits
    with "No module named sgtree" whenever the parent was started by a task
    that does not already export it.
    """
    src_root = str(pathlib.Path(__file__).resolve().parents[2])
    env = dict(os.environ)
    existing = env.get("PYTHONPATH", "")
    if src_root not in existing.split(os.pathsep):
        env["PYTHONPATH"] = os.pathsep.join(p for p in (src_root, existing) if p)
    return env


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

    # Avoid pandas' nullable/string extension dtype here; DuckDB register()
    # expects a plain object-backed string column.
    request_df = pd.DataFrame(
        {"assembly_accession": pd.Series(sorted(set(assembly_accessions)), dtype="object")}
    )
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
    if key == "order_name" and key not in meta:
        key = "order"
    value = str(meta.get(key, "")).strip()
    return "" if value.lower() in {"nan", "none", "<na>"} else value


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


def classify_taxonomic_distance(
    recipient_meta: dict[str, str],
    donor_meta: dict[str, str],
) -> dict[str, str]:
    """Classify the donor-vs-recipient source rank for a contamination event."""
    differs: dict[str, str] = {}
    first_known_difference = ""
    unknown_seen = False

    for rank_label, rank_key in TAXONOMY_DISTANCE_RANKS:
        recipient_value = _taxonomy_field(recipient_meta, rank_key)
        donor_value = _taxonomy_field(donor_meta, rank_key)
        if not recipient_value or not donor_value:
            differs[rank_label] = "unknown"
            unknown_seen = True
        elif recipient_value != donor_value:
            differs[rank_label] = "yes"
            if not first_known_difference:
                first_known_difference = rank_label
        else:
            differs[rank_label] = "no"

    if first_known_difference:
        level = first_known_difference
        label = TAXONOMY_DISTANCE_LABELS[level]
    elif unknown_seen:
        level = "unknown"
        label = "unknown"
    else:
        level = "same_taxonomy"
        label = "same_taxonomy"

    return {
        "contamination_source_taxonomic_level": level,
        "contamination_source_taxonomic_label": label,
        **{
            f"contamination_source_differs_{rank_label}": differs[rank_label]
            for rank_label, _rank_key in TAXONOMY_DISTANCE_RANKS
        },
    }


def _event_taxonomic_distance_fields(
    recipient_meta: dict[str, str],
    donor_meta: dict[str, str],
) -> dict[str, str]:
    return classify_taxonomic_distance(recipient_meta, donor_meta)


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

        duplicate_target = int(spec["pair_blocks"]) * int(spec["markers_per_block"]) * 2
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
                    **_event_taxonomic_distance_fields(recipient_taxonomy[recipient], donor_taxonomy[donor]),
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
                    **_event_taxonomic_distance_fields(recipient_taxonomy[recipient], donor_taxonomy[donor]),
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


def _replacement_outcome(
    run_dir: Path,
    recipient: str,
    marker: str,
    contaminant_record_id: str,
    native_record_id: str,
) -> str:
    aligned_path = run_dir / "aligned_final" / f"{marker}.faa"
    if not aligned_path.exists():
        return "unknown"
    with open(aligned_path) as handle:
        ids = [record.id for record in SeqIO.parse(handle, "fasta")]
    if not ids:
        return "unknown"
    if contaminant_record_id in ids:
        return "contaminant_retained"
    if native_record_id in ids:
        return "native_retained"
    if not any(record_id.startswith(f"{recipient}|") for record_id in ids):
        return "marker_dropped"
    return "unknown"


def _sorted_leaf_names(tree_path: Path) -> list[str]:
    tree = Tree(str(tree_path))
    return sorted(leaf.name for leaf in tree.iter_leaves())


def _format_taxa_list(names: list[str]) -> str:
    return ",".join(sorted(names))


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
        marker = path.stem.split("_")[-2]
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

    reference_taxa = _expected_reference_taxa(manifest, scenario_meta, truth_tree)
    final_taxa = _sorted_leaf_names(result_tree)
    final_taxa_set = set(final_taxa)
    replacement_outcomes = []
    for row in replacement_events.itertuples(index=False):
        outcome = _replacement_outcome(
            run_dir,
            str(row.recipient_genome),
            str(row.marker),
            str(row.contaminant_record_id),
            str(row.native_record_id),
        )
        if outcome in {"marker_dropped", "native_retained"} and str(row.recipient_genome) not in final_taxa_set:
            outcome = "recipient_lost"
        replacement_outcomes.append(outcome)
    replacement_contaminant_retained = int(sum(outcome == "contaminant_retained" for outcome in replacement_outcomes))
    replacement_contaminant_removed = int(
        sum(outcome in {"marker_dropped", "native_retained"} for outcome in replacement_outcomes)
    )
    total_contaminants = int(len(duplicate_events) + len(replacement_events))
    total_removed = int(contaminant_correct + replacement_contaminant_removed)
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
        "replacement_unknown": int(sum(outcome == "unknown" for outcome in replacement_outcomes)),
        "replacement_recipient_lost": int(sum(outcome == "recipient_lost" for outcome in replacement_outcomes)),
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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate and run SGTree contamination benchmarks")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gen = subparsers.add_parser("generate", help="Generate a benchmark dataset")
    gen.add_argument("--source-dir", default="benchmarking/testgenomes/Chloroflexi")
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
    gentax.add_argument("--taxonomy-db", default=DEFAULT_TAXONOMY_DB,
                        required=DEFAULT_TAXONOMY_DB is None)
    gentax.add_argument("--lineage-label", default=None)
    gentax.add_argument("--donor-lineage-label", default=None)
    gentax.add_argument("--n-genomes", type=int, default=50)
    gentax.add_argument("--n-markers", type=int, default=10)
    gentax.add_argument("--min-marker-presence-fraction", type=float, default=0.8)
    gentax.add_argument("--seed", type=int, default=42)
    gentax.add_argument("--num-cpus", type=int, default=8)

    prepburk = subparsers.add_parser("prepare-burkholderiaceae", help="Materialize the 50-genome Burkholderiaceae ANI benchmark panel")
    prepburk.add_argument("--outdir", default="benchmarking/testgenomes/Burkholderiaceae50")
    prepburk.add_argument("--lookup", default="benchmarking/testgenomes/burkholderiaceae50.lookup")
    prepburk.add_argument("--taxonomy-tsv", default="benchmarking/testgenomes/burkholderiaceae50_taxonomy.tsv")
    prepburk.add_argument("--selection-tsv", default="benchmarking/testgenomes/burkholderiaceae50_selection.tsv")
    prepburk.add_argument("--taxonomy-db", default=DEFAULT_TAXONOMY_DB,
                          required=DEFAULT_TAXONOMY_DB is None)
    prepburk.add_argument("--prefix", default="BURK__")
    prepburk.add_argument("--overwrite", action="store_true")

    run = subparsers.add_parser("run", help="Run legacy vs coordinate benchmark comparisons")
    run.add_argument("--benchmark-dir", default="runs/benchmarks/dev_chloroflexi")
    run.add_argument("--num-cpus", type=int, default=8)
    run.add_argument("--selection-max-rounds", type=int, default=5)

    return parser.parse_args()


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
    elif args.command == "run":
        run_benchmark_suite(
            benchmark_dir=Path(args.benchmark_dir),
            num_cpus=args.num_cpus,
            selection_max_rounds=args.selection_max_rounds,
        )


if __name__ == "__main__":
    main()
