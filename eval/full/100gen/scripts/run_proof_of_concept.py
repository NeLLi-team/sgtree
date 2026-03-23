#!/usr/bin/env python3
"""Proof-of-concept benchmark: Alphaproteobacteria and Bacteroidota.

Generates benchmark panels across multiple taxonomic scopes and seed
replicates, runs SGTree in 3 contig variants (solo, merged_contam,
native_contig), and evaluates results.

Full matrix:
  2 lineages × 3 scopes × 3 seeds = 18 panels
  18 panels × 3 scenarios × 3 variants = 162 SGTree runs

Usage:
    pixi run python runs/contig_variant_proof_of_concept/run_proof_of_concept.py \
        --step prepare     # Step 1: download genomes
        --step generate    # Step 2: generate benchmark panels
        --step run         # Step 3: run SGTree on all variants
        --step evaluate    # Step 4: evaluate results
        --step all         # Run all steps
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import duckdb
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNS_DIR = Path(__file__).resolve().parent
TAXONOMY_DB = Path("/home/fschulz/dev/nelli-genomes-db/resources/database/gtdb_genomes.duckdb")
MODELS_PATH = PROJECT_ROOT / "resources" / "models" / "UNI56.hmm"
NUM_CPUS = 8

SCOPES = ["genus", "family", "order", "phylum"]

# For phylum scope, each lineage uses the other as donor source.
CROSS_LINEAGE_DONOR = {
    "alpha": "bacteroidota",
    "bacteroidota": "alpha",
}
SEEDS = [101, 202, 303]

DOUBLED_SCENARIOS = {
    "duplicate_only": {
        "pair_blocks": 4,
        "markers_per_block": 2,
        "replacement_events": 0,
        "native_degrade_fraction": 0.08,
    },
    "replacement_only": {
        "pair_blocks": 0,
        "markers_per_block": 0,
        "replacement_events": 8,
        "native_degrade_fraction": 0.12,
    },
    "combined": {
        "pair_blocks": 4,
        "markers_per_block": 2,
        "replacement_events": 8,
        "native_degrade_fraction": 0.12,
    },
}

CLEANUP_PROFILES = {
    "duplicate_only": {
        "name": "duplicate_cleanup",
        "marker_selection": True,
        "singles": False,
        "singles_mode": "delta_rf",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
    },
    "replacement_only": {
        "name": "singles_delta_rf",
        "marker_selection": True,
        "singles": True,
        "singles_mode": "delta_rf",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
    },
    "combined": {
        "name": "duplicate_plus_singles_delta_rf",
        "marker_selection": True,
        "singles": True,
        "singles_mode": "delta_rf",
        "selection_mode": "coordinate",
        "selection_global_rounds": 2,
    },
}

LINEAGES = {
    "alpha": {
        "prefix": "ALPHA__",
        "label": "alpha",
        "query_filter": "class = 'Alphaproteobacteria'",
        "n_genomes": 100,
        "description": "Alphaproteobacteria (100 genomes)",
    },
    "bacteroidota": {
        "prefix": "BACTE__",
        "label": "bacteroidota",
        "query_filter": "phylum = 'Bacteroidota'",
        "n_genomes": 100,
        "description": "Bacteroidota (100 genomes)",
    },
}


def panel_name(scope: str, seed: int) -> str:
    return f"{scope}__seed_{seed}"


def iter_panels():
    """Yield (lineage_key, cfg, scope, seed, panel_dir) for all panels."""
    for lineage_key, cfg in LINEAGES.items():
        for scope in SCOPES:
            for seed in SEEDS:
                pname = panel_name(scope, seed)
                panel_dir = RUNS_DIR / "benchmarks" / lineage_key / pname
                yield lineage_key, cfg, scope, seed, panel_dir


def select_diverse_genomes(
    db_path: Path,
    query_filter: str,
    n_genomes: int,
    prefix: str,
    seed: int = 42,
) -> pd.DataFrame:
    """Select diverse genomes from the taxonomy DB, maximizing genus coverage."""
    db = duckdb.connect(str(db_path), read_only=True)

    shared_cols = (
        "assembly_accession, organism_name, ftp_path, genome_size_bp, "
        "assembly_level, refseq_category, phylum, class, order_name, "
        "family, genus, species"
    )
    candidates = db.execute(f"""
        WITH ranked AS (
            SELECT *,
                ROW_NUMBER() OVER (
                    PARTITION BY assembly_accession
                    ORDER BY CASE WHEN taxonomy_source = 'gtdb' THEN 0 ELSE 1 END
                ) AS rn
            FROM (
                SELECT {shared_cols}, 'gtdb' AS taxonomy_source FROM gtdb_genomes
                UNION ALL
                SELECT {shared_cols}, 'non_gtdb' AS taxonomy_source FROM non_gtdb_genomes
            )
            WHERE {query_filter}
            AND ftp_path IS NOT NULL AND ftp_path != ''
        )
        SELECT assembly_accession, organism_name, ftp_path, genome_size_bp,
               assembly_level, refseq_category, phylum, class, order_name,
               family, genus, species, taxonomy_source
        FROM ranked WHERE rn = 1
    """).fetchdf()
    db.close()

    if len(candidates) < n_genomes:
        raise ValueError(
            f"Only {len(candidates)} candidates for '{query_filter}', need {n_genomes}"
        )

    level_rank = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    cat_rank = {"reference genome": 0, "representative genome": 1}
    candidates["level_rank"] = candidates["assembly_level"].map(level_rank).fillna(4)
    candidates["cat_rank"] = candidates["refseq_category"].map(cat_rank).fillna(2)

    genera = sorted(candidates["genus"].unique())
    selected_indices = []
    genera_pools = {}
    for g in genera:
        pool = candidates[candidates["genus"] == g].sort_values(
            ["level_rank", "cat_rank", "assembly_accession"]
        )
        genera_pools[g] = list(pool.index)

    import random
    rng = random.Random(seed)
    rng_genera = list(genera)
    rng.shuffle(rng_genera)

    round_robin = rng_genera.copy()
    while len(selected_indices) < n_genomes and round_robin:
        next_round = []
        for g in round_robin:
            if len(selected_indices) >= n_genomes:
                break
            pool = genera_pools.get(g, [])
            if pool:
                selected_indices.append(pool.pop(0))
                if pool:
                    next_round.append(g)
                genera_pools[g] = pool
        round_robin = next_round

    selection = candidates.loc[selected_indices].copy()
    selection["genome_id"] = selection["assembly_accession"].apply(
        lambda acc: f"{prefix}{acc.replace('_', '-').replace('.', '-')}"
    )
    selection["taxonomy_lookup"] = selection.apply(
        lambda row: "|".join([
            str(row["phylum"]), str(row["class"]), str(row["order_name"]),
            str(row["family"]), str(row["genus"]), str(row["species"]),
        ]),
        axis=1,
    )

    print(f"Selected {len(selection)} genomes across {selection['genus'].nunique()} genera, "
          f"{selection['family'].nunique()} families")
    return selection


def download_genomes(selection: pd.DataFrame, output_dir: Path) -> None:
    """Download FNA files from NCBI FTP."""
    output_dir.mkdir(parents=True, exist_ok=True)
    existing = {f.stem for f in output_dir.glob("*.fna")}

    for _, row in selection.iterrows():
        genome_id = row["genome_id"]
        if genome_id in existing:
            continue

        ftp_path = row["ftp_path"]
        base_name = ftp_path.rsplit("/", 1)[-1]
        url = f"{ftp_path}/{base_name}_genomic.fna.gz"

        out_path = output_dir / f"{genome_id}.fna"
        print(f"  Downloading {genome_id}...")

        try:
            import urllib.request
            import gzip
            import ssl

            for attempt in range(3):
                try:
                    ctx = ssl.create_default_context()
                    req = urllib.request.Request(url)
                    with urllib.request.urlopen(req, context=ctx, timeout=120) as resp:
                        compressed = resp.read()
                    break
                except Exception:
                    if attempt < 2:
                        ctx = ssl._create_unverified_context()
                        try:
                            with urllib.request.urlopen(req, context=ctx, timeout=120) as resp:
                                compressed = resp.read()
                            break
                        except Exception:
                            continue
                    raise

            raw = gzip.decompress(compressed)
            lines = raw.decode("utf-8", errors="replace").splitlines()
            rewritten = []
            for line in lines:
                if line.startswith(">"):
                    contig_header = line[1:].split()[0]
                    rewritten.append(f">{genome_id}|{contig_header}")
                else:
                    rewritten.append(line)
            out_path.write_text("\n".join(rewritten) + "\n")

        except Exception as e:
            print(f"  FAILED {genome_id}: {e}")
            if out_path.exists():
                out_path.unlink()


def write_metadata(selection: pd.DataFrame, output_dir: Path) -> None:
    """Write lookup and taxonomy files."""
    lookup_path = output_dir / "lookup.tsv"
    taxonomy_path = output_dir / "taxonomy.tsv"

    with open(lookup_path, "w") as f:
        for _, row in selection.iterrows():
            f.write(f"{row['genome_id']}\t{row['taxonomy_lookup']}\n")

    selection[["genome_id", "assembly_accession", "organism_name", "taxonomy_source",
               "phylum", "class", "order_name", "family", "genus", "species"]].to_csv(
        taxonomy_path, sep="\t", index=False
    )
    print(f"  Wrote {lookup_path} and {taxonomy_path}")


def step_prepare() -> None:
    """Step 1: Select and download genomes for both lineages."""
    for lineage_key, cfg in LINEAGES.items():
        print(f"\n{'='*80}")
        print(f"Preparing {cfg['description']}")
        print(f"{'='*80}")

        source_dir = RUNS_DIR / "source_genomes" / lineage_key
        fna_dir = source_dir / "fna"

        selection = select_diverse_genomes(
            TAXONOMY_DB,
            cfg["query_filter"],
            cfg["n_genomes"],
            cfg["prefix"],
        )

        download_genomes(selection, fna_dir)
        write_metadata(selection, source_dir)

        n_downloaded = len(list(fna_dir.glob("*.fna")))
        print(f"  {n_downloaded}/{cfg['n_genomes']} genomes downloaded")


def _ensure_shared_staging(lineage_key: str) -> Path:
    """Ensure the shared staging directory exists for a lineage.

    All panels for the same lineage share one staging dir (gene calling,
    hmmsearch, alignment, tree building for all 100 genomes). This avoids
    repeating the expensive staging step for each scope/seed combo.

    Returns the path to the shared staging directory.
    """
    shared_stage = RUNS_DIR / "benchmarks" / lineage_key / "_shared_staging"
    if (shared_stage / "proteomes").exists() and (shared_stage / "tree.nwk").exists():
        return shared_stage

    # Check if any existing panel already has a completed staging dir
    lineage_dir = RUNS_DIR / "benchmarks" / lineage_key
    if lineage_dir.exists():
        for candidate in sorted(lineage_dir.iterdir()):
            stage_dir = candidate / "stage_truth_full_model_clean"
            if stage_dir.is_dir() and (stage_dir / "tree.nwk").exists():
                # Copy this staging as the shared one
                print(f"  Reusing staging from {candidate.name}")
                if shared_stage.exists():
                    shutil.rmtree(shared_stage)
                shutil.copytree(stage_dir, shared_stage)
                return shared_stage

    # No existing staging found — will be created by first generate call
    return shared_stage


def _link_staging_into_panel(shared_stage: Path, panel_dir: Path) -> None:
    """Symlink the shared staging dir into a panel directory."""
    panel_dir.mkdir(parents=True, exist_ok=True)
    dest = panel_dir / "stage_truth_full_model_clean"
    if dest.exists() or dest.is_symlink():
        if dest.is_symlink():
            dest.unlink()
        elif dest.is_dir():
            shutil.rmtree(dest)
    dest.symlink_to(shared_stage.resolve())


def run_sgtree(
    genomedir: Path, modeldir: Path, outdir: Path, *, profile: dict
) -> None:
    """Run SGTree with the given cleanup profile."""
    cmd = [
        sys.executable, "-m", "sgtree",
        str(genomedir), str(modeldir),
        "--num_cpus", str(NUM_CPUS),
        "--percent_models", "100",
        "--save_dir", str(outdir),
        "--selection_mode", profile["selection_mode"],
        "--selection_max_rounds", "5",
        "--selection_global_rounds", str(profile["selection_global_rounds"]),
        "--keep_intermediates", "yes",
    ]
    if profile["marker_selection"]:
        cmd.extend(["--marker_selection", "yes"])
    if profile["singles"]:
        cmd.extend(["--singles", "yes"])
    cmd.extend(["--singles_mode", profile["singles_mode"]])
    subprocess.run(cmd, check=False)


def _clean_intermediates(result_dir: Path) -> None:
    """Remove bulky intermediate files from a completed run to save disk."""
    for subdir_name in ["alignments", "hmmsearch", "proteomes", "temp"]:
        subdir = result_dir / subdir_name
        if subdir.is_dir():
            shutil.rmtree(subdir)
    # Also clean rewritten proteome dirs for non-solo variants
    variant_proteomes = result_dir.parent / "proteomes"
    if variant_proteomes.is_dir():
        shutil.rmtree(variant_proteomes)


def rewrite_proteome_merged_contam(src_dir: Path, dst_dir: Path) -> None:
    """Variant A: merge all contaminant contigs per genome into one shared contig."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    for faa in sorted(src_dir.glob("*.faa")):
        genome = faa.stem
        merged_contig = f"contig__contam_merged__{genome}"
        lines = faa.read_text().splitlines()
        out_lines = []
        for line in lines:
            if line.startswith(">") and "contig__contam__" in line:
                parts = line[1:].split("|")
                parts[1] = merged_contig
                out_lines.append(">" + "|".join(parts))
            else:
                out_lines.append(line)
        (dst_dir / faa.name).write_text("\n".join(out_lines) + "\n")


def rewrite_proteome_native_contig(
    src_dir: Path, dst_dir: Path, events: pd.DataFrame
) -> None:
    """Variant B: move each contaminant marker onto its native target contig."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    contam_to_native: dict[tuple[str, str], str] = {}
    for row in events.itertuples(index=False):
        if hasattr(row, "contaminant_contig_id") and hasattr(row, "native_contig_id"):
            contam_contig = str(row.contaminant_contig_id)
            native_contig = str(row.native_contig_id)
            genome = str(row.recipient_genome)
            contam_to_native[(genome, contam_contig)] = native_contig

    for faa in sorted(src_dir.glob("*.faa")):
        genome = faa.stem
        lines = faa.read_text().splitlines()
        out_lines = []
        for line in lines:
            if line.startswith(">") and "contig__contam__" in line:
                parts = line[1:].split("|")
                old_contig = parts[1]
                native = contam_to_native.get((genome, old_contig))
                if native is not None:
                    parts[1] = native
                out_lines.append(">" + "|".join(parts))
            else:
                out_lines.append(line)
        (dst_dir / faa.name).write_text("\n".join(out_lines) + "\n")


def _enrich_events_with_contig_ids(events: pd.DataFrame) -> pd.DataFrame:
    """Add contaminant_contig_id and native_contig_id if missing."""
    if "contaminant_contig_id" in events.columns:
        return events
    from sgtree.id_schema import parse_sequence_id
    contam_contigs = []
    native_contigs = []
    for _, row in events.iterrows():
        _, contam_contig, _ = parse_sequence_id(str(row["contaminant_record_id"]))
        _, native_contig, _ = parse_sequence_id(str(row["native_record_id"]))
        contam_contigs.append(contam_contig)
        native_contigs.append(native_contig)
    events = events.copy()
    events["contaminant_contig_id"] = contam_contigs
    events["native_contig_id"] = native_contigs
    return events


def step_generate() -> None:
    """Step 2: Generate benchmark panels for all scopes and seeds."""
    sys.path.insert(0, str(PROJECT_ROOT))
    from sgtree import benchmarks as bm
    original_scenarios = bm.DEFAULT_SCENARIOS
    bm.DEFAULT_SCENARIOS = DOUBLED_SCENARIOS

    try:
        for lineage_key, cfg in LINEAGES.items():
            print(f"\n{'='*80}")
            print(f"Staging {cfg['description']}")
            print(f"{'='*80}")

            source_dir = RUNS_DIR / "source_genomes" / lineage_key / "fna"
            shared_stage = _ensure_shared_staging(lineage_key)

            for scope in SCOPES:
                for seed in SEEDS:
                    pname = panel_name(scope, seed)
                    panel_dir = RUNS_DIR / "benchmarks" / lineage_key / pname

                    if (panel_dir / "benchmark_manifest.json").exists():
                        print(f"  SKIP {lineage_key}/{pname}: already exists")
                        continue

                    print(f"\n  GENERATE {lineage_key}/{pname}")

                    # Link shared staging so it's not re-run
                    if shared_stage.exists() and (shared_stage / "tree.nwk").exists():
                        _link_staging_into_panel(shared_stage, panel_dir)

                    if scope == "phylum":
                        donor_key = CROSS_LINEAGE_DONOR[lineage_key]
                        donor_cfg = LINEAGES[donor_key]
                        donor_source = RUNS_DIR / "source_genomes" / donor_key / "fna"
                        donor_label = donor_cfg["label"]
                        # Link the donor lineage's shared staging to avoid re-staging
                        donor_shared = _ensure_shared_staging(donor_key)
                        if donor_shared.exists() and (donor_shared / "tree.nwk").exists():
                            donor_stage_dest = panel_dir / "stage_donor_full_model_clean"
                            panel_dir.mkdir(parents=True, exist_ok=True)
                            if donor_stage_dest.exists() or donor_stage_dest.is_symlink():
                                if donor_stage_dest.is_symlink():
                                    donor_stage_dest.unlink()
                                elif donor_stage_dest.is_dir():
                                    shutil.rmtree(donor_stage_dest)
                            donor_stage_dest.symlink_to(donor_shared.resolve())
                    else:
                        donor_source = None
                        donor_label = None

                    bm.generate_taxonomic_benchmark_dataset(
                        truth_source_dir=source_dir,
                        donor_source_dir=donor_source,
                        models_path=MODELS_PATH,
                        outdir=panel_dir,
                        taxonomic_scope=scope,
                        taxonomy_db_path=TAXONOMY_DB,
                        lineage_label=cfg["label"],
                        donor_lineage_label=donor_label,
                        n_genomes=cfg["n_genomes"],
                        n_markers=10,
                        min_marker_presence_fraction=0.8,
                        seed=seed,
                        num_cpus=NUM_CPUS,
                    )

                    # After first successful generate, ensure shared staging exists
                    if not shared_stage.exists() or not (shared_stage / "tree.nwk").exists():
                        real_stage = panel_dir / "stage_truth_full_model_clean"
                        if real_stage.is_dir() and not real_stage.is_symlink():
                            shutil.copytree(real_stage, shared_stage)
    finally:
        bm.DEFAULT_SCENARIOS = original_scenarios


def step_run() -> None:
    """Step 3: Run SGTree on all 3 variants for each panel."""
    for lineage_key, cfg, scope, seed, panel_dir in iter_panels():
        pname = panel_name(scope, seed)
        if not (panel_dir / "benchmark_manifest.json").exists():
            print(f"SKIP {lineage_key}/{pname}: no benchmark manifest")
            continue

        model_dir = panel_dir / "truth_markers.hmm"
        scenarios_dir = panel_dir / "scenarios"

        for scenario_dir in sorted(scenarios_dir.iterdir()):
            if not scenario_dir.is_dir():
                continue
            scenario_name = scenario_dir.name
            if scenario_name not in CLEANUP_PROFILES:
                continue

            profile = CLEANUP_PROFILES[scenario_name]
            profile_name = profile["name"]
            src_proteomes = scenario_dir / "proteomes"
            events_path = scenario_dir / "events.tsv"

            if not events_path.exists():
                print(f"  SKIP {lineage_key}/{pname}/{scenario_name}: no events.tsv")
                continue

            events = pd.read_csv(events_path, sep="\t")
            events = _enrich_events_with_contig_ids(events)

            for variant_label, rewrite_fn in [
                ("solo_contig", None),
                ("merged_contam", rewrite_proteome_merged_contam),
                ("native_contig", rewrite_proteome_native_contig),
            ]:
                result_dir = (
                    RUNS_DIR / "results" / lineage_key / pname / scenario_name
                    / variant_label / f"{scenario_name}__{profile_name}"
                )

                if (result_dir / "tree_final.nwk").exists() or (result_dir / "tree.nwk").exists():
                    print(f"  SKIP {lineage_key}/{pname}/{scenario_name}/{variant_label}: done")
                    continue

                if variant_label == "solo_contig":
                    proteome_dir = src_proteomes
                else:
                    proteome_dir = (
                        RUNS_DIR / "results" / lineage_key / pname / scenario_name
                        / variant_label / "proteomes"
                    )
                    if variant_label == "merged_contam":
                        rewrite_fn(src_proteomes, proteome_dir)
                    else:
                        rewrite_fn(src_proteomes, proteome_dir, events)

                print(f"  RUN  {lineage_key}/{pname}/{scenario_name}/{variant_label}")
                run_sgtree(proteome_dir, model_dir, result_dir, profile=profile)

                has_tree = (result_dir / "tree_final.nwk").exists()
                status = "OK" if has_tree else "INCOMPLETE"
                print(f"  {status} {lineage_key}/{pname}/{scenario_name}/{variant_label}")

                # Clean intermediates to save disk space
                _clean_intermediates(result_dir)


def step_evaluate() -> None:
    """Step 4: Evaluate all results."""
    sys.path.insert(0, str(PROJECT_ROOT))
    from sgtree.benchmarks import _rf_norm, _parse_singleton_result_file

    all_results = []

    for lineage_key, cfg, scope, seed, panel_dir in iter_panels():
        pname = panel_name(scope, seed)
        if not (panel_dir / "benchmark_manifest.json").exists():
            continue

        scenarios_dir = panel_dir / "scenarios"
        truth_tree = panel_dir / "truth_run" / "tree.nwk"

        for scenario_dir in sorted(scenarios_dir.iterdir()):
            if not scenario_dir.is_dir():
                continue
            scenario_name = scenario_dir.name
            if scenario_name not in CLEANUP_PROFILES:
                continue

            events_path = scenario_dir / "events.tsv"
            if not events_path.exists():
                continue
            events = pd.read_csv(events_path, sep="\t")
            profile_name = CLEANUP_PROFILES[scenario_name]["name"]

            dup_events = events[events["event_type"].isin(["duplicate", "triplicate"])]
            rep_events = events[events["event_type"] == "replacement"]

            for variant in ["solo_contig", "merged_contam", "native_contig"]:
                result_dir = (
                    RUNS_DIR / "results" / lineage_key / pname / scenario_name
                    / variant / f"{scenario_name}__{profile_name}"
                )

                final_tree = result_dir / "tree_final.nwk"
                initial_tree = result_dir / "tree.nwk"
                if not initial_tree.exists():
                    initial_tree = result_dir / "tree_round_1.nwk"
                if not final_tree.exists():
                    if initial_tree.exists():
                        final_tree = initial_tree
                    else:
                        continue
                if not final_tree.exists():
                    continue

                try:
                    initial_rf = _rf_norm(truth_tree, initial_tree)
                    final_rf = _rf_norm(truth_tree, final_tree)
                except Exception:
                    continue

                # RF status for duplicate removal
                rf_status: dict[tuple[str, str], str] = {}
                rf_path = result_dir / "marker_selection_rf_values.txt"
                if rf_path.exists():
                    for line in rf_path.read_text().splitlines():
                        line = line.strip()
                        if not line or line.startswith("#") or line.startswith("Protein"):
                            continue
                        parts = line.split()
                        if len(parts) >= 4 and parts[-1] in ("Kept", "Removed"):
                            protein_id = parts[0].replace("/", "|")
                            marker = parts[1]
                            rf_status[(protein_id, marker)] = parts[-1]

                rf_by_gene: dict[tuple[str, str], str] = {}
                for (pid, marker), st in rf_status.items():
                    gene = pid.split("|")[-1] if "|" in pid else pid
                    rf_by_gene[(gene, marker)] = st

                dup_removed = 0
                for row in dup_events.itertuples(index=False):
                    cid = str(row.contaminant_record_id).replace("/", "|")
                    if rf_status.get((cid, row.marker)) == "Removed":
                        dup_removed += 1
                    else:
                        gene = cid.split("|")[-1] if "|" in cid else cid
                        if rf_by_gene.get((gene, row.marker)) == "Removed":
                            dup_removed += 1

                # Singleton calls
                removed_dir = result_dir / "removed"
                calls = []
                if removed_dir.exists():
                    for marker_file in sorted(removed_dir.iterdir()):
                        if marker_file.is_file():
                            calls.extend(_parse_singleton_result_file(marker_file))

                contam_calls = [c for c in calls if c.get("is_contaminant_leaf", False)]
                native_calls = [c for c in calls if not c.get("is_contaminant_leaf", False)]

                cc = sum(1 for c in contam_calls if c.get("singleton_class") == "contamination_candidate")
                ch = sum(1 for c in contam_calls if c.get("singleton_class") == "hgt_candidate")
                cp = sum(1 for c in contam_calls if c.get("decision") == "pruned")
                nc = sum(1 for c in native_calls if c.get("singleton_class") == "contamination_candidate")
                nh = sum(1 for c in native_calls if c.get("singleton_class") == "hgt_candidate")

                total_contam = len(dup_events) + len(rep_events)
                total_removed = dup_removed + cp

                precision = cc / (cc + nc) if (cc + nc) > 0 else float("nan")

                all_results.append({
                    "lineage": lineage_key,
                    "scope": scope,
                    "seed": seed,
                    "panel": pname,
                    "scenario": scenario_name,
                    "variant": variant,
                    "initial_rf": round(initial_rf, 4),
                    "final_rf": round(final_rf, 4),
                    "rf_delta": round(initial_rf - final_rf, 4),
                    "dup_events": len(dup_events),
                    "dup_removed": dup_removed,
                    "rep_events": len(rep_events),
                    "contam_pruned": cp,
                    "total_contam": total_contam,
                    "total_removed": total_removed,
                    "removal_pct": round(total_removed / total_contam * 100, 1) if total_contam > 0 else 0,
                    "contam_as_contam": cc,
                    "contam_as_hgt": ch,
                    "contam_detect_pct": round(cc / (cc + ch) * 100, 1) if (cc + ch) > 0 else float("nan"),
                    "native_as_hgt": nh,
                    "native_as_contam": nc,
                    "native_keep_pct": round(nh / (nh + nc) * 100, 1) if (nh + nc) > 0 else float("nan"),
                    "contam_precision_pct": round(precision * 100, 1),
                })

    if not all_results:
        print("No results found!")
        return

    df = pd.DataFrame(all_results)
    out_path = RUNS_DIR / "proof_of_concept_results.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nFull results: {len(df)} rows saved to {out_path}")

    # Summary by scope × variant (pooled across lineages, seeds, scenarios)
    print(f"\n{'='*80}")
    print("SUMMARY BY SCOPE x VARIANT (pooled across lineages, seeds, scenarios)")
    print(f"{'='*80}")
    for scope in SCOPES:
        for variant in ["solo_contig", "merged_contam", "native_contig"]:
            sub = df[(df["scope"] == scope) & (df["variant"] == variant)]
            if sub.empty:
                continue
            _print_summary(f"{scope} / {variant}", sub)

    # Summary by lineage × scope × variant
    print(f"\n{'='*80}")
    print("SUMMARY BY LINEAGE x SCOPE x VARIANT")
    print(f"{'='*80}")
    for lineage in LINEAGES:
        for scope in SCOPES:
            for variant in ["solo_contig", "merged_contam", "native_contig"]:
                sub = df[
                    (df["lineage"] == lineage)
                    & (df["scope"] == scope)
                    & (df["variant"] == variant)
                ]
                if sub.empty:
                    continue
                _print_summary(f"{lineage} / {scope} / {variant}", sub)


def _print_summary(label: str, sub: pd.DataFrame) -> None:
    total_rem = sub["total_removed"].sum()
    total_add = sub["total_contam"].sum()
    cc = sub["contam_as_contam"].sum()
    ch = sub["contam_as_hgt"].sum()
    nc = sub["native_as_contam"].sum()
    nh = sub["native_as_hgt"].sum()
    print(f"\n--- {label} ({len(sub)} runs) ---")
    print(f"  Mean RF delta: {sub['rf_delta'].mean():.4f} +/- {sub['rf_delta'].std():.4f}")
    if total_add > 0:
        print(f"  Contaminants removed: {total_rem}/{total_add} ({total_rem/total_add*100:.1f}%)")
    print(f"  Dup removed: {sub['dup_removed'].sum()}/{sub['dup_events'].sum()}")
    if (cc + ch) > 0:
        print(f"  Contam singleton detection: {cc}/{cc+ch} ({cc/(cc+ch)*100:.1f}%)")
    if (nh + nc) > 0:
        print(f"  Native keep rate: {nh}/{nh+nc} ({nh/(nh+nc)*100:.1f}%)")
    if (cc + nc) > 0:
        print(f"  Precision: {cc}/{cc+nc} ({cc/(cc+nc)*100:.1f}%)")


def main() -> None:
    parser = argparse.ArgumentParser(description="Contig variant proof-of-concept benchmark")
    parser.add_argument("--step", required=True,
                        choices=["prepare", "generate", "run", "evaluate", "all"])
    args = parser.parse_args()

    if args.step in ("prepare", "all"):
        step_prepare()
    if args.step in ("generate", "all"):
        step_generate()
    if args.step in ("run", "all"):
        step_run()
    if args.step in ("evaluate", "all"):
        step_evaluate()


if __name__ == "__main__":
    main()
