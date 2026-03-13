from __future__ import annotations

import os
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from sgtree import ani as ani_core
from sgtree.config import Config
from sgtree.id_schema import sanitize_token
from sgtree.input_stage import NUCLEOTIDE_EXTENSIONS, detect_input_format


QUERY_MANIFEST_NAME = "query_manifest.tsv"
REF_MANIFEST_NAME = "ref_manifest.tsv"


def _fasta_stats(path: Path) -> tuple[int, int]:
    contigs = 0
    total_bases = 0
    with path.open() as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contigs += 1
            total_bases += len(record.seq)
    return contigs, total_bases


def _list_nucleotide_files(input_dir: str) -> list[Path]:
    directory = Path(input_dir)
    if not directory.is_dir():
        raise ValueError(f"ANI clustering currently requires a genome directory, got: {input_dir}")
    paths: list[Path] = []
    for ext in NUCLEOTIDE_EXTENSIONS:
        paths.extend(sorted(directory.glob(f"*{ext}")))
    unique = sorted({path.resolve() for path in paths if path.is_file()}, key=lambda path: path.name)
    if not unique:
        raise FileNotFoundError(f"No nucleotide genome files found under {input_dir}")
    return unique


def _discover_inputs(input_dir: str, source: str) -> list[dict[str, object]]:
    if detect_input_format(input_dir) != "fna":
        raise ValueError(
            f"ANI clustering requires nucleotide genomes for {source}, but {input_dir} did not resolve to FNA input"
        )
    rows: list[dict[str, object]] = []
    for path in _list_nucleotide_files(input_dir):
        contigs, total_bases = _fasta_stats(path)
        genome_id = sanitize_token(path.stem, path.stem)
        rows.append(
            {
                "genome_id": genome_id,
                "input_format": "fna",
                "source_file": str(path),
                "assembly_path": str(path),
                "staged_proteome_path": "",
                "contigs": contigs,
                "total_bases": total_bases,
                "source_role": source,
            }
        )
    return rows


def _write_manifest(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    if frame.empty:
        frame = pd.DataFrame(
            columns=[
                "genome_id",
                "input_format",
                "source_file",
                "assembly_path",
                "staged_proteome_path",
                "contigs",
                "total_bases",
                "source_role",
            ]
        )
    frame.sort_values(["genome_id"]).to_csv(path, sep="\t", index=False)


def _manifest_paths(cfg: Config) -> tuple[Path, Path | None]:
    query_manifest = Path(cfg.ani_dir) / QUERY_MANIFEST_NAME
    ref_manifest = Path(cfg.ani_dir) / REF_MANIFEST_NAME
    if cfg.original_ref or cfg.ref:
        return query_manifest, ref_manifest
    return query_manifest, None


def prepare_ani_cluster_inputs(cfg: Config) -> pd.DataFrame:
    os.makedirs(cfg.ani_dir, exist_ok=True)
    query_rows = _discover_inputs(cfg.original_genomedir or cfg.genomedir, "query")
    ref_source_dir = cfg.original_ref or cfg.ref
    ref_rows = _discover_inputs(ref_source_dir, "ref") if ref_source_dir else []

    query_manifest, ref_manifest = _manifest_paths(cfg)
    _write_manifest(query_rows, query_manifest)
    if ref_manifest is not None:
        _write_manifest(ref_rows, ref_manifest)

    ani_core.run_ani_clustering(
        query_manifest=query_manifest,
        ref_manifest=ref_manifest,
        outdir=cfg.ani_dir,
        ani_threshold=cfg.ani_threshold / 100.0 if cfg.ani_threshold > 1.0 else cfg.ani_threshold,
        inflation=cfg.ani_mcl_inflation,
        backend=cfg.ani_backend,
        cpus=cfg.num_cpus,
    )

    members = pd.read_csv(cfg.ani_cluster_members_path, sep="\t")
    source_by_genome = {
        str(row["genome_id"]): str(row["source_role"])
        for row in query_rows + ref_rows
    }
    path_by_genome = {
        str(row["genome_id"]): str(row["assembly_path"])
        for row in query_rows + ref_rows
    }
    members = members.copy()
    members["source"] = members["source_role"].astype(str)
    members["genome_path"] = members["genome_id"].astype(str).map(path_by_genome)
    members["representative_genome_id"] = members["representative_genome"].astype(str)
    members["representative_source"] = members["representative_genome"].astype(str).map(source_by_genome)
    members["is_representative"] = members["kept_for_species_tree"].astype(str).str.lower() == "yes"
    members.to_csv(cfg.ani_cluster_members_path, sep="\t", index=False)

    for directory in (cfg.ani_selected_query_dir, cfg.ani_selected_ref_dir):
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)

    kept = members[members["is_representative"]].copy()
    for row in kept.itertuples(index=False):
        target_dir = cfg.ani_selected_ref_dir if str(row.source) == "ref" else cfg.ani_selected_query_dir
        target = Path(target_dir) / Path(str(row.genome_path)).name
        if target.exists() or target.is_symlink():
            target.unlink()
        try:
            target.symlink_to(Path(str(row.genome_path)).resolve())
        except OSError:
            shutil.copy2(str(row.genome_path), target)

    cfg.genomedir = cfg.ani_selected_query_dir
    cfg.ref = cfg.ani_selected_ref_dir if ref_rows else None
    cfg.ref_concat = os.path.join(cfg.outdir, "reference_cache")

    print(
        f"-... ANI clustering retained {len(kept)} representative genomes from {len(members)} total inputs "
        f"at >= {cfg.ani_threshold:.2f}% ANI"
    )
    return members


def build_cluster_snp_trees(cfg: Config) -> pd.DataFrame:
    if not Path(cfg.ani_cluster_members_path).exists():
        return pd.DataFrame()
    query_manifest, ref_manifest = _manifest_paths(cfg)
    summary_path = ani_core.build_snp_trees(
        clusters_path=cfg.ani_cluster_members_path,
        query_manifest=query_manifest,
        ref_manifest=ref_manifest,
        outdir=cfg.snp_trees_dir,
        min_cluster_size=cfg.snp_tree_min_cluster_size,
        cpus=cfg.num_cpus,
    )
    return pd.read_csv(summary_path, sep="\t")
