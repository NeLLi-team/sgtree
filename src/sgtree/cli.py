"""Command-line entry point for the SGTree analysis pipeline."""

from __future__ import annotations

import argparse
import datetime
import logging
import os
import shutil
import time
import traceback
from pathlib import Path
from typing import TYPE_CHECKING

from sgtree import (
    align,
    ani_clustering,
    cleanup,
    duplicates,
    extract,
    marker_selection,
    phylogeny,
    reference,
    search,
    sgtree_logging,
    supermatrix,
)
from sgtree.config import Config

if TYPE_CHECKING:
    import pandas as pd

os.environ["QT_QPA_PLATFORM"] = "offscreen"

logger = logging.getLogger("sgtree")

_RUN_MARKER_NAME = cleanup.RUN_MARKER_NAME
_RUN_MARKER_CONTENT = cleanup.RUN_MARKER_CONTENT


def _local_now() -> datetime.datetime:
    """Return local wall-clock time without timezone text in legacy outputs."""
    return datetime.datetime.now(datetime.UTC).astimezone().replace(tzinfo=None)


def _parse_bool(*, value: str | bool, flag: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"yes", "true", "1"}:
        return True
    if normalized in {"no", "false", "0"}:
        return False
    raise ValueError(f"{flag} expects one of: yes/no, true/false, 1/0")


def _resolve_required_path_arg(
    parser: argparse.ArgumentParser,
    *,
    positional: str | None,
    flagged: str | None,
    name: str,
) -> str:
    if flagged is None and positional is None:
        parser.error(f"the following arguments are required: {name}")
    if flagged is not None and positional is not None and flagged != positional:
        parser.error(
            f"{name} specified with conflicting values: "
            f"positional={positional!r}, flag={flagged!r}"
        )
    value = flagged if flagged is not None else positional
    assert value is not None
    return value


def _normalize_tree_method(value: str) -> str:
    normalized = value.strip().lower()
    if normalized == "fasttree":
        return "veryfasttree"
    return normalized


def _normalize_cli_path(value: str) -> str:
    """Remove redundant separators without collapsing the filesystem root."""
    return os.path.normpath(value)


def _path_label(value: str) -> str:
    return Path(_normalize_cli_path(value)).name or "root"


def _run_marker_path(cfg: Config) -> Path:
    return Path(cfg.outdir) / _RUN_MARKER_NAME


def _has_run_marker(cfg: Config) -> bool:
    marker = _run_marker_path(cfg)
    if marker.is_symlink() or not marker.is_file():
        return False
    return marker.read_text(encoding="utf-8") == _RUN_MARKER_CONTENT


def _mark_run_directory(cfg: Config) -> None:
    marker = _run_marker_path(cfg)
    if _has_run_marker(cfg):
        return
    if marker.exists() or marker.is_symlink():
        raise FileExistsError(f"refusing to overwrite reserved run marker: {marker}")
    with marker.open("x", encoding="utf-8") as handle:
        handle.write(_RUN_MARKER_CONTENT)


def _clean_previous_run(cfg: Config) -> None:
    """Remove the outputs of an earlier run from cfg.outdir.

    Unreserved files survive. Reserved paths in an unmarked, incomplete output
    directory cause an error because SGTree cannot establish their ownership.
    """
    output_dir = Path(cfg.outdir)
    generated_paths = [
        path for path in output_dir.iterdir() if cleanup.is_generated_entry(path.name)
    ]
    has_legacy_tree = any(
        (output_dir / name).exists() for name in ("tree.nwk", "tree_final.nwk")
    )
    if not _has_run_marker(cfg) and not has_legacy_tree:
        if generated_paths:
            names = ", ".join(sorted(path.name for path in generated_paths))
            raise FileExistsError(
                f"output directory contains reserved SGTree paths but no "
                f"{_RUN_MARKER_NAME} ownership marker: {names}"
            )
        return

    linked_paths = [path for path in generated_paths if path.is_symlink()]
    if linked_paths:
        raise FileExistsError(
            f"refusing to remove linked generated path: {linked_paths[0]}"
        )

    for path in generated_paths:
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def _build_parser() -> argparse.ArgumentParser:
    """Build the SGTree command-line parser."""
    parser = argparse.ArgumentParser(
        description="SGTree - Species tree from marker gene phylogenies",
        epilog=(
            "Provide the required inputs either as positional arguments or with "
            "--genomedir/--modeldir."
        ),
    )

    parser.add_argument(
        "genomedir",
        nargs="?",
        type=str,
        help=(
            "directory containing .faa proteomes or .fna assemblies "
            "(or a concatenated fasta)"
        ),
    )
    parser.add_argument(
        "modeldir",
        nargs="?",
        type=str,
        help=(
            "path to marker-set .hmm file (legacy: directory of per-marker .hmm files)"
        ),
    )
    parser.add_argument(
        "--genomedir",
        dest="genomedir_flag",
        metavar="GENOMEDIR",
        type=str,
        default=None,
        help=(
            "directory containing .faa proteomes or .fna assemblies "
            "(flag alias for the required input)"
        ),
    )
    parser.add_argument(
        "--modeldir",
        dest="modeldir_flag",
        metavar="MODELDIR",
        type=str,
        default=None,
        help=("path to marker-set .hmm file (flag alias for the required input)"),
    )
    parser.add_argument(
        "--ref_concat",
        type=str,
        default=None,
        help="path to store reference directory concat files",
    )
    parser.add_argument("--num_cpus", type=int, default=8, help="number of CPUs to use")
    parser.add_argument(
        "--marker_selection",
        type=str,
        default="no",
        help="run marker selection (yes/no)",
    )
    parser.add_argument(
        "--ref", type=str, default=None, help="reference genomes directory"
    )
    parser.add_argument(
        "--percent_models",
        type=int,
        default=10,
        help="minimum percentage of models a genome must have",
    )
    parser.add_argument(
        "--save_dir",
        "--outdir",
        dest="save_dir",
        type=str,
        default=None,
        help="output directory name",
    )
    parser.add_argument(
        "--singles",
        type=str,
        default="no",
        help=(
            "run mode-dependent singleton diagnostics or experimental pruning (yes/no)"
        ),
    )
    parser.add_argument(
        "--singles-mode",
        "--singles_mode",
        dest="singles_mode",
        type=str,
        default="delta_rf",
        choices=[
            "delta_rf",
            "topoknn",
            "hybrid",
            "composite",
            "contig_consensus",
            "recipient_consensus",
            "neighbor_clade",
            "neighbor_ml",
            "gcp",
            "loo_profile",
        ],
        help="singleton mode; loo_profile reports evidence without pruning",
    )
    parser.add_argument(
        "--lflt",
        type=int,
        default=0,
        help="remove sequences shorter than N%% of median length",
    )
    parser.add_argument(
        "--num_nei", type=int, default=0, help="singleton neighborhood size (0=auto)"
    )
    parser.add_argument(
        "--singles_min_rfdist",
        type=float,
        default=0.25,
        help=(
            "minimum marker-tree/global-tree RF distance required to activate "
            "singleton filtering"
        ),
    )
    parser.add_argument(
        "--aln",
        type=str,
        default="hmmalign",
        choices=["hmmalign", "mafft", "mafft-linsi", "famsa"],
        help="alignment method: hmmalign, mafft, mafft-linsi, or famsa",
    )
    parser.add_argument(
        "--tree_method",
        type=str,
        default="veryfasttree",
        choices=["veryfasttree", "fasttree", "iqtree"],
        help=(
            "tree builder for species and marker trees: veryfasttree "
            "(legacy alias: fasttree) or iqtree"
        ),
    )
    parser.add_argument(
        "--iqtree_fast",
        type=str,
        default="yes",
        help="when --tree_method iqtree, use IQ-TREE -fast (yes/no)",
    )
    parser.add_argument(
        "--iqtree_model",
        type=str,
        default="LG+F+I+G4",
        help="IQ-TREE model string (used when --tree_method iqtree)",
    )
    parser.add_argument(
        "--hmmsearch_cutoff",
        type=str,
        default="cut_ga",
        choices=["cut_ga", "cut_tc", "cut_nc", "evalue"],
        help="hmmsearch threshold mode",
    )
    parser.add_argument(
        "--hmmsearch_evalue",
        type=float,
        default=1e-5,
        help="hmmsearch E-value threshold when --hmmsearch_cutoff evalue",
    )
    parser.add_argument(
        "--selection_mode",
        type=str,
        default="coordinate",
        choices=["legacy", "coordinate"],
        help="duplicate-resolution mode for marker selection",
    )
    parser.add_argument(
        "--selection_max_rounds",
        type=int,
        default=5,
        help="max coordinate-descent rounds for --selection_mode coordinate",
    )
    parser.add_argument(
        "--selection_global_rounds",
        type=int,
        default=1,
        help="max guide-tree rebuild rounds for marker selection",
    )
    parser.add_argument(
        "--lock_references",
        type=str,
        default="no",
        help="keep reference duplicate resolution score-locked (yes/no)",
    )
    parser.add_argument(
        "--max_sdup",
        type=int,
        default=-1,
        help="max copies allowed for any single marker per genome (-1 disables)",
    )
    parser.add_argument(
        "--max_dupl",
        type=float,
        default=-1.0,
        help="max fraction of markers allowed in duplicate per genome (-1 disables)",
    )
    parser.add_argument(
        "--is_ref", type=str, default="no", help="internal flag, not for user use"
    )
    parser.add_argument(
        "--keep_intermediates",
        type=str,
        default="no",
        help="keep intermediate directories/files instead of archiving them (yes/no)",
    )
    parser.add_argument(
        "--ani_cluster",
        type=str,
        default="no",
        help=(
            "collapse query+reference genomes by ANI before species-tree inference "
            "(yes/no)"
        ),
    )
    parser.add_argument(
        "--snp",
        type=str,
        default="no",
        help="build cluster-level SNP trees after ANI clustering (yes/no)",
    )
    parser.add_argument(
        "--ani_threshold",
        type=float,
        default=95.0,
        help="minimum ANI threshold for graph edges before MCL clustering",
    )
    parser.add_argument(
        "--ani_backend",
        type=str,
        default="auto",
        choices=["auto", "skani", "minimap2"],
        help="ANI backend preference",
    )
    parser.add_argument(
        "--ani_mcl_inflation",
        type=float,
        default=2.0,
        help="MCL inflation parameter used after ANI edge filtering",
    )
    parser.add_argument(
        "--snp_tree_min_cluster_size",
        type=int,
        default=3,
        help="build SNP trees only for ANI clusters with at least this many genomes",
    )

    return parser


def _validate_args(
    args: argparse.Namespace,
    *,
    ani_cluster: bool,
    snp: bool,
) -> None:
    validations = (
        (args.num_cpus >= 1, "--num_cpus must be >= 1"),
        (
            0 <= args.percent_models <= 100,
            "--percent_models must be between 0 and 100",
        ),
        (0 <= args.lflt <= 100, "--lflt must be between 0 and 100"),
        (args.max_sdup >= -1, "--max_sdup must be >= -1"),
        (
            args.max_dupl == -1.0 or 0.0 <= args.max_dupl <= 1.0,
            "--max_dupl must be between 0 and 1, or -1 to disable",
        ),
        (
            args.hmmsearch_cutoff != "evalue" or args.hmmsearch_evalue > 0,
            "--hmmsearch_evalue must be > 0 when cutoff mode is evalue",
        ),
        (args.num_nei >= 0, "--num_nei must be >= 0"),
        (
            0.0 <= args.singles_min_rfdist <= 1.0,
            "--singles_min_rfdist must be between 0 and 1",
        ),
        (args.selection_max_rounds >= 1, "--selection_max_rounds must be >= 1"),
        (
            args.selection_global_rounds >= 1,
            "--selection_global_rounds must be >= 1",
        ),
        (
            0.0 < args.ani_threshold <= 100.0,
            "--ani_threshold must be > 0 and <= 100",
        ),
        (args.ani_mcl_inflation > 0, "--ani_mcl_inflation must be > 0"),
        (
            args.snp_tree_min_cluster_size >= 2,
            "--snp_tree_min_cluster_size must be >= 2",
        ),
        (not snp or ani_cluster, "--snp requires --ani_cluster yes"),
    )
    for valid, message in validations:
        if not valid:
            raise ValueError(message)


def _config_from_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> Config:
    start_time = str(_local_now())
    genomedir = _normalize_cli_path(
        _resolve_required_path_arg(
            parser,
            positional=args.genomedir,
            flagged=args.genomedir_flag,
            name="genomedir",
        )
    )
    modeldir = _normalize_cli_path(
        _resolve_required_path_arg(
            parser,
            positional=args.modeldir,
            flagged=args.modeldir_flag,
            name="modeldir",
        )
    )

    if args.save_dir:
        outdir = _normalize_cli_path(args.save_dir)
    else:
        ref_part = _path_label(args.ref) if args.ref else "no_ref_directory"
        run_name = (
            f"SG_{_path_label(genomedir)}_{ref_part}_{_path_label(modeldir)}_"
            f"{start_time.replace(':', '-').replace(' ', '_').split('.')[0]}"
        )
        outdir = str(Path.cwd() / "runs" / "python" / run_name)

    ref_concat = args.ref_concat or str(Path.cwd() / "runs" / "reference_cache")
    ref_concat = _normalize_cli_path(ref_concat)
    iqtree_fast = _parse_bool(value=args.iqtree_fast, flag="--iqtree_fast")
    marker_selection_enabled = _parse_bool(
        value=args.marker_selection,
        flag="--marker_selection",
    )
    singles = _parse_bool(value=args.singles, flag="--singles")
    is_ref = _parse_bool(value=args.is_ref, flag="--is_ref")
    lock_references = _parse_bool(
        value=args.lock_references,
        flag="--lock_references",
    )
    keep_intermediates = _parse_bool(
        value=args.keep_intermediates,
        flag="--keep_intermediates",
    )
    ani_cluster = _parse_bool(value=args.ani_cluster, flag="--ani_cluster")
    snp = _parse_bool(value=args.snp, flag="--snp")
    _validate_args(args, ani_cluster=ani_cluster, snp=snp)

    return Config(
        genomedir=genomedir,
        modeldir=modeldir,
        outdir=outdir,
        num_cpus=args.num_cpus,
        percent_models=args.percent_models,
        input_format="auto",
        lflt_fraction=float(args.lflt) / 100,
        aln_method=args.aln,
        tree_method=_normalize_tree_method(args.tree_method),
        iqtree_fast=iqtree_fast,
        iqtree_model=args.iqtree_model,
        hmmsearch_cutoff=args.hmmsearch_cutoff,
        hmmsearch_evalue=args.hmmsearch_evalue,
        selection_mode=args.selection_mode,
        selection_max_rounds=args.selection_max_rounds,
        selection_global_rounds=args.selection_global_rounds,
        lock_references=lock_references,
        max_sdup=args.max_sdup,
        max_dupl=args.max_dupl,
        ref=_normalize_cli_path(args.ref) if args.ref else None,
        ref_concat=ref_concat,
        marker_selection=marker_selection_enabled,
        singles=singles,
        singles_mode=args.singles_mode,
        num_nei=args.num_nei,
        singles_min_rfdist=args.singles_min_rfdist,
        keep_intermediates=keep_intermediates,
        is_ref=is_ref,
        start_time=start_time,
        ani_cluster=ani_cluster,
        snp=snp,
        ani_threshold=args.ani_threshold,
        ani_backend=args.ani_backend,
        ani_mcl_inflation=args.ani_mcl_inflation,
        snp_tree_min_cluster_size=args.snp_tree_min_cluster_size,
    )


def parse_args() -> Config:
    """Parse command-line arguments and return a validated configuration."""
    parser = _build_parser()
    return _config_from_args(parser, parser.parse_args())


def _tree_label(cfg: Config) -> str:
    return "VeryFastTree" if cfg.tree_method == "veryfasttree" else cfg.tree_method


def _print_configuration(cfg: Config) -> None:
    print(
        f"arguments:\n"
        f" proteomes, fasta {cfg.genomedir}\n"
        f" models, hmm {cfg.modeldir}\n"
        f" working directory {cfg.outdir}\n"
        f" number of CPUs {cfg.num_cpus}\n"
        f" tree method {cfg.tree_method}\n"
        f" iqtree fast {'yes' if cfg.iqtree_fast else 'no'}\n"
        f" iqtree model {cfg.iqtree_model}\n"
        f" hmmsearch cutoff {cfg.hmmsearch_cutoff}\n"
        f" hmmsearch evalue {cfg.hmmsearch_evalue}\n"
        f" selection mode {cfg.selection_mode}\n"
        f" selection max rounds {cfg.selection_max_rounds}\n"
        f" selection global rounds {cfg.selection_global_rounds}\n"
        f" lock references {'yes' if cfg.lock_references else 'no'}\n"
        f" minimum percentage of models {cfg.percent_models}\n"
        f" max single-marker copies {cfg.max_sdup}\n"
        f" max duplicated-marker fraction {cfg.max_dupl}\n"
        f" singleton mode {cfg.singles_mode}\n"
        f" singleton neighborhood size {cfg.num_nei} (0=auto)\n"
        f" singleton min RF distance {cfg.singles_min_rfdist}\n"
        f" ANI clustering {'yes' if cfg.ani_cluster else 'no'}\n"
        f" SNP trees {'yes' if cfg.snp else 'no'}\n"
        f" ANI threshold {cfg.ani_threshold}\n"
        f" ANI backend {cfg.ani_backend}\n"
        f" ANI MCL inflation {cfg.ani_mcl_inflation}\n"
        f" SNP tree min cluster size {cfg.snp_tree_min_cluster_size}\n"
        f" reference directory {cfg.ref}\n"
        f" keep intermediates {'yes' if cfg.keep_intermediates else 'no'}\n"
        f" --marker_selection {'yes' if cfg.marker_selection else 'no'}\n"
    )
    if cfg.ref:
        print(f"--ref_concat {cfg.ref_dir_path()}\n")
    else:
        print("--ref_concat no reference directory\n")
    print("=" * 80)


def _search_and_extract(
    cfg: Config,
    timings: dict[str, tuple[datetime.datetime, float]],
) -> pd.DataFrame:
    cfg.model_count = search.concat_inputs(cfg)

    started_at = _local_now()
    search_time = search.run_hmmsearch(cfg)
    timings["running hmmsearch"] = (started_at, search_time)

    started_at = _local_now()
    timer_start = time.time()
    finaldf, _ = search.parse_hmmsearch(cfg)
    df, df_fordups = search.build_working_df(cfg, finaldf)
    extract.extract_hits(cfg, df)
    extract.write_extracted_sequences(cfg)
    extract_time = time.time() - timer_start
    print(f"\nextraction of best hits done - runtime: {extract_time:.1f} seconds")
    print("=" * 80)
    timings["extracting best hits"] = (started_at, extract_time)
    return df_fordups


def _align_and_remove_duplicates(
    cfg: Config,
    df_fordups: pd.DataFrame,
    timings: dict[str, tuple[datetime.datetime, float]],
) -> None:
    started_at = _local_now()
    timer_start = time.time()
    align.run_alignment(cfg)
    alignment_time = time.time() - timer_start
    print(f"\nalignment done - runtime: {alignment_time:.1f} seconds")
    print("=" * 80 + "\n")
    timings["running alignment"] = (started_at, alignment_time)

    timer_start = time.time()
    print("- ... eliminating duplicates\n")
    duplicates.eliminate_duplicates(cfg, df_fordups)
    duplicate_time = time.time() - timer_start
    print(
        "aligned to aln_SpecTree runtime, duplicates eliminated - "
        f"runtime: {duplicate_time:.1f} seconds"
    )
    print("=" * 80 + "\n")


def _trim_build_matrix_and_tree(
    cfg: Config,
    timings: dict[str, tuple[datetime.datetime, float]],
) -> None:
    started_at = _local_now()
    timer_start = time.time()
    print("- ...running trimal")
    supermatrix.run_trimal(cfg, cfg.aln_spectree_dir, cfg.trimmed_dir)
    trim_time = time.time() - timer_start
    print(f"\ntrimming done - runtime: {trim_time:.1f} seconds")
    print("=" * 80 + "\n")
    timings["running trimal"] = (started_at, trim_time)

    started_at = _local_now()
    timer_start = time.time()
    print("- ...creating supermatrix")
    table_path = str(Path(cfg.tables_dir) / "table_df_concatenated_w_X")
    concat_path = str(Path(cfg.concat_dir) / "concatenated.faa")
    supermatrix.build_supermatrix(
        cfg.trimmed_dir,
        cfg.concat_dir,
        table_path,
        concat_path,
    )
    concat_time = time.time() - timer_start
    print(f"\nsupermatrix created - runtime: {concat_time:.1f} seconds")
    print("=" * 80 + "\n")
    timings["creating supermatrix"] = (started_at, concat_time)

    started_at = _local_now()
    timer_start = time.time()
    tree_label = _tree_label(cfg)
    print(f"- ...running {tree_label}")
    tree_path = str(Path(cfg.outdir) / "tree.nwk")
    phylogeny.run_species_tree(cfg, concat_path, tree_path)
    tree_time = time.time() - timer_start
    print(f"\n{tree_label} done - total runtime: {tree_time:.1f} seconds")
    print("=" * 80 + "\n")
    timings[f"running {cfg.tree_method}"] = (started_at, tree_time)


def _run_basic_pipeline(cfg: Config) -> None:
    timings: dict[str, tuple[datetime.datetime, float]] = {}
    try:
        df_fordups = _search_and_extract(cfg, timings)
        _align_and_remove_duplicates(cfg, df_fordups, timings)
        _trim_build_matrix_and_tree(cfg, timings)
        sgtree_logging.write_logfile(cfg, timings)
    except Exception as exc:
        logger.exception("sgtree basic pipeline failed")
        print(f"ERROR: {exc.__doc__}\n {exc}")
        traceback.print_exc()
        raise


def _prepare_marker_trees(
    cfg: Config,
    timings: dict[str, tuple[datetime.datetime, float]],
) -> datetime.datetime:
    started_at = _local_now()
    timer_start = time.time()
    print("- ...running trimal (on cleaned proteomes)")
    trimmed_prot_dir = str(Path(cfg.outdir) / "trimmed_protTrees")
    supermatrix.run_trimal_simple(cfg, cfg.aligned_dir, trimmed_prot_dir)
    trim_time = time.time() - timer_start
    print(f"\ntrimming done - runtime: {trim_time:.1f} seconds")
    print("=" * 80 + "\n")
    timings["running trimal"] = (started_at, trim_time)

    marker_selection_started_at = _local_now()
    timer_start = time.time()
    tree_label = _tree_label(cfg)
    print(f"- ...running {tree_label}, making protein trees for marker selection:")
    treeout_dir = str(Path(cfg.outdir) / "treeouts_protTrees")
    phylogeny.run_fasttree_per_marker(cfg, trimmed_prot_dir, treeout_dir)
    tree_time = time.time() - timer_start
    print(f"\n{tree_label} done - total runtime: {tree_time:.1f} seconds")
    print("=" * 80 + "\n")
    return marker_selection_started_at


def _run_selection_rounds(cfg: Config, ls_refs: list[str] | None) -> float:
    timer_start = time.time()
    print("- ...starting marker selection (Noperm):\n")
    current_species_tree = str(Path(cfg.outdir) / "tree.nwk")
    previous_kept: dict[tuple[str, str], str] | None = None
    total_rounds = max(1, cfg.selection_global_rounds)

    for round_index in range(1, total_rounds + 1):
        kept = marker_selection.run_noperm(
            cfg,
            ls_refs,
            species_tree_path=current_species_tree,
            initial_kept=previous_kept,
        )
        converged = previous_kept is not None and kept == previous_kept
        use_singles = cfg.singles and (round_index == total_rounds or converged)
        rf_source = str(Path(cfg.outdir) / "marker_selection_rf_values.txt")
        if total_rounds > 1:
            rf_round = str(
                Path(cfg.outdir) / f"marker_selection_rf_values_round{round_index}.txt"
            )
            shutil.copyfile(rf_source, rf_round)

        if use_singles:
            marker_selection.remove_singles(
                cfg,
                species_tree_path=current_species_tree,
            )

        cleaned_seq_dir = marker_selection.write_cleaned_sequences(
            cfg,
            use_singles=use_singles,
        )
        print(f"- ...running {cfg.aln_method} for final alignment:")
        aligned_final_dir = str(Path(cfg.outdir) / "aligned_final")
        align.run_alignment(
            cfg,
            extracted_seqs_dir=cleaned_seq_dir,
            aligned_dir=aligned_final_dir,
        )

        print("- ...running trimal for final alignment:")
        trimmed_final_dir = str(Path(cfg.outdir) / "trimmed_final")
        supermatrix.run_trimal_simple(cfg, aligned_final_dir, trimmed_final_dir)
        print("\ntrimming done")
        print("=" * 80 + "\n")

        print("- ...creating supermatrix")
        table_path = str(Path(cfg.tables_dir) / "table_df_concatenated_w_X_final")
        concat_final_dir = str(Path(cfg.outdir) / "concat_final")
        concat_final_path = str(Path(concat_final_dir) / "concatenated.faa")
        supermatrix.build_supermatrix(
            trimmed_final_dir,
            concat_final_dir,
            table_path,
            concat_final_path,
        )
        print("\nsupermatrix created")
        print("=" * 80 + "\n")

        print(f"- ...running {cfg.tree_method} for tree_final.nwk")
        tree_final_path = str(Path(cfg.outdir) / "tree_final.nwk")
        phylogeny.run_species_tree(cfg, concat_final_path, tree_final_path)
        if total_rounds > 1:
            round_tree_path = str(Path(cfg.outdir) / f"tree_round_{round_index}.nwk")
            shutil.copyfile(tree_final_path, round_tree_path)

        if converged:
            print(f"- ...marker selection converged after round {round_index}")
            break
        previous_kept = kept
        current_species_tree = tree_final_path

    return time.time() - timer_start


def _run_marker_selection(cfg: Config, ls_refs: list[str] | None) -> None:
    if not cfg.marker_selection:
        return

    timings: dict[str, tuple[datetime.datetime, float]] = {}
    try:
        started_at = _prepare_marker_trees(cfg, timings)
        marker_selection_time = _run_selection_rounds(cfg, ls_refs)
        print(f"Marker selection runtime {marker_selection_time:.1f}")
        print("=" * 80 + "\n")
        timings["Marker selection"] = (started_at, marker_selection_time)
        sgtree_logging.append_logfile(cfg, timings)
    except Exception as exc:
        logger.exception("sgtree marker-selection phase failed")
        print(f"ERROR in marker selection: {exc.__doc__}\n {exc}")
        traceback.print_exc()
        raise


def _finalize_outputs(cfg: Config) -> None:
    if cfg.ani_cluster and cfg.snp and not cfg.is_ref:
        ani_clustering.build_cluster_snp_trees(cfg)

    if not cfg.is_ref and not cfg.keep_intermediates:
        print("- ...finalizing outputs (archiving intermediates)")
        if cfg.marker_selection:
            cleanup.cleanup_marker_selection(cfg.outdir)
        else:
            cleanup.cleanup_basic(cfg.outdir)
        print("- ...finalization complete")

    tree_name = "tree_final.nwk" if cfg.marker_selection else "tree.nwk"
    final_tree_path = Path(cfg.outdir) / tree_name
    if final_tree_path.exists():
        print(f"Final tree: {final_tree_path.resolve()}")


def main() -> None:
    """Run the configured SGTree workflow."""
    cfg = parse_args()
    Path(cfg.outdir).mkdir(parents=True, exist_ok=True)
    cfg.print_banner()
    print(
        f"-... Reference directory and arguments\n"
        f" Reference directory located at {cfg.ref_concat}\n ...starting sgtree..."
    )

    _clean_previous_run(cfg)
    _mark_run_directory(cfg)
    if cfg.ani_cluster and not cfg.is_ref:
        ani_clustering.prepare_ani_cluster_inputs(cfg)

    ls_refs = reference.prepare_reference(cfg)
    _print_configuration(cfg)
    _run_basic_pipeline(cfg)
    _run_marker_selection(cfg, ls_refs)
    _finalize_outputs(cfg)
    print(f"START: {cfg.start_time} END {_local_now()}")


if __name__ == "__main__":
    main()
