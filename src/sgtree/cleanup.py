"""Archive SGTree intermediates without touching unrelated output files."""

from __future__ import annotations

import logging
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

logger = logging.getLogger(__name__)

RUN_MARKER_NAME = ".sgtree-run"
RUN_MARKER_CONTENT = "SGTree managed output directory\n"

# Top-level entries under the output directory that SGTree itself generates.
GENERATED_ENTRY_NAMES = frozenset(
    {
        # directories
        "models",
        "models_split",
        "proteomes",
        "staged_proteomes",
        "tables",
        "extracted",
        "extracted_seqs",
        "aligned",
        "aln_SpecTree",
        "trimmed_SpeciesTree",
        "concat",
        "ref_and_query_proteomes",
        "ani",
        "snp_trees",
        "temp",
        "protTrees",
        "removed",
        "trimmed_protTrees",
        "treeouts_protTrees",
        "aligned_final",
        "trimmed_final",
        "concat_final",
        "extracted_final",
        # files
        "tree.nwk",
        "tree_final.nwk",
        "hits.hmmout",
        "hits.hmmout.del.ls",
        "hits.hmmout.lfilt",
        "table_elim_dups",
        "marker_count_matrix.csv",
        "marker_selection_rf_values.txt",
        "proteomes_header_map.tsv",
        "genome_manifest.tsv",
        "gene_calls.tsv",
        "log_genomes_removed.txt",
        "singleton_candidates.tsv",
        "cache_meta.json",
    }
)

_LOGFILE_PATTERN = re.compile(
    r"logfile_\d{4}_\d{2}_\d{2}_(?:\d{2}:\d{2}:\d{2}|\d{2}_\d{2}_\d{2})\.txt"
)
_ROUND_TREE_PATTERN = re.compile(r"tree_round_\d+\.nwk")
_ROUND_RF_PATTERN = re.compile(r"marker_selection_rf_values_round\d+\.txt")
_IQTREE_PREFIXES = ("tree.nwk.iqtree.", "tree_final.nwk.iqtree.")
_IQTREE_SUFFIXES = frozenset(
    {
        "bionj",
        "boottrees",
        "ckp.gz",
        "contree",
        "iqtree",
        "log",
        "mldist",
        "model.gz",
        "splits.nex",
        "treefile",
        "ufboot",
        "uniqueseq.phy",
    }
)


def is_generated_entry(basename: str) -> bool:
    """Return whether a top-level name is reserved for an SGTree output."""
    if basename in GENERATED_ENTRY_NAMES:
        return True
    if _LOGFILE_PATTERN.fullmatch(basename):
        return True
    if _ROUND_TREE_PATTERN.fullmatch(basename):
        return True
    if _ROUND_RF_PATTERN.fullmatch(basename):
        return True
    return any(
        basename.removeprefix(prefix) in _IQTREE_SUFFIXES
        for prefix in _IQTREE_PREFIXES
        if basename.startswith(prefix)
    )


def _path_contains_symlink(path: Path) -> bool:
    """Return whether path is a link or contains one below a directory."""
    if path.is_symlink():
        return True
    if not path.is_dir():
        return False
    return any(entry.is_symlink() for entry in path.rglob("*"))


def _write_archive(source: Path, destination: Path) -> None:
    with tempfile.TemporaryDirectory(
        dir=source.parent,
        prefix=f".{source.name}.",
    ) as temp_dir:
        temporary_archive = Path(temp_dir) / destination.name
        with zipfile.ZipFile(temporary_archive, "w", zipfile.ZIP_DEFLATED) as zip_file:
            zip_file.write(source, source.name)
            if source.is_dir():
                for entry in source.rglob("*"):
                    zip_file.write(entry, entry.relative_to(source.parent))
        temporary_archive.replace(destination)


def _archive_directory(dirpath: str) -> Path | None:
    """Archive a directory recursively and remove it after a complete ZIP exists."""
    source = Path(dirpath)
    if _path_contains_symlink(source):
        logger.warning("leaving linked path unarchived: %s", source)
        return None

    archive = Path(f"{source}.zip")
    if archive.exists() or archive.is_symlink():
        raise FileExistsError(f"refusing to overwrite {archive}")

    _write_archive(source, archive)
    shutil.rmtree(source)
    return archive


def _zip_file_in_place(filepath: str) -> None:
    """Atomically replace one file with a ZIP archive under the same name."""
    source = Path(filepath)
    if source.is_symlink():
        logger.warning("leaving linked path unarchived: %s", source)
        return

    _write_archive(source, source)


def _prepare_temp_directory(outdir: Path) -> Path:
    temp_dir = outdir / "temp"
    if temp_dir.is_symlink():
        raise FileExistsError(f"refusing to use linked temp directory: {temp_dir}")
    temp_dir.mkdir(exist_ok=True)
    if not temp_dir.is_dir():
        raise NotADirectoryError(f"cleanup temp path is not a directory: {temp_dir}")
    return temp_dir


def _move_to_directory(source: Path, destination_dir: Path) -> None:
    if _path_contains_symlink(source):
        logger.warning("leaving linked path in place: %s", source)
        return
    destination = destination_dir / source.name
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"refusing to overwrite {destination}")
    source.rename(destination)


def _archive_generated_entries(
    outdir: Path,
    temp_dir: Path,
    *,
    keep_extensions: set[str],
    keep_names: set[str],
) -> list[Path]:
    paths = [
        path
        for path in outdir.iterdir()
        if is_generated_entry(path.name)
        and path.suffix.removeprefix(".") not in keep_extensions
        and path.name not in keep_names
    ]
    for path in paths:
        if not path.is_dir() or _path_contains_symlink(path):
            continue
        for destination in (Path(f"{path}.zip"), temp_dir / f"{path.name}.zip"):
            if destination.exists() or destination.is_symlink():
                raise FileExistsError(f"refusing to overwrite {destination}")

    archives: list[Path] = []
    for path in paths:
        if path.is_dir():
            archive = _archive_directory(str(path))
            if archive is not None:
                archives.append(archive)
        else:
            _zip_file_in_place(str(path))
    return archives


def _check_temp_destinations(
    outdir: Path,
    temp_dir: Path,
    names: tuple[str, ...],
) -> None:
    for name in names:
        source = outdir / name
        destination = temp_dir / name
        if (source.exists() or source.is_symlink()) and (
            destination.exists() or destination.is_symlink()
        ):
            raise FileExistsError(f"refusing to overwrite {destination}")


def _finalize_temp_entries(
    outdir: Path,
    temp_dir: Path,
    archives: list[Path],
    names: tuple[str, ...],
) -> None:
    for archive in archives:
        _move_to_directory(archive, temp_dir)
    for name in names:
        source = outdir / name
        if source.exists() or source.is_symlink():
            _move_to_directory(source, temp_dir)


def cleanup_basic(outdir: str) -> None:
    """Archive intermediate files after a basic run."""
    output_dir = Path(outdir)
    temp_dir = _prepare_temp_directory(output_dir)
    moved_names = (
        "models",
        "proteomes",
        "table_elim_dups",
        "hits.hmmout",
        "ref_and_query_proteomes",
    )
    _check_temp_destinations(output_dir, temp_dir, moved_names)
    archives = _archive_generated_entries(
        output_dir,
        temp_dir,
        keep_extensions={"txt", "png", "csv", "nwk"},
        keep_names={
            "aligned_final",
            "concat",
            "tree.nwk",
            "hits.hmmout",
            "ref_and_query_proteomes",
            "marker_count_matrix.csv",
            "ani",
            "snp_trees",
            "genome_manifest.tsv",
            "proteomes_header_map.tsv",
            "gene_calls.tsv",
            "singleton_candidates.tsv",
            "temp",
        },
    )
    _finalize_temp_entries(
        output_dir,
        temp_dir,
        archives,
        moved_names,
    )


def cleanup_marker_selection(outdir: str) -> None:
    """Archive intermediate files after a marker-selection run."""
    output_dir = Path(outdir)
    temp_dir = _prepare_temp_directory(output_dir)
    moved_names = (
        "models",
        "proteomes",
        "table_elim_dups",
        "tree.nwk",
        "hits.hmmout",
        "ref_and_query_proteomes",
    )
    _check_temp_destinations(output_dir, temp_dir, moved_names)
    archives = _archive_generated_entries(
        output_dir,
        temp_dir,
        keep_extensions={"txt", "png", "nwk"},
        keep_names={
            "tree_final.nwk",
            "hits.hmmout",
            "marker_count_matrix.csv",
            "ref_and_query_proteomes",
            "concat_final",
            "marker_selection_rf_values.txt",
            "ani",
            "snp_trees",
            "genome_manifest.tsv",
            "proteomes_header_map.tsv",
            "gene_calls.tsv",
            "singleton_candidates.tsv",
            "temp",
        },
    )
    tree_path = output_dir / "tree.nwk"
    if tree_path.exists() or tree_path.is_symlink():
        _zip_file_in_place(str(tree_path))
    _finalize_temp_entries(
        output_dir,
        temp_dir,
        archives,
        moved_names,
    )
