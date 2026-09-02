import os
import fnmatch
import glob
import shutil
import zipfile


# Top-level entries under the output directory that SGTree itself generates.
# Cleanup archives only these, and cli._clean_previous_run removes only these,
# so anything else a user keeps in the output directory is left alone.
GENERATED_ENTRY_NAMES = frozenset({
    # directories
    "models", "models_split", "proteomes", "staged_proteomes", "tables",
    "extracted", "extracted_seqs", "aligned", "aln_SpecTree",
    "trimmed_SpeciesTree", "concat", "ref_and_query_proteomes", "ani",
    "snp_trees", "temp", "protTrees", "removed", "trimmed_protTrees",
    "treeouts_protTrees", "aligned_final", "trimmed_final", "concat_final",
    "extracted_final",
    # files
    "tree.nwk", "tree_final.nwk", "hits.hmmout", "table_elim_dups",
    "marker_count_matrix.csv", "marker_selection_rf_values.txt",
    "proteomes_header_map.tsv", "genome_manifest.tsv", "gene_calls.tsv",
    "log_genomes_removed.txt", "singleton_candidates.tsv",
})

# Generated entries whose names carry a timestamp, a round number, or a
# tool-specific suffix.
GENERATED_ENTRY_GLOBS = (
    "logfile_[0-9]*.txt",
    "tree_round_[0-9]*.nwk",
    "marker_selection_rf_values_round[0-9]*.txt",
    "hits.hmmout.*",
    "tree.nwk.iqtree.*",
    "tree_final.nwk.iqtree.*",
)


def is_generated_entry(basename: str) -> bool:
    """True when SGTree wrote this top-level entry of the output directory."""
    if basename in GENERATED_ENTRY_NAMES:
        return True
    return any(fnmatch.fnmatch(basename, pattern) for pattern in GENERATED_ENTRY_GLOBS)


def _archive_directory(dirpath: str) -> None:
    """Archive a directory to <dirpath>.zip and remove the original.

    Walks the tree, so nested directories keep their contents. Uses zipfile
    rather than the zip binary, which is not part of the environment and which
    fails on an empty directory.
    """
    archive = dirpath + ".zip"
    if os.path.exists(archive):
        # never overwrite a file the user put there under the archive's name
        raise FileExistsError(f"refusing to overwrite {archive}")
    parent = os.path.dirname(dirpath.rstrip(os.sep)) or "."
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(dirpath):
            zf.write(root, os.path.relpath(root, parent))
            for name in files:
                path = os.path.join(root, name)
                zf.write(path, os.path.relpath(path, parent))
    shutil.rmtree(dirpath)


def _zip_file_in_place(filepath: str) -> None:
    """Compress a single file into a zip archive, replacing the original."""
    zip_path = filepath + ".zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(filepath, os.path.basename(filepath))
    os.remove(filepath)
    os.rename(zip_path, filepath)


def cleanup_basic(outdir: str):
    """Archive intermediate files for basic run (no marker selection)."""
    keep_extensions = {"txt", "png", "csv", "nwk"}
    keep_names = {"aligned_final", "concat", "tree.nwk", "hits.hmmout",
                  "ref_and_query_proteomes", "marker_count_matrix.csv",
                  "ani", "snp_trees", "genome_manifest.tsv",
                  "proteomes_header_map.tsv", "gene_calls.tsv",
                  "singleton_candidates.tsv", "temp"}

    for filepath in glob.glob(os.path.join(outdir, "*")):
        basename = os.path.basename(filepath)

        # never touch anything SGTree did not write into the output directory
        if not is_generated_entry(basename):
            continue

        # skip files we want to keep
        ext = basename.split(".")[-1]
        if ext in keep_extensions:
            continue
        if basename.startswith("logfile"):
            continue
        if basename in keep_names:
            continue

        if os.path.isdir(filepath):
            _archive_directory(filepath)
        else:
            _zip_file_in_place(filepath)

    # organize into temp directories
    os.makedirs(os.path.join(outdir, "temp"), exist_ok=True)
    for f in glob.glob(os.path.join(outdir, "*.zip")):
        # only archives cleanup just made, never a user's own zip file
        if not is_generated_entry(os.path.basename(f)[: -len(".zip")]):
            continue
        shutil.move(f, os.path.join(outdir, "temp"))

    for name in ("models", "proteomes", "table_elim_dups", "hits.hmmout", "ref_and_query_proteomes"):
        src = os.path.join(outdir, name)
        if os.path.exists(src):
            shutil.move(src, os.path.join(outdir, "temp"))


def cleanup_marker_selection(outdir: str):
    """Archive intermediate files for marker selection run."""
    keep_names = {
        "tree_final.nwk", "hits.hmmout", "marker_count_matrix.csv",
        "ref_and_query_proteomes", "concat_final",
        "marker_selection_rf_values.txt",
        "ani", "snp_trees", "genome_manifest.tsv",
        "proteomes_header_map.tsv", "gene_calls.tsv",
        "singleton_candidates.tsv", "temp",
    }

    for filepath in glob.glob(os.path.join(outdir, "*")):
        basename = os.path.basename(filepath)

        # never touch anything SGTree did not write into the output directory
        if not is_generated_entry(basename):
            continue

        ext = basename.split(".")[-1]
        if ext in ("txt", "png"):
            continue
        if basename.startswith("logfile"):
            continue
        if basename in keep_names:
            continue
        if basename.startswith("tree_round_") and basename.endswith(".nwk"):
            continue

        if os.path.isdir(filepath):
            _archive_directory(filepath)
        else:
            _zip_file_in_place(filepath)

    os.makedirs(os.path.join(outdir, "temp"), exist_ok=True)
    for f in glob.glob(os.path.join(outdir, "*.zip")):
        # only archives cleanup just made, never a user's own zip file
        if not is_generated_entry(os.path.basename(f)[: -len(".zip")]):
            continue
        shutil.move(f, os.path.join(outdir, "temp"))

    for name in ("models", "proteomes", "table_elim_dups", "tree.nwk",
                  "hits.hmmout", "ref_and_query_proteomes"):
        src = os.path.join(outdir, name)
        if os.path.exists(src):
            shutil.move(src, os.path.join(outdir, "temp"))
