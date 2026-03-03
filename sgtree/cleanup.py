import os
import glob
import shutil
import subprocess
import zipfile


def cleanup_basic(outdir: str):
    """Archive intermediate files for basic run (no marker selection)."""
    keep_extensions = {"txt", "png", "csv", "nwk"}
    keep_names = {"aligned_final", "concat", "tree.nwk", "hits.hmmout",
                  "ref_and_query_proteomes", "marker_count_matrix.csv"}

    for filepath in glob.glob(os.path.join(outdir, "*")):
        basename = os.path.basename(filepath)

        # skip files we want to keep
        ext = basename.split(".")[-1]
        if ext in keep_extensions:
            continue
        if basename.startswith("logfile"):
            continue
        if basename in keep_names:
            continue

        if os.path.isdir(filepath):
            subprocess.run(
                ["zip", filepath] + glob.glob(os.path.join(filepath, "*")),
                stdout=subprocess.DEVNULL,
            )
            shutil.rmtree(filepath)
        else:
            with zipfile.ZipFile(filepath, "w") as myzip:
                myzip.write(filepath)

    # organize into temp directories
    os.makedirs(os.path.join(outdir, "temp"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "temp", "itol"), exist_ok=True)

    for f in glob.glob(os.path.join(outdir, "*.zip")):
        shutil.move(f, os.path.join(outdir, "temp"))

    for name in ("models", "proteomes", "table_elim_dups", "hits.hmmout", "ref_and_query_proteomes"):
        src = os.path.join(outdir, name)
        if os.path.exists(src):
            shutil.move(src, os.path.join(outdir, "temp"))

    for name in ("color.txt", "marker_count.txt"):
        src = os.path.join(outdir, name)
        if os.path.exists(src):
            shutil.move(src, os.path.join(outdir, "temp", "itol"))


def cleanup_marker_selection(outdir: str):
    """Archive intermediate files for marker selection run."""
    keep_names = {
        "tree_final.nwk", "hits.hmmout", "marker_count_matrix.csv",
        "ref_and_query_proteomes", "concat_final",
        "marker_selection_rf_values.txt",
    }

    for filepath in glob.glob(os.path.join(outdir, "*")):
        basename = os.path.basename(filepath)

        ext = basename.split(".")[-1]
        if ext in ("txt", "png"):
            continue
        if basename.startswith("logfile"):
            continue
        if basename in keep_names:
            continue

        if os.path.isdir(filepath):
            subprocess.run(
                ["zip", filepath] + glob.glob(os.path.join(filepath, "*")),
                stdout=subprocess.DEVNULL,
            )
            shutil.rmtree(filepath)
        else:
            with zipfile.ZipFile(filepath, "w") as myzip:
                myzip.write(filepath)

    os.makedirs(os.path.join(outdir, "temp"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "temp", "itol"), exist_ok=True)

    for f in glob.glob(os.path.join(outdir, "*.zip")):
        shutil.move(f, os.path.join(outdir, "temp"))

    for name in ("models", "proteomes", "table_elim_dups", "tree.nwk",
                  "hits.hmmout", "ref_and_query_proteomes"):
        src = os.path.join(outdir, name)
        if os.path.exists(src):
            shutil.move(src, os.path.join(outdir, "temp"))

    for name in ("color.txt", "marker_counts.txt"):
        src = os.path.join(outdir, name)
        if os.path.exists(src):
            shutil.move(src, os.path.join(outdir, "temp", "itol"))
