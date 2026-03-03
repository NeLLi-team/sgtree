import os
import sys
import glob
import subprocess
import shutil
import zipfile

from sgtree.config import Config


def check_duplicate_proteomes(genomedir: str, refdir: str):
    """Check for duplicate proteomes between query and reference directories."""
    if os.path.isdir(genomedir):
        ls_genomes = [
            os.path.basename(f).split(".")[0]
            for f in glob.glob(os.path.join(genomedir, "*"))
        ]
    else:
        ls_genomes = []
        with open(genomedir) as fi:
            for ln in fi:
                if ln.startswith(">"):
                    gid = ln[1:].split("|")[0]
                    if gid not in ls_genomes:
                        ls_genomes.append(gid)

    if os.path.isdir(refdir):
        ls_ref = [
            os.path.basename(f).split(".")[0]
            for f in glob.glob(os.path.join(refdir, "*"))
        ]
    else:
        ls_ref = []
        with open(refdir) as fi:
            for ln in fi:
                if ln.startswith(">"):
                    gid = ln[1:].split("|")[0]
                    if gid not in ls_ref:
                        ls_ref.append(gid)

    common = [item for item in ls_genomes if item in ls_ref]
    if common:
        print(
            "WARNING: Duplicate proteomes found in reference directory, "
            "this will cause errors for --marker_selection yes (Noperm)\n"
            "please delete the following from either directory/concat file and try again."
        )
        for each in common:
            print("please delete", each)
        sys.exit(1)

    return ls_ref


def prepare_reference(cfg: Config) -> list[str] | None:
    """Prepare reference directory, running sgtree on references if needed.

    Returns list of reference filenames (e.g. ['IMG123.faa']) or None.
    """
    if cfg.ref is None:
        print("no reference directory")
        return None

    ref = cfg.ref.rstrip("/")
    ls_ref = check_duplicate_proteomes(cfg.genomedir, ref)
    ls_refs = [x + ".faa" for x in ls_ref]
    print(ls_refs)

    ref_dir = cfg.ref_dir_path()
    hmmout_path = os.path.join(ref_dir, "hits.hmmout")
    table_path = os.path.join(ref_dir, "table_elim_dups")

    if os.path.exists(hmmout_path):
        if os.path.exists(table_path):
            print(f"already have files for references at {ref_dir}/")
        else:
            print(
                f"WARNING: Malformed reference directory found please delete directory\n"
                f" {ref_dir}\n"
                " and try again with properly formatted input data (or more genomes/models)."
            )
            sys.exit(1)
    else:
        cmd = [
            sys.executable, "-m", "sgtree", ref, cfg.modeldir,
            "--num_cpus", str(cfg.num_cpus),
            "--aln", cfg.aln_method,
            "--tree_method", cfg.tree_method,
            "--iqtree_fast", "yes" if cfg.iqtree_fast else "no",
            "--iqtree_model", cfg.iqtree_model,
            "--hmmsearch_cutoff", cfg.hmmsearch_cutoff,
            "--hmmsearch_evalue", str(cfg.hmmsearch_evalue),
            "--max_sdup", str(cfg.max_sdup),
            "--max_dupl", str(cfg.max_dupl),
            "--save_dir", ref_dir,
            "--is_ref", "yes",
        ]
        print("- ... Creating new reference directory\n", cmd)
        subprocess.run(cmd, stdout=subprocess.PIPE, check=True)

        # archive intermediate files
        for filepath in glob.glob(os.path.join(ref_dir, "*")):
            basename = os.path.basename(filepath)
            if os.path.isdir(filepath):
                if basename in ("tables", "concat", "extracted_seqs"):
                    continue
                subprocess.run(
                    ["zip", filepath] + glob.glob(os.path.join(filepath, "*")),
                    stdout=subprocess.DEVNULL,
                )
                shutil.rmtree(filepath)
            else:
                keep_files = (
                    "marker_count_matrix.csv", "proteomes",
                    "hits.hmmout", "table_elim_dups",
                )
                if basename in keep_files:
                    continue
                with zipfile.ZipFile(filepath, "w") as myzip:
                    myzip.write(filepath)

        os.makedirs(os.path.join(ref_dir, "temp"), exist_ok=True)

        for f in glob.glob(os.path.join(ref_dir, "*.zip")):
            shutil.move(f, os.path.join(ref_dir, "temp"))
        for f in glob.glob(os.path.join(ref_dir, "*.txt")):
            shutil.move(f, os.path.join(ref_dir, "temp"))

        for name in ("ref_and_query_proteomes", "models", "tree.nwk"):
            src = os.path.join(ref_dir, name)
            if os.path.exists(src):
                shutil.move(src, os.path.join(ref_dir, "temp"))

        for f in glob.glob(os.path.join(ref_dir, "temp", "*")):
            with zipfile.ZipFile(f, "w") as myzip:
                myzip.write(f)

    return ls_refs
