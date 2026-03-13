import os
import sys
import glob
import csv
import json
import subprocess
import shutil
import zipfile

from sgtree.config import Config


REF_CACHE_META_FILE = "cache_meta.json"
REF_CACHE_SCHEMA_VERSION = 1
REQUIRED_TABLE_COLUMNS = ("savedname", "namemodel", "score_bits")


def _cache_signature(cfg: Config) -> dict:
    return {
        "schema_version": REF_CACHE_SCHEMA_VERSION,
        "modeldir": os.path.abspath(cfg.modeldir),
        "percent_models": cfg.percent_models,
        "max_sdup": cfg.max_sdup,
        "max_dupl": cfg.max_dupl,
        "hmmsearch_cutoff": cfg.hmmsearch_cutoff,
        "hmmsearch_evalue": cfg.hmmsearch_evalue,
        "aln_method": "hmmalign",
        "ref_lflt": 0,
        "tree_method": cfg.tree_method,
        "iqtree_fast": cfg.iqtree_fast,
        "iqtree_model": cfg.iqtree_model,
    }


def _cache_meta_path(ref_dir: str) -> str:
    return os.path.join(ref_dir, REF_CACHE_META_FILE)


def _write_cache_meta(cfg: Config, ref_dir: str):
    with open(_cache_meta_path(ref_dir), "w") as f:
        json.dump(_cache_signature(cfg), f, sort_keys=True, indent=2)
        f.write("\n")


def _load_cache_meta(ref_dir: str) -> dict | None:
    meta_path = _cache_meta_path(ref_dir)
    if not os.path.exists(meta_path):
        return None
    with open(meta_path) as f:
        return json.load(f)


def _validate_table_schema(table_path: str) -> tuple[bool, str]:
    try:
        with open(table_path, newline="") as f:
            header = next(csv.reader(f), [])
    except OSError as e:
        return False, f"failed to read table schema ({e})"

    missing = [col for col in REQUIRED_TABLE_COLUMNS if col not in header]
    if missing:
        return False, f"table_elim_dups missing required columns: {', '.join(missing)}"
    return True, ""


def _validate_reference_cache(cfg: Config, ref_dir: str) -> tuple[bool, str]:
    hmmout_path = os.path.join(ref_dir, "hits.hmmout")
    table_path = os.path.join(ref_dir, "table_elim_dups")

    if not os.path.exists(hmmout_path):
        return False, "missing hits.hmmout"
    if not os.path.exists(table_path):
        return False, "missing table_elim_dups"

    schema_ok, schema_reason = _validate_table_schema(table_path)
    if not schema_ok:
        return False, schema_reason

    meta = _load_cache_meta(ref_dir)
    if meta is None:
        return False, f"missing {REF_CACHE_META_FILE}"

    expected = _cache_signature(cfg)
    for key, expected_val in expected.items():
        if meta.get(key) != expected_val:
            return False, f"cache metadata mismatch for '{key}'"

    return True, ""


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
    cache_valid, cache_reason = _validate_reference_cache(cfg, ref_dir)
    if cache_valid:
        print(f"already have files for references at {ref_dir}/")
        return ls_refs

    if os.path.exists(ref_dir):
        print(f"reference cache invalid ({cache_reason}); rebuilding {ref_dir}/")
        shutil.rmtree(ref_dir)

    cmd = [
        sys.executable, "-m", "sgtree", ref, cfg.modeldir,
        "--num_cpus", str(cfg.num_cpus),
        "--percent_models", str(cfg.percent_models),
        "--lflt", "0",
        "--aln", "hmmalign",
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
                "hits.hmmout", "table_elim_dups", "genome_manifest.tsv",
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

    _write_cache_meta(cfg, ref_dir)

    return ls_refs
