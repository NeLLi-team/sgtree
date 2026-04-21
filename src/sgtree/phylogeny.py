import os
import subprocess
import glob
import shutil

from sgtree._subprocess import run_capture, run_check
from sgtree.config import Config
from sgtree.parallel import map_threaded


def _fasttree_executable() -> str:
    for candidate in ("VeryFastTree", "veryfasttree"):
        if shutil.which(candidate):
            return candidate
    raise FileNotFoundError("Could not find VeryFastTree in PATH")


# Per-executable cache of whether the binary accepts the `-threads` flag.
# Populated lazily from the first real run; subsequent calls skip the probe.
_THREADS_SUPPORT: dict[str, bool] = {}


def _run_fasttree_with_optional_threads(threaded_cmd: list[str], fallback_cmd: list[str]) -> None:
    executable = threaded_cmd[0]
    if _THREADS_SUPPORT.get(executable) is False:
        run_check(fallback_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return

    threaded = run_capture(threaded_cmd)
    if threaded.returncode == 0:
        _THREADS_SUPPORT[executable] = True
        return

    message = (threaded.stderr or "") + "\n" + (threaded.stdout or "")
    if "-threads" not in message:
        raise subprocess.CalledProcessError(
            returncode=threaded.returncode,
            cmd=threaded_cmd,
            output=threaded.stdout,
            stderr=threaded.stderr,
        )
    _THREADS_SUPPORT[executable] = False
    run_check(fallback_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_fasttree(input_fasta: str, output_tree: str, threads: int):
    """Run VeryFastTree to build an approximate ML tree."""
    executable = _fasttree_executable()
    threaded_cmd = [executable, "-threads", str(max(1, threads)), "-quiet", "-out", output_tree, input_fasta]
    fallback_cmd = [executable, "-quiet", "-out", output_tree, input_fasta]
    _run_fasttree_with_optional_threads(threaded_cmd, fallback_cmd)


def _run_iqtree(input_fasta: str, output_tree: str, cpus: int, model: str, fast: bool):
    """Run IQ-TREE and copy resulting treefile to output_tree."""
    prefix = output_tree + ".iqtree"
    cmd = [
        "iqtree",
        "--quiet",
        "--prefix", prefix,
        "-m", model,
        "-T", str(max(1, cpus)),
    ]
    if fast:
        cmd.append("-fast")
    cmd.extend(["-s", input_fasta])
    run_check(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    treefile = prefix + ".treefile"
    if not os.path.exists(treefile):
        raise FileNotFoundError(f"IQ-TREE did not produce expected treefile: {treefile}")
    shutil.copyfile(treefile, output_tree)


def run_species_tree(cfg: Config, input_fasta: str, output_tree: str):
    """Run selected tree method for the species tree."""
    if cfg.tree_method == "iqtree":
        _run_iqtree(input_fasta, output_tree, cfg.num_cpus, cfg.iqtree_model, cfg.iqtree_fast)
    else:
        run_fasttree(input_fasta, output_tree, cfg.num_cpus)


def run_snp_tree(cfg: Config, input_fasta: str, output_tree: str):
    """Run a nucleotide tree for cluster-level SNP alignments."""
    if cfg.tree_method == "iqtree":
        prefix = output_tree + ".iqtree"
        cmd = [
            "iqtree",
            "--quiet",
            "--prefix", prefix,
            "-m", "GTR+G",
            "-st", "DNA",
            "-T", "1",
        ]
        if cfg.iqtree_fast:
            cmd.append("-fast")
        cmd.extend(["-s", input_fasta])
        run_check(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        treefile = prefix + ".treefile"
        if not os.path.exists(treefile):
            raise FileNotFoundError(f"IQ-TREE did not produce expected treefile: {treefile}")
        shutil.copyfile(treefile, output_tree)
    else:
        executable = _fasttree_executable()
        threaded_cmd = [executable, "-threads", "1", "-nt", "-gtr", "-quiet", "-out", output_tree, input_fasta]
        fallback_cmd = [executable, "-nt", "-gtr", "-quiet", "-out", output_tree, input_fasta]
        _run_fasttree_with_optional_threads(threaded_cmd, fallback_cmd)


def _build_tree_worker(args):
    """Worker: build a single protein tree with selected tree method."""
    filepath, treeout_dir, tree_method, iqtree_model, iqtree_fast = args
    tree_out = os.path.join(
        treeout_dir,
        os.path.basename(filepath) + "_tree.out",
    )
    if tree_method == "iqtree":
        _run_iqtree(filepath, tree_out, 1, iqtree_model, iqtree_fast)
    else:
        executable = _fasttree_executable()
        threaded_cmd = [executable, "-threads", "1", "-quiet", "-out", tree_out, filepath]
        fallback_cmd = [executable, "-quiet", "-out", tree_out, filepath]
        _run_fasttree_with_optional_threads(threaded_cmd, fallback_cmd)

def run_fasttree_per_marker(cfg: Config, trimmed_dir: str, treeout_dir: str):
    """Build individual protein trees for each trimmed marker alignment."""
    os.makedirs(treeout_dir, exist_ok=True)

    files = glob.glob(os.path.join(trimmed_dir, "*"))
    args = [
        (f, treeout_dir, cfg.tree_method, cfg.iqtree_model, cfg.iqtree_fast)
        for f in files
    ]
    if not args:
        return

    workers = max(1, min(cfg.num_cpus, len(args)))
    map_threaded(_build_tree_worker, args, workers)
