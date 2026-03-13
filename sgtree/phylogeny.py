import os
import subprocess
import glob
import shutil

from sgtree.config import Config
from sgtree.parallel import map_threaded


def run_fasttree(input_fasta: str, output_tree: str):
    """Run FastTree to build an approximate ML tree."""
    cmd = ["FastTree", "-quiet", "-out", output_tree, input_fasta]
    subprocess.run(cmd, stdout=subprocess.PIPE, check=True)


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
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    treefile = prefix + ".treefile"
    if not os.path.exists(treefile):
        raise FileNotFoundError(f"IQ-TREE did not produce expected treefile: {treefile}")
    shutil.copyfile(treefile, output_tree)


def run_species_tree(cfg: Config, input_fasta: str, output_tree: str):
    """Run selected tree method for the species tree."""
    if cfg.tree_method == "iqtree":
        _run_iqtree(input_fasta, output_tree, cfg.num_cpus, cfg.iqtree_model, cfg.iqtree_fast)
    else:
        run_fasttree(input_fasta, output_tree)


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
            "-T", str(max(1, cfg.num_cpus)),
        ]
        if cfg.iqtree_fast:
            cmd.append("-fast")
        cmd.extend(["-s", input_fasta])
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        treefile = prefix + ".treefile"
        if not os.path.exists(treefile):
            raise FileNotFoundError(f"IQ-TREE did not produce expected treefile: {treefile}")
        shutil.copyfile(treefile, output_tree)
    else:
        cmd = ["FastTree", "-nt", "-gtr", "-quiet", "-out", output_tree, input_fasta]
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)


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
        cmd = ["FastTree", "-quiet", "-out", tree_out, filepath]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
        print(result.stdout.decode("utf-8"))

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
