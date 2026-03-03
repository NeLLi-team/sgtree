import os
import subprocess
import glob
import multiprocessing as mp

from sgtree.config import Config


def run_fasttree(input_fasta: str, output_tree: str):
    """Run FastTree to build an approximate ML tree."""
    cmd = ["FastTree", "-quiet", "-out", output_tree, input_fasta]
    subprocess.run(cmd, stdout=subprocess.PIPE, check=True)


def _build_tree_worker(args):
    """Worker: build a single protein tree with FastTree."""
    filepath, treeout_dir = args
    tree_out = os.path.join(
        treeout_dir,
        os.path.basename(filepath) + "_tree.out",
    )
    cmd = ["FastTree", "-quiet", "-out", tree_out, filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    print(result.stdout.decode("utf-8"))


def run_fasttree_per_marker(cfg: Config, trimmed_dir: str, treeout_dir: str):
    """Build individual protein trees for each trimmed marker alignment."""
    os.makedirs(treeout_dir, exist_ok=True)

    files = glob.glob(os.path.join(trimmed_dir, "*"))
    args = [(f, treeout_dir) for f in files]

    pool = mp.Pool(4)
    pool.map(_build_tree_worker, args)
    pool.close()
    pool.join()
