#!/usr/bin/env python
"""Regression parity checks for modular Python and Nextflow SGTree workflows."""
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from ete3 import Tree


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "regression"


def run(cmd):
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=ROOT, check=True)


def assert_exists(path: Path):
    if not path.exists():
        raise AssertionError(f"Missing expected output: {path}")


def assert_any(paths):
    for path in paths:
        if path.exists():
            return
    rendered = ", ".join(str(p) for p in paths)
    raise AssertionError(f"Missing expected output (any of): {rendered}")


def rf_norm(tree_a: Path, tree_b: Path) -> float:
    t1 = Tree(str(tree_a))
    t2 = Tree(str(tree_b))
    rf, maxrf, *_ = t1.robinson_foulds(t2, unrooted_trees=True)
    return (rf / maxrf) if maxrf else 0.0


def load_rf_status(path: Path) -> dict:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    return {
        (row.ProteinID, row.MarkerGene): row.Status
        for row in df.itertuples(index=False)
    }


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    return df.sort_index().sort_index(axis=1).fillna(0)


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    py_basic = OUT / "py_basic"
    py_full = OUT / "py_full"
    py_full_singles = OUT / "py_full_singles"
    nf_basic = OUT / "nf_basic"
    nf_full = OUT / "nf_full"
    nf_full_singles = OUT / "nf_full_singles"
    logs_dir = OUT / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    run([
        "python", "-m", "sgtree",
        "testgenomes/Chloroflexi", "hmms/UNI56",
        "--num_cpus", "8",
        "--save_dir", str(py_basic),
    ])
    run([
        "python", "-m", "sgtree",
        "testgenomes/Chloroflexi", "hmms/UNI56",
        "--num_cpus", "8",
        "--marker_selection", "yes",
        "--ref", "testgenomes/chlorref",
        "--save_dir", str(py_full),
    ])
    run([
        "python", "-m", "sgtree",
        "testgenomes/Chloroflexi", "hmms/UNI56",
        "--num_cpus", "8",
        "--marker_selection", "yes",
        "--singles", "yes",
        "--ref", "testgenomes/chlorref",
        "--save_dir", str(py_full_singles),
    ])
    run([
        "./nextflow", "-log", str(logs_dir / "nf_basic.log"), "run", "main.nf",
        "-work-dir", str(OUT / "work_basic"),
        "--genomedir", "testgenomes/Chloroflexi",
        "--modeldir", "hmms/UNI56",
        "--outdir", str(nf_basic),
    ])
    run([
        "./nextflow", "-log", str(logs_dir / "nf_full.log"), "run", "main.nf",
        "-work-dir", str(OUT / "work_full"),
        "--genomedir", "testgenomes/Chloroflexi",
        "--modeldir", "hmms/UNI56",
        "--outdir", str(nf_full),
        "--marker_selection", "true",
        "--ref", "testgenomes/chlorref",
    ])
    run([
        "./nextflow", "-log", str(logs_dir / "nf_full_singles.log"), "run", "main.nf",
        "-work-dir", str(OUT / "work_full_singles"),
        "--genomedir", "testgenomes/Chloroflexi",
        "--modeldir", "hmms/UNI56",
        "--outdir", str(nf_full_singles),
        "--marker_selection", "true",
        "--singles", "yes",
        "--ref", "testgenomes/chlorref",
    ])

    # Expected output files
    assert_exists(py_basic / "tree.nwk")
    assert_exists(py_basic / "marker_count_matrix.csv")
    assert_any([py_basic / "marker_count.txt", py_basic / "temp" / "itol" / "marker_count.txt"])

    assert_exists(py_full / "tree_final.nwk")
    assert_exists(py_full / "tree_final.png")
    assert_exists(py_full / "marker_selection_rf_values.txt")
    assert_exists(py_full / "marker_count_matrix.csv")
    assert_any([py_full / "marker_counts.txt", py_full / "temp" / "itol" / "marker_counts.txt"])
    assert_exists(py_full_singles / "tree_final.nwk")
    assert_exists(py_full_singles / "tree_final.png")
    assert_exists(py_full_singles / "marker_selection_rf_values.txt")
    assert_exists(py_full_singles / "marker_count_matrix.csv")
    assert_any([py_full_singles / "marker_counts.txt", py_full_singles / "temp" / "itol" / "marker_counts.txt"])

    assert_exists(nf_basic / "tree.nwk")
    assert_exists(nf_basic / "marker_count_matrix.csv")
    assert_exists(nf_basic / "marker_count.txt")
    assert_exists(nf_basic / "color.txt")

    assert_exists(nf_full / "tree_final.nwk")
    assert_exists(nf_full / "tree_final.png")
    assert_exists(nf_full / "marker_selection_rf_values.txt")
    assert_exists(nf_full / "marker_count_matrix.csv")
    assert_exists(nf_full / "marker_counts.txt")
    assert_exists(nf_full / "color.txt")
    assert_exists(nf_full_singles / "tree_final.nwk")
    assert_exists(nf_full_singles / "tree_final.png")
    assert_exists(nf_full_singles / "marker_selection_rf_values.txt")
    assert_exists(nf_full_singles / "marker_count_matrix.csv")
    assert_exists(nf_full_singles / "marker_counts.txt")
    assert_exists(nf_full_singles / "color.txt")

    # Topology parity checks
    if rf_norm(py_basic / "tree.nwk", nf_basic / "tree.nwk") != 0:
        raise AssertionError("Topology mismatch: python basic vs nextflow basic")
    if rf_norm(py_full / "tree_final.nwk", nf_full / "tree_final.nwk") != 0:
        raise AssertionError("Topology mismatch: python full vs nextflow full")
    if rf_norm(py_full_singles / "tree_final.nwk", nf_full_singles / "tree_final.nwk") != 0:
        raise AssertionError("Topology mismatch: python singles vs nextflow singles")

    # Matrix parity checks (ignore row/column ordering)
    if not load_matrix(py_basic / "marker_count_matrix.csv").equals(load_matrix(nf_basic / "marker_count_matrix.csv")):
        raise AssertionError("Marker matrix mismatch: python basic vs nextflow basic")
    if not load_matrix(py_full / "marker_count_matrix.csv").equals(load_matrix(nf_full / "marker_count_matrix.csv")):
        raise AssertionError("Marker matrix mismatch: python full vs nextflow full")
    if not load_matrix(py_full_singles / "marker_count_matrix.csv").equals(load_matrix(nf_full_singles / "marker_count_matrix.csv")):
        raise AssertionError("Marker matrix mismatch: python singles vs nextflow singles")

    # RF status parity checks
    py_rf = load_rf_status(py_full / "marker_selection_rf_values.txt")
    nf_rf = load_rf_status(nf_full / "marker_selection_rf_values.txt")
    if py_rf != nf_rf:
        only_py = sorted(set(py_rf.items()) - set(nf_rf.items()))
        only_nf = sorted(set(nf_rf.items()) - set(py_rf.items()))
        raise AssertionError(
            "RF status mismatch between python and nextflow.\n"
            f"Only in python: {only_py[:5]}\n"
            f"Only in nextflow: {only_nf[:5]}"
        )

    py_rf_singles = load_rf_status(py_full_singles / "marker_selection_rf_values.txt")
    nf_rf_singles = load_rf_status(nf_full_singles / "marker_selection_rf_values.txt")
    if py_rf_singles != nf_rf_singles:
        only_py = sorted(set(py_rf_singles.items()) - set(nf_rf_singles.items()))
        only_nf = sorted(set(nf_rf_singles.items()) - set(py_rf_singles.items()))
        raise AssertionError(
            "RF status mismatch between python and nextflow (singles).\n"
            f"Only in python: {only_py[:5]}\n"
            f"Only in nextflow: {only_nf[:5]}"
        )

    print("Regression parity checks passed.")


if __name__ == "__main__":
    main()
