#!/usr/bin/env python
"""Thin CLI wrapper around sgtree.supermatrix.build_supermatrix."""
import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from sgtree.supermatrix import build_supermatrix  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trimmed_dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args()

    out_dir = str(Path(args.out).resolve().parent)
    build_supermatrix(args.trimmed_dir, out_dir, args.table, args.out)


if __name__ == "__main__":
    main()
