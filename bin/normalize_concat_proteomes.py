#!/usr/bin/env python
"""Normalize and concatenate input proteomes into a SGTree-safe FASTA."""

import argparse
import sys
from pathlib import Path

# Ensure local package imports work when invoked as a script.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.fasta_normalize import normalize_and_concat_proteomes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input proteome directory (*.faa) or fasta file")
    parser.add_argument("--out", required=True, help="Normalized output FASTA")
    parser.add_argument("--map", default=None, help="Optional header mapping TSV")
    args = parser.parse_args()

    stats = normalize_and_concat_proteomes(args.input, args.out, args.map)
    print(
        "normalized proteomes: "
        f"genomes={stats['genomes']} records={stats['records']} "
        f"invalid_chars_replaced={stats['invalid_chars_replaced']}"
    )


if __name__ == "__main__":
    main()
