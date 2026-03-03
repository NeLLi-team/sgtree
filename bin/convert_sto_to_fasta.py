#!/usr/bin/env python
"""Convert Stockholm alignment to FASTA, cleaning headers."""
import argparse

from Bio import AlignIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sto", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    aln = AlignIO.read(args.sto, "stockholm")
    AlignIO.write(aln, args.out, "fasta")

    # clean up fasta headers (preserve full pipe-delimited ID)
    lines = []
    with open(args.out) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                lines.append("|".join(line.split("|")[0:]))
            else:
                lines.append(line)

    with open(args.out, "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
