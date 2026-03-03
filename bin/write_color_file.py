#!/usr/bin/env python
"""Write iTOL color strip annotation file from query and optional reference proteomes."""
import argparse
import glob
import os

from Bio import SeqIO


def _ids_from_input(path: str):
    if os.path.isdir(path):
        return [os.path.basename(f).split(".")[0] for f in glob.glob(os.path.join(path, "*"))]

    genome_ids = []
    with open(path) as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            gid = rec.id.split("|")[0]
            if gid not in genome_ids:
                genome_ids.append(gid)
    return genome_ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomedir", required=True)
    parser.add_argument("--refdir", default="")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    query_ids = _ids_from_input(args.genomedir)
    ref_ids = _ids_from_input(args.refdir) if args.refdir else []

    with open(args.out, "w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR SPACE\n")
        f.write("DATASET_LABEL label1\n")
        f.write("COLOR #ff0000\n")
        f.write("COLOR_BRANCHES 0\n")
        f.write("DATA\n")
        for gid in query_ids:
            f.write(f"{gid} #FF0000\n")
        for gid in ref_ids:
            f.write(f"{gid} #C0C0C0\n")


if __name__ == "__main__":
    main()
