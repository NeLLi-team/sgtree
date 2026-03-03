#!/usr/bin/env python
"""Extract sequences for a single marker from combined proteomes."""
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id_list", required=True)
    parser.add_argument("--proteomes", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.id_list) as f:
        wanted_ids = set(f.read().strip().split("\n"))

    wanted = [rec for rec in SeqIO.parse(args.proteomes, "fasta") if rec.id in wanted_ids]
    SeqIO.write(wanted, args.out, "fasta")


if __name__ == "__main__":
    main()
