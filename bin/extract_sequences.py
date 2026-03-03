#!/usr/bin/env python
"""Extract sequences for a single marker from combined proteomes."""
import argparse

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--id_list", required=True)
    parser.add_argument("--proteomes", required=True)
    parser.add_argument("--index_db", default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.id_list) as f:
        wanted_ids = [line.strip() for line in f if line.strip()]

    if args.index_db:
        idx = SeqIO.index_db(args.index_db)
        try:
            with open(args.out, "w") as out:
                for seq_id in wanted_ids:
                    if seq_id in idx:
                        SeqIO.write(idx[seq_id], out, "fasta")
        finally:
            idx.close()
        return

    wanted_ids_set = set(wanted_ids)
    wanted = [rec for rec in SeqIO.parse(args.proteomes, "fasta") if rec.id in wanted_ids_set]
    SeqIO.write(wanted, args.out, "fasta")


if __name__ == "__main__":
    main()
