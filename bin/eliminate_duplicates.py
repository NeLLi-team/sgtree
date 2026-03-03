#!/usr/bin/env python
"""Eliminate duplicate hits per genome for a single marker alignment."""
import argparse

import pandas as pd
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--table_elim_dups", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df_fordups = pd.read_csv(args.table_elim_dups).set_index("savedname")
    record_dict = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))
    all_ids = list(record_dict.keys())

    # find duplicates (same genome, multiple hits)
    dups = {}
    for key in all_ids:
        genome = key.split("|")[0]
        other_genomes = [k.split("|")[0] for k in all_ids if k != key]
        if genome in other_genomes:
            if genome in dups:
                dups[genome] = dups[genome] + "," + key
            else:
                dups[genome] = key

    # split and get scores
    for key in dups:
        dups[key] = dups[key].split(",")
        scored = []
        for v in dups[key]:
            lookup = v.replace("|", "/")
            row = df_fordups.loc[lookup]
            scored.append(f"{row.name}:{row.iloc[7]}")
        dups[key] = scored

    # for each set of duplicates, keep the best score
    ids_to_remove = set()
    for key, scored_ids in dups.items():
        scores = [float(s.split(":")[1]) for s in scored_ids]
        best_idx = scores.index(max(scores))
        for i, scored_id in enumerate(scored_ids):
            if i != best_idx:
                raw_id = scored_id.split(":")[0].replace("/", "|")
                ids_to_remove.add(raw_id)

    kept_ids = [k for k in all_ids if k not in ids_to_remove]

    with open(args.out, "w") as out:
        for seq_id in kept_ids:
            SeqIO.write(record_dict[seq_id], out, "fasta")


if __name__ == "__main__":
    main()
