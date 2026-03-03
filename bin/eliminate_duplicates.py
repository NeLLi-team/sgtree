#!/usr/bin/env python
"""Eliminate duplicate hits per genome for a single marker alignment."""
import argparse

import pandas as pd
from Bio import SeqIO


SCORE_COLUMNS = ("score_bits", "7")


def resolve_score_column(df_fordups: pd.DataFrame) -> str:
    score_col = next((col for col in SCORE_COLUMNS if col in df_fordups.columns), None)
    if score_col is None:
        raise ValueError(
            f"Missing score column in duplicate table; expected one of: {', '.join(SCORE_COLUMNS)}"
        )
    return score_col


def pick_best_scored_id(scored_ids: list[str]) -> str:
    if not scored_ids:
        raise ValueError("Expected at least one scored identifier")
    return min(
        scored_ids,
        key=lambda entry: (-float(entry.rsplit(":", 1)[1]), entry.rsplit(":", 1)[0]),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--table_elim_dups", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df_fordups = pd.read_csv(args.table_elim_dups).set_index("savedname")
    score_col = resolve_score_column(df_fordups)
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
            scored.append(f"{row.name}:{float(row[score_col])}")
        dups[key] = scored

    # for each set of duplicates, keep the best score
    ids_to_remove = set()
    for key, scored_ids in dups.items():
        best_scored_id = pick_best_scored_id(scored_ids)
        for scored_id in scored_ids:
            if scored_id != best_scored_id:
                raw_id = scored_id.split(":", 1)[0].replace("/", "|")
                ids_to_remove.add(raw_id)

    kept_ids = [k for k in all_ids if k not in ids_to_remove]

    with open(args.out, "w") as out:
        for seq_id in kept_ids:
            SeqIO.write(record_dict[seq_id], out, "fasta")


if __name__ == "__main__":
    main()
