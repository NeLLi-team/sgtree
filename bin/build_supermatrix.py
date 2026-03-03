#!/usr/bin/env python
"""Build concatenated alignment (supermatrix) from trimmed per-marker alignments."""
import argparse
import glob
import os

import pandas as pd
from Bio import SeqIO


def fill_nan_gaps(df_conc):
    """Replace NaN cells with X characters matching the alignment width."""
    # forward pass
    for j in range(1, df_conc.shape[1]):
        len_string = 0
        for i in range(df_conc.shape[0]):
            val = df_conc.iloc[i, j]
            if isinstance(val, float):
                if len_string > 0:
                    df_conc.iloc[i, j] = "X" * len_string
            else:
                len_string = len(str(val).replace("\n", ""))

    # backward pass
    for j in range(1, df_conc.shape[1]):
        if isinstance(df_conc.iloc[0, j], float):
            ref_len = 0
            for k in range(df_conc.shape[0] - 1, -1, -1):
                val = df_conc.iloc[k, j]
                if not isinstance(val, float):
                    ref_len = len(str(val).replace("\n", ""))
                    break
            for i in range(df_conc.shape[0]):
                if isinstance(df_conc.iloc[i, j], float):
                    df_conc.iloc[i, j] = "X" * ref_len
                else:
                    break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trimmed_dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args()

    df_conc = pd.DataFrame(columns=["SeqID"])
    for filepath in sorted(glob.glob(os.path.join(args.trimmed_dir, "*.faa"))):
        record_dict = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))
        record_dict = {k: v.format("fasta").split("\n", 1)[1] for k, v in record_dict.items()}
        new_dict = {}
        for key in record_dict:
            new_dict[key.split("|")[0]] = record_dict[key]
        new_df = pd.DataFrame(
            list(new_dict.items()),
            columns=["SeqID", os.path.basename(filepath)],
        )
        df_conc = pd.merge(new_df, df_conc, how="outer")

    marker_cols = sorted(col for col in df_conc.columns if col != "SeqID")
    df_conc = df_conc[["SeqID"] + marker_cols].sort_values("SeqID")

    fill_nan_gaps(df_conc)
    df_conc.to_csv(args.table)

    # rebuild from saved CSV and write concatenated FASTA
    df_conc = pd.read_csv(args.table)
    df_conc = df_conc.set_index("SeqID").sort_index()
    record_dict = df_conc.T.to_dict("list")
    record_dict = {k: v[1:] for k, v in record_dict.items()}
    record_dict = {k: "".join(str(x) for x in v) for k, v in record_dict.items()}
    record_dict = {k: v.replace("\n", "") for k, v in record_dict.items()}

    with open(args.out, "w") as f:
        for k in sorted(record_dict):
            v = record_dict[k]
            f.write(f">{k}\n{v}\n")


if __name__ == "__main__":
    main()
