#!/usr/bin/env python
"""Parse hmmsearch domtblout, build working df, extract per-marker ID lists."""
import argparse
import os
import shutil

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmout", required=True)
    parser.add_argument("--proteomes", required=True)
    parser.add_argument("--model_count", required=True, type=int)
    parser.add_argument("--percent_models", required=True, type=int)
    parser.add_argument("--lflt", required=True, type=int)
    parser.add_argument("--max_sdup", type=int, default=-1)
    parser.add_argument("--max_dupl", type=float, default=-1.0)
    parser.add_argument("--ref_merged_final", default=None)
    parser.add_argument("--ref_proteomes", default=None)
    args = parser.parse_args()

    hmmout = args.hmmout
    min_models_fraction = args.percent_models / 100
    lflt_fraction = args.lflt / 100

    # optional length filtering
    if lflt_fraction > 0:
        import subprocess
        df_lflt = pd.read_csv(hmmout, comment="#", sep=r'\s+', header=None)
        filtered_out = df_lflt[df_lflt[2] < df_lflt.groupby(3)[2].transform('median') * lflt_fraction][0]
        list_torm_path = "del.ls"
        filtered_out.to_csv(list_torm_path, index=False, header=False)
        out_lengthfilter = "hits.lfilt"
        with open(out_lengthfilter, "w") as l_out:
            subprocess.run(["grep", "-vwFf", list_torm_path, hmmout], stdout=l_out)
        shutil.move(out_lengthfilter, hmmout)

    # parse domtblout
    finaldf = pd.read_csv(hmmout, comment="#", sep=r'\s+', header=None, dtype=str)

    # count markers per genome
    genome_count = len(finaldf[0].apply(lambda x: x.split("|")[0]).unique())
    dict_counts = {}
    seen = set()
    for _, row in finaldf.iterrows():
        name = row[0].split("|")[0]
        model = row[3]
        proteinid = row[0]
        key = f"{proteinid}_{model}"
        if key not in seen:
            seen.add(key)
            if name not in dict_counts:
                dict_counts[name] = {}
            if model in dict_counts[name]:
                dict_counts[name][model] += 1
            else:
                dict_counts[name][model] = 1

    # filter incomplete genomes
    incomplete_genomes = {
        g for g in dict_counts
        if len(dict_counts[g]) < (args.model_count * min_models_fraction)
    }

    removed_reasons = {
        g: [f"minmarker:{len(dict_counts[g]) / args.model_count:.4f}"]
        for g in incomplete_genomes
    }
    removed_genomes = set(incomplete_genomes)

    if args.max_sdup >= 0:
        high_single_dup = {
            g for g, counts in dict_counts.items()
            if counts and max(counts.values()) > args.max_sdup
        }
        for g in high_single_dup:
            removed_reasons.setdefault(g, []).append(f"maxsdup:{max(dict_counts[g].values())}")
        removed_genomes.update(high_single_dup)

    if args.max_dupl >= 0:
        high_dup_fraction = set()
        for g, counts in dict_counts.items():
            dup_fraction = sum(1 for v in counts.values() if v > 1) / args.model_count
            if dup_fraction > args.max_dupl:
                high_dup_fraction.add(g)
                removed_reasons.setdefault(g, []).append(f"maxdupl:{dup_fraction:.4f}")
        removed_genomes.update(high_dup_fraction)

    rows_to_drop = []
    for idx, row in finaldf.iterrows():
        name = row[0].split("|")[0]
        if name in removed_genomes:
            rows_to_drop.append(idx)
    finaldf = finaldf.drop(rows_to_drop)

    with open("log_genomes_removed.txt", "w") as f:
        for genome in sorted(removed_genomes):
            reasons = ";".join(removed_reasons.get(genome, ["filtered"]))
            f.write(f"{genome}\t{reasons}\n")

    # write marker count matrix
    count_mat = pd.DataFrame.from_dict(dict_counts).fillna(0)
    count_mat.to_csv("marker_count_matrix.csv")

    # build working df
    finaldf_copy = finaldf.copy()
    finaldf_copy["savedname"] = finaldf_copy[0].apply(lambda c: c.replace("|", "/"))
    df = finaldf_copy.map(lambda x: x.split("|")[0] if isinstance(x, str) else x)
    df['namemodel'] = df[0] + "/" + df[3]
    df = df.drop_duplicates(subset='savedname', keep='first')

    # keep max 5 duplicates per namemodel
    helperdf = pd.concat(g for _, g in df.groupby("namemodel") if 5 >= len(g) >= 1)
    newdf = df.drop_duplicates(subset='namemodel', keep=False)
    df = pd.merge(helperdf, newdf, how='outer')

    # merge with reference if available
    if args.ref_merged_final is not None:
        df_ref = pd.read_csv(args.ref_merged_final)
        df = pd.concat([df, df_ref])

    df_fordups = df.set_index(df["savedname"])
    df_fordups.to_csv("table_elim_dups")

    # extract per-marker ID lists
    ls_seq_model = [
        (df.iloc[i]["savedname"].replace("/", "|"), df.iloc[i]["namemodel"].split("/")[1])
        for i in range(len(df))
    ]

    model_seqs = {}
    for seq_id, model in ls_seq_model:
        if model not in model_seqs:
            model_seqs[model] = []
        model_seqs[model].append(seq_id)

    os.makedirs("extracted", exist_ok=True)
    for model, seqs in model_seqs.items():
        with open(os.path.join("extracted", model), "w") as f:
            f.write("\n".join(seqs) + "\n")

    # build combined proteomes (query + reference)
    with open(args.proteomes) as fp:
        data = fp.read()
    if args.ref_proteomes is not None:
        with open(args.ref_proteomes) as fp:
            data += "\n" + fp.read()
    with open("combined_proteomes.faa", "w") as fp:
        fp.write(data)


if __name__ == "__main__":
    main()
