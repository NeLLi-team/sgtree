#!/usr/bin/env python
"""Generate iTOL heatmap annotation file from marker count matrix."""
import argparse

import pandas as pd
from ete3 import Tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    itol_tree = Tree(args.tree)
    lst_nodes = [node for node in next(itol_tree.copy().traverse())]
    treetaxa = [n.name for n in lst_nodes]

    mkr_count_df = pd.read_csv(args.matrix)
    ls_order_cols = list(mkr_count_df["Unnamed: 0"])
    target_list_d = mkr_count_df.set_index("Unnamed: 0").to_dict()
    target_list_d = {k: list(v.values()) for k, v in target_list_d.items()}

    with open(args.out, "w") as f:
        f.write(
            "DATASET_HEATMAP\n"
            "SEPARATOR COMMA\n"
            "DATASET_LABEL,conservedOGs\n"
            "COLOR,#ff0000\n"
            "COLOR_BRANCHES,0\n"
            "LEGEND_TITLE,conservedOGs\n"
            f"FIELD_LABELS,{','.join(ls_order_cols)}\n"
            "MARGIN,0\n"
            "STRIP_WIDTH,25\n"
            "SHOW_INTERNAL,0\n"
            "COLOR_NAN,#000000\n"
            "AUTO_LEGEND,1\n"
            "COLOR_MIN,#FFFFFF\n"
            "COLOR_MAX,#b3b3b3\n"
            "MAXIMUM_SIZE,10\n"
            "DASHED_LINES,1\n"
            "BORDER_WIDTH,0\n"
            "BORDER_COLOR,#0000ff\n"
            "DATA\n"
        )
        for taxon in treetaxa:
            if taxon in target_list_d:
                values = ",".join(str(v) for v in target_list_d[taxon])
                f.write(f"{taxon},{values}\n")


if __name__ == "__main__":
    main()
