import os
import glob

import pandas as pd
from ete3 import Tree
from Bio import SeqIO

from sgtree.config import Config


def _ids_from_input(path: str) -> list[str]:
    if os.path.isdir(path):
        return [os.path.basename(f).split(".")[0] for f in glob.glob(os.path.join(path, "*"))]

    genome_ids = []
    with open(path) as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            gid = rec.id.split("|")[0]
            if gid not in genome_ids:
                genome_ids.append(gid)
    return genome_ids


def write_color_file(cfg: Config):
    """Write iTOL color strip annotation file."""
    filepath = os.path.join(cfg.outdir, "color.txt")
    with open(filepath, "w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR SPACE\n")
        f.write("DATASET_LABEL label1\n")
        f.write("COLOR #ff0000\n")
        f.write("COLOR_BRANCHES 0\n")
        f.write("DATA\n")
        for gid in _ids_from_input(cfg.genomedir):
            f.write(f"{gid} #FF0000\n")
        if cfg.ref is not None:
            for gid in _ids_from_input(cfg.ref):
                f.write(f"{gid} #C0C0C0\n")


def write_heatmap(cfg: Config, tree_path: str, outsuffix: str):
    """Generate iTOL heatmap annotation file from marker count matrix."""
    itol_tree = Tree(tree_path)
    lst_nodes = [node for node in next(itol_tree.copy().traverse())]
    treetaxa = [n.name for n in lst_nodes]

    mkr_count_df = pd.read_csv(os.path.join(cfg.outdir, "marker_count_matrix.csv"))
    ls_order_cols = list(mkr_count_df["Unnamed: 0"])
    target_list_d = mkr_count_df.set_index("Unnamed: 0").to_dict()
    target_list_d = {k: list(v.values()) for k, v in target_list_d.items()}

    outpath = os.path.join(cfg.outdir, outsuffix)
    with open(outpath, "w") as f:
        f.write(
            "DATASET_HEATMAP\n"
            "SEPARATOR COMMA\n"
            f"DATASET_LABEL,conservedOGs\n"
            "COLOR,#ff0000\n"
            "COLOR_BRANCHES,0\n"
            f"LEGEND_TITLE,conservedOGs\n"
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
