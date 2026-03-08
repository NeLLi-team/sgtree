from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from Bio import Phylo


def load_colors(color_file: str) -> dict[str, str]:
    color_dict: dict[str, str] = {}
    with open(color_file) as handle:
        for line in handle.readlines()[6:]:
            parts = line.strip().split(" ")
            if len(parts) >= 2:
                color_dict[parts[0]] = parts[1]
    return color_dict


def render_tree_file(tree_path: str, color_file: str, out_png: str) -> None:
    color_dict = load_colors(color_file)
    tree = Phylo.read(tree_path, "newick")

    terminals = tree.get_terminals()
    if len(terminals) > 2:
        try:
            tree.root_at_midpoint()
        except Exception:
            pass

    fig_height = max(2.5, min(0.32 * max(1, len(terminals)) + 1.5, 24))
    fig_width = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda clade: clade.name if clade.is_terminal() else None,
    )

    for text in ax.texts:
        label = text.get_text().strip()
        if label in color_dict:
            text.set_color(color_dict[label])
            text.set_fontweight("bold")

    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_tree(cfg, color_file: str):
    tree_path = os.path.join(cfg.outdir, "tree_final.nwk")
    if not os.path.exists(tree_path):
        return
    render_tree_file(tree_path, color_file, os.path.join(cfg.outdir, "tree_final.png"))
