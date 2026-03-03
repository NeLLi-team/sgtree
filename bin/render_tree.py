#!/usr/bin/env python
"""Render a Newick tree as PNG with color-coded leaf labels."""
import argparse
import os

from ete3 import Tree, TreeStyle, faces

os.environ["QT_QPA_PLATFORM"] = "offscreen"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--color", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    color_dict = {}
    with open(args.color) as f:
        for line in f.readlines()[6:]:
            parts = line.strip().split(" ")
            if len(parts) >= 2:
                color_dict[parts[0]] = parts[1]

    t = Tree(args.tree)

    def add_branchcolor(node):
        if not node.is_leaf():
            return
        color = color_dict.get(node.name, "#000000")
        taxon_face = faces.TextFace(node.name, fgcolor=color, fsize=10, fstyle="bold")
        faces.add_face_to_node(taxon_face, node, 1)
        node.img_style["fgcolor"] = color
        node.img_style["vt_line_color"] = color
        node.img_style["hz_line_color"] = color

    ts = TreeStyle()
    ts.mode = "r"
    ts.min_leaf_separation = 0
    ts.layout_fn = add_branchcolor

    if len(t.get_leaves()) > 2:
        t.set_outgroup(t.get_midpoint_outgroup())

    t.render(args.out, w=183, units="mm", tree_style=ts)


if __name__ == "__main__":
    main()
