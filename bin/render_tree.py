#!/usr/bin/env python
"""Render a Newick tree as PNG with color-coded leaf labels.

If ETE treeview dependencies are unavailable (common on headless clusters),
write a placeholder PNG and warning sidecar instead of failing the pipeline.
"""
import argparse
import base64
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

PLACEHOLDER_PNG_B64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO7+9pUAAAAASUVORK5CYII="
)


def write_placeholder_png(out_png: str, reason: str) -> None:
    with open(out_png, "wb") as handle:
        handle.write(base64.b64decode(PLACEHOLDER_PNG_B64))
    warning_path = f"{out_png}.warning.txt"
    with open(warning_path, "w") as handle:
        handle.write(f"{reason}\n")
    print(f"warning: {reason}", file=sys.stderr)
    print(f"warning: wrote placeholder image {out_png}", file=sys.stderr)
    print(f"warning: details in {warning_path}", file=sys.stderr)


def load_colors(path: str) -> dict[str, str]:
    color_dict = {}
    with open(path) as handle:
        for line in handle.readlines()[6:]:
            parts = line.strip().split(" ")
            if len(parts) >= 2:
                color_dict[parts[0]] = parts[1]
    return color_dict


def render_tree(tree_path: str, color_path: str, out_png: str) -> None:
    from ete3 import Tree, TreeStyle, faces

    color_dict = load_colors(color_path)
    tree = Tree(tree_path)

    def add_branchcolor(node):
        if not node.is_leaf():
            return
        color = color_dict.get(node.name, "#000000")
        taxon_face = faces.TextFace(node.name, fgcolor=color, fsize=10, fstyle="bold")
        faces.add_face_to_node(taxon_face, node, 1)
        node.img_style["fgcolor"] = color
        node.img_style["vt_line_color"] = color
        node.img_style["hz_line_color"] = color

    tree_style = TreeStyle()
    tree_style.mode = "r"
    tree_style.min_leaf_separation = 0
    tree_style.layout_fn = add_branchcolor

    if len(tree.get_leaves()) > 2:
        tree.set_outgroup(tree.get_midpoint_outgroup())

    tree.render(out_png, w=183, units="mm", tree_style=tree_style)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--color", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    try:
        render_tree(args.tree, args.color, args.out)
    except Exception as exc:
        write_placeholder_png(args.out, f"tree rendering skipped: {exc}")


if __name__ == "__main__":
    main()
