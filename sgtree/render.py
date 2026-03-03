import os

from ete3 import Tree, TreeStyle, faces


def render_tree(cfg, color_file: str):
    """Render tree_final.nwk as PNG with color-coded branches."""
    tree_path = os.path.join(cfg.outdir, "tree_final.nwk")
    if not os.path.exists(tree_path):
        return

    # parse color file
    color_dict = {}
    try:
        with open(color_file) as f:
            content = f.readlines()
            for line in content[6:]:
                parts = line.strip().split(" ")
                if len(parts) >= 2:
                    color_dict[parts[0]] = parts[1]
    except Exception:
        print("ERROR creating colordict")
        return

    def add_branchcolor(node):
        if node.is_leaf():
            color = color_dict.get(node.name, "#000000")
            taxon_face = faces.TextFace(node.name, fgcolor=color, fsize=10, fstyle="bold")
            faces.add_face_to_node(taxon_face, node, 1)
            node.img_style["fgcolor"] = color
            node.img_style["vt_line_color"] = color
            node.img_style["hz_line_color"] = color

    tfin = Tree(tree_path)
    tfin.dist = 0
    tfin.allow_face_overlap = True
    ts = TreeStyle()
    ts.mode = "r"
    ts.min_leaf_separation = 0
    ts.layout_fn = add_branchcolor
    tfin.set_outgroup(tfin.get_midpoint_outgroup())
    tfin.render(os.path.join(cfg.outdir, "tree_final.png"), w=183, units="mm", tree_style=ts)
