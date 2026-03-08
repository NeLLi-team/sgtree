#!/usr/bin/env python
"""Remove singleton markers poorly placed by comparing neighbor topology."""
import argparse
import sys
from pathlib import Path

from ete3 import Tree

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.marker_selection import (
    _load_score_table,
    _tree_to_genome_level,
    choose_tree_by_rf,
    prune_singletons,
)


def score_func(ls_mod, ls_in):
    score = 0
    for i, init in enumerate(ls_in):
        if init in ls_mod:
            score += len(ls_in) - abs(ls_in.index(init) - ls_mod.index(init)) - i
    return score


def rf_distance_between(species_tree: Tree, marker_tree: Tree) -> float:
    species = species_tree.copy()
    marker = _tree_to_genome_level(marker_tree)
    species.prune([leaf.name for leaf in marker.iter_leaves()])
    rf, maxrf, *_ = species.robinson_foulds(marker, unrooted_trees=True)
    return rf / maxrf if maxrf else 0.0


def choose_tree_by_rf(species_tree: Tree, original_tree: Tree, candidate_tree: Tree) -> Tree:
    before = rf_distance_between(species_tree, original_tree)
    after = rf_distance_between(species_tree, candidate_tree)
    return candidate_tree if after < before else original_tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--species_tree", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--num_nei", type=int, default=0, help="Override neighbor count (0=auto)")
    parser.add_argument(
        "--singles_min_rfdist",
        type=float,
        default=0.25,
        help="Minimum marker-tree/global-tree RF distance required to activate singleton filtering",
    )
    parser.add_argument(
        "--mode",
        default="neighbor",
        choices=["neighbor", "delta_rf", "backbone", "ensemble"],
        help="singleton filtering mode",
    )
    parser.add_argument(
        "--table_elim_dups",
        default=None,
        help="optional duplicate-score table for ensemble singleton scoring",
    )
    args = parser.parse_args()

    tf = Tree(args.tree)
    td = _tree_to_genome_level(tf)
    td_leaves = list(td.iter_leaves())
    lst_nodes = [node for node in next(tf.copy().traverse())]
    ti = Tree(args.species_tree)
    ti.prune([n.name for n in td_leaves])

    rf, maxrf, *_ = ti.robinson_foulds(td, unrooted_trees=True)
    maxrf = maxrf + 0.0001
    rdist = rf / maxrf
    if rdist < args.singles_min_rfdist:
        tf.write(format=1, outfile=args.out)
        return
    if args.mode in {"delta_rf", "backbone", "ensemble"}:
        score_table = None
        score_col = None
        if args.mode == "ensemble":
            if args.table_elim_dups is None:
                raise ValueError("--table_elim_dups is required for --mode ensemble")
            score_table, score_col = _load_score_table(args.table_elim_dups)
        if args.num_nei > 0:
            num_nei = min(args.num_nei, max(1, len(lst_nodes) - 1))
        else:
            num_nei = max(2, min(5, len(lst_nodes) - 1))
        chosen = prune_singletons(
            species_tree=ti,
            working_tree=tf,
            mode=args.mode,
            k=num_nei,
            score_table=score_table,
            score_col=score_col,
        )
        choose_tree_by_rf(ti, tf, chosen).write(format=1, outfile=args.out)
        return
    max_neighbors = max(1, len(lst_nodes) - 1)
    if args.num_nei > 0:
        num_nei = min(args.num_nei, max_neighbors)
    else:
        num_nei = max(1, round(len(lst_nodes) * (1 - rdist)))
    total_score = num_nei ** 2
    cutoff = total_score / 15

    # get neighbors from protein tree
    dict_neighbors = {}
    for node in lst_nodes:
        ori_leaf = node
        while len(node.get_leaves()) - 1 < num_nei:
            node = node.up
        dict_neighbors[ori_leaf.name] = {}
        for leaf in node.get_leaves():
            if leaf.name != ori_leaf.name:
                dict_neighbors[ori_leaf.name][leaf.name] = ori_leaf.get_distance(leaf.name)
        dict_neighbors[ori_leaf.name] = sorted(
            dict_neighbors[ori_leaf.name],
            key=lambda leaf: (dict_neighbors[ori_leaf.name][leaf], leaf),
        )[:num_nei]

    # get neighbors from species tree
    lst_nodesi = [node for node in next(ti.copy().traverse())]
    dict_neighborsi = {}
    for nodei in lst_nodesi:
        ori_leafi = nodei
        while len(nodei.get_leaves()) - 1 < num_nei:
            nodei = nodei.up
        dict_neighborsi[ori_leafi.name] = {}
        for leaf in nodei.get_leaves():
            if leaf.name != ori_leafi.name:
                dict_neighborsi[ori_leafi.name][leaf.name] = ori_leafi.get_distance(leaf.name)
        dict_neighborsi[ori_leafi.name] = sorted(
            dict_neighborsi[ori_leafi.name],
            key=lambda leaf: (dict_neighborsi[ori_leafi.name][leaf], leaf),
        )[:num_nei]

    # score and filter
    flagged = []
    for node in lst_nodes:
        ls_model = [n.split("|")[0] for n in dict_neighbors[node.name]]
        ls_init = dict_neighborsi[node.name.split("|")[0]]
        score = score_func(ls_model, ls_init)
        if score > cutoff:
            flagged.append(node)

    flagged_names = [x.name for x in flagged]
    t_prot = tf.copy()
    t_prot.prune(flagged_names)
    choose_tree_by_rf(ti, tf, t_prot).write(format=1, outfile=args.out)


if __name__ == "__main__":
    main()
