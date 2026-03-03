#!/usr/bin/env python
"""Remove singleton markers poorly placed by comparing neighbor topology."""
import argparse

from ete3 import Tree


def score_func(ls_mod, ls_in):
    score = 0
    for i, init in enumerate(ls_in):
        if init in ls_mod:
            score += len(ls_in) - abs(ls_in.index(init) - ls_mod.index(init)) - i
    return score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--species_tree", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--num_nei", type=int, default=0, help="Override neighbor count (0=auto)")
    args = parser.parse_args()

    tf = Tree(args.tree)
    td = Tree(args.tree)

    # RF comparison is genome-level; deduplicate leaves that collapse to the
    # same genome id (can happen for reference genomes).
    seen_genomes = set()
    keep_leaves = []
    for leaf in td.iter_leaves():
        genome = leaf.name.split("|")[0]
        if genome not in seen_genomes:
            seen_genomes.add(genome)
            keep_leaves.append(leaf.name)
    td.prune(keep_leaves)

    td_leaves = list(td.iter_leaves())
    for node in td_leaves:
        node.name = node.name.split("|")[0]

    lst_nodes = [node for node in next(tf.copy().traverse())]
    ti = Tree(args.species_tree)
    ti.prune([n.name for n in td_leaves])

    rf, maxrf, *_ = ti.robinson_foulds(td, unrooted_trees=True)
    maxrf = maxrf + 0.0001
    rdist = rf / maxrf
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
    t_prot.write(format=1, outfile=args.out)


if __name__ == "__main__":
    main()
