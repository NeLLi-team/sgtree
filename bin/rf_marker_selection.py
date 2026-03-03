#!/usr/bin/env python
"""RF-distance based duplicate resolution for one marker tree."""
import argparse

import pandas as pd
from ete3 import Tree


def get_ascore(identifier, table_path):
    dfa = pd.read_csv(table_path).set_index("savedname")
    row = dfa.loc[identifier.replace("|", "/")]
    return row.name + ":" + str(row.iloc[8])


def best_score(scored_list):
    scores = [float(s.split(":")[1]) for s in scored_list]
    best_idx = scores.index(max(scores))
    return scored_list[best_idx]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker_tree", required=True)
    parser.add_argument("--species_tree", required=True)
    parser.add_argument("--table_elim_dups", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--rf_out", required=True)
    parser.add_argument("--ref_list", default=None)
    args = parser.parse_args()

    ls_refs = args.ref_list.split(",") if args.ref_list else None

    t = Tree(args.marker_tree)
    lst_nodes = [node.name for node in next(t.copy().traverse())]

    # find duplicates
    seen = {}
    for x in lst_nodes:
        genome = x.split("|")[0]
        if genome not in seen:
            seen[genome] = x
        else:
            seen[genome] = seen[genome] + "," + x
    dups = {k: v.split(",") for k, v in seen.items()}

    # score all entries
    for key in dups:
        dups[key] = [get_ascore(v, args.table_elim_dups) for v in dups[key]]

    speciestree = Tree(args.species_tree)
    bad_nodes = []

    rf_lines = []
    for key, value in dups.items():
        if ls_refs is not None and key.split("|")[0] + ".faa" in ls_refs:
            continue
        if len(value) == 1:
            continue

        best_each = None
        best_protein = None
        rf_results = []
        candidates = []

        for each in value:
            ls_nodes_local = [each]
            dcopy = {k: v for k, v in dups.items() if k != key}
            for k, v in dcopy.items():
                bs = best_score(v)
                if bs is not None:
                    ls_nodes_local.append(bs)

            ls_nodes_local = [e for e in ls_nodes_local if e is not None]
            alt = [x.split(":")[0].replace("/", "|") for x in ls_nodes_local]

            speciestree_copy = speciestree.copy()
            speciestree_copy.prune([n.split("|")[0] for n in alt])

            t_prot = t.copy()
            t_prot.prune(alt)
            t_protcopy = t_prot.copy()

            for node in next(t_protcopy.traverse()):
                node.name = node.name.split("|")[0]

            rf, maxrf, *_ = speciestree_copy.robinson_foulds(t_protcopy, unrooted_trees=True)
            maxrf = maxrf + 0.0001
            rf_dist = rf / maxrf

            rf_results.append((each.split(":")[0], rf_dist))

            candidates.append((rf_dist, each))

        if candidates:
            # Deterministic tie-breaker: lowest RF, then lexicographically smallest protein id.
            _, best_each = min(candidates, key=lambda item: (item[0], item[1].split(":")[0]))
            best_protein = best_each.split(":")[0]

        for protein, rf_d in rf_results:
            status = "Kept" if protein == best_protein else "Removed"
            rf_lines.append(f"{protein} {args.marker} {rf_d:.6f} {status}")

        if best_each is not None:
            bad_nodes_local = [
                n.split(":")[0].replace("/", "|") for n in value if n != best_each
            ]
            bad_nodes.extend(bad_nodes_local)

    # write RF values
    with open(args.rf_out, "w") as f:
        for line in rf_lines:
            f.write(line + "\n")

    # write cleaned tree
    cleaned_nodes = [n for n in lst_nodes if n not in bad_nodes]
    t_final = t.copy()
    t_final.prune(cleaned_nodes)
    t_final.write(format=1, outfile=args.out)


if __name__ == "__main__":
    main()
