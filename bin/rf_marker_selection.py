#!/usr/bin/env python
"""RF-distance based duplicate resolution for one marker tree."""
import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ete3 import Tree

from sgtree.marker_selection import resolve_marker_tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker_tree", required=True)
    parser.add_argument("--species_tree", required=True)
    parser.add_argument("--table_elim_dups", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--rf_out", required=True)
    parser.add_argument("--ref_list", default=None)
    parser.add_argument("--selection_mode", default="coordinate", choices=["legacy", "coordinate"])
    parser.add_argument("--selection_max_rounds", type=int, default=5)
    parser.add_argument("--lock_references", default="no")
    args = parser.parse_args()

    ls_refs = args.ref_list.split(",") if args.ref_list else None
    lock_references = str(args.lock_references).strip().lower() in {"yes", "true", "1"}

    cleaned_nodes, records = resolve_marker_tree(
        marker_tree_path=args.marker_tree,
        species_tree_path=args.species_tree,
        table_path=args.table_elim_dups,
        marker_name=args.marker,
        ls_refs=ls_refs,
        selection_mode=args.selection_mode,
        max_rounds=args.selection_max_rounds,
        lock_references=lock_references,
    )

    with open(args.rf_out, "w") as handle:
        for record in records:
            handle.write(
                f"{record['protein_id']} {record['marker']} "
                f"{record['rf_distance']:.6f} {record['status']}\n"
            )

    tree = Tree(args.marker_tree)
    tree.prune(cleaned_nodes)
    tree.write(format=1, outfile=args.out)


if __name__ == "__main__":
    main()
