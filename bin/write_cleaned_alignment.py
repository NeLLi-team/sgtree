#!/usr/bin/env python
"""Write cleaned alignment based on nodes remaining in a newick tree."""
import argparse

from Bio import SeqIO
from ete3 import Tree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--newick", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    t = Tree(args.newick)
    lst_nodes = [node.name for node in t.traverse("postorder")]

    record_dict = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))

    with open(args.out, "w") as out:
        for key in list(record_dict.keys()):
            if key in lst_nodes:
                SeqIO.write(record_dict[key], out, "fasta")


if __name__ == "__main__":
    main()
