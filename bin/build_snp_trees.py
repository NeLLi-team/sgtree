#!/usr/bin/env python
"""Build cluster-specific SNP trees from ANI clusters."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.ani import build_snp_trees


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--query-manifest", required=True)
    parser.add_argument("--ref-manifest", default=None)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-cluster-size", type=int, default=3)
    parser.add_argument("--num-cpus", type=int, default=1)
    parser.add_argument("--models", default=None)
    parser.add_argument("--hmmsearch-cutoff", default="cut_ga")
    parser.add_argument("--hmmsearch-evalue", type=float, default=1e-5)
    parser.add_argument("--min-contig-ani", type=float, default=0.95)
    args = parser.parse_args()

    build_snp_trees(
        clusters_path=args.clusters,
        query_manifest=args.query_manifest,
        ref_manifest=args.ref_manifest,
        outdir=args.outdir,
        min_cluster_size=args.min_cluster_size,
        cpus=args.num_cpus,
        models_path=args.models,
        hmmsearch_cutoff=args.hmmsearch_cutoff,
        hmmsearch_evalue=args.hmmsearch_evalue,
        min_contig_ani=args.min_contig_ani,
    )


if __name__ == "__main__":
    main()
