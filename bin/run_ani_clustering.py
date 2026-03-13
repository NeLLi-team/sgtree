#!/usr/bin/env python
"""Run ANI clustering from staged genome manifests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.ani import run_ani_clustering


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query-manifest", required=True)
    parser.add_argument("--ref-manifest", default=None)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--ani-threshold", type=float, default=95.0)
    parser.add_argument("--ani-backend", default="auto")
    parser.add_argument("--ani-mcl-inflation", type=float, default=2.0)
    parser.add_argument("--num-cpus", type=int, default=1)
    args = parser.parse_args()

    run_ani_clustering(
        query_manifest=args.query_manifest,
        ref_manifest=args.ref_manifest,
        outdir=args.outdir,
        ani_threshold=args.ani_threshold / 100.0 if args.ani_threshold > 1.0 else args.ani_threshold,
        inflation=args.ani_mcl_inflation,
        backend=args.ani_backend,
        cpus=args.num_cpus,
    )


if __name__ == "__main__":
    main()
