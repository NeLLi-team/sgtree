#!/usr/bin/env python
"""Split a combined HMM file into one profile file per marker."""
import argparse
import os

from pyhmmer import plan7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split marker-set HMM file into per-marker files")
    parser.add_argument("--models", required=True, help="Combined marker-set HMM file")
    parser.add_argument("--outdir", required=True, help="Output directory for per-marker .hmm files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    seen = set()
    with plan7.HMMFile(args.models) as hmm_file:
        for hmm_profile in hmm_file:
            name = hmm_profile.name
            if isinstance(name, bytes):
                name = name.decode("utf-8", errors="replace")
            marker = str(name)
            if marker in seen:
                raise ValueError(f"Duplicate marker name in HMM set: {marker}")
            seen.add(marker)
            out_path = os.path.join(args.outdir, f"{marker}.hmm")
            with open(out_path, "wb") as handle:
                hmm_profile.write(handle)

    if not seen:
        raise ValueError(f"No HMM profiles found in marker set: {args.models}")


if __name__ == "__main__":
    main()
