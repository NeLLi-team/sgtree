#!/usr/bin/env python
"""Run HMMER-style domain search using pyhmmer and write domtblout."""
import argparse

from pyhmmer import easel, hmmer, plan7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pyhmmer hmmsearch equivalent")
    parser.add_argument("--models", required=True, help="Combined marker-set HMM file")
    parser.add_argument("--proteomes", required=True, help="Combined proteomes FASTA file")
    parser.add_argument("--out", required=True, help="Output domtblout file")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs to use")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with plan7.HMMFile(args.models) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise ValueError(f"No HMM profiles found in marker set: {args.models}")

    with easel.SequenceFile(args.proteomes, digital=True, alphabet=hmms[0].alphabet) as seq_file:
        with open(args.out, "wb") as out_handle:
            for i, hits in enumerate(
                hmmer.hmmsearch(
                    hmms,
                    seq_file,
                    cpus=max(1, args.cpus),
                    bit_cutoffs="gathering",
                )
            ):
                hits.write(out_handle, format="domains", header=(i == 0))


if __name__ == "__main__":
    main()
