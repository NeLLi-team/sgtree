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
    parser.add_argument(
        "--hmmsearch_cutoff",
        default="cut_ga",
        choices=["cut_ga", "cut_tc", "cut_nc", "evalue"],
        help="Threshold mode: model cutoffs or plain E-value",
    )
    parser.add_argument(
        "--hmmsearch_evalue",
        type=float,
        default=1e-5,
        help="E-value threshold when --hmmsearch_cutoff evalue",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with plan7.HMMFile(args.models) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise ValueError(f"No HMM profiles found in marker set: {args.models}")

    base_opts = {}
    if args.hmmsearch_cutoff == "cut_ga":
        base_opts["bit_cutoffs"] = "gathering"
    elif args.hmmsearch_cutoff == "cut_tc":
        base_opts["bit_cutoffs"] = "trusted"
    elif args.hmmsearch_cutoff == "cut_nc":
        base_opts["bit_cutoffs"] = "noise"
    else:
        base_opts["E"] = args.hmmsearch_evalue
        base_opts["domE"] = args.hmmsearch_evalue

    requested_cpus = max(1, args.cpus)

    def _run_search(cpus: int):
        search_opts = dict(base_opts)
        search_opts["cpus"] = cpus
        with easel.SequenceFile(args.proteomes, digital=True, alphabet=hmms[0].alphabet) as seq_file:
            with open(args.out, "wb") as out_handle:
                for i, hits in enumerate(hmmer.hmmsearch(hmms, seq_file, **search_opts)):
                    hits.write(out_handle, format="domains", header=(i == 0))

    try:
        _run_search(requested_cpus)
    except PermissionError as e:
        if requested_cpus == 1:
            raise
        print(
            f"warning: pyhmmer multiprocessing unavailable ({e}); "
            "retrying hmmsearch with cpus=1"
        )
        _run_search(1)


if __name__ == "__main__":
    main()
