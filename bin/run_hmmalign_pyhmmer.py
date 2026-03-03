#!/usr/bin/env python
"""Run pyhmmer hmmalign for one marker and write aligned FASTA."""
import argparse

from pyhmmer import easel, hmmer, plan7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pyhmmer hmmalign for a marker")
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--model", help="Single-marker HMM file")
    model_group.add_argument("--models", help="Combined marker-set HMM file")
    parser.add_argument("--marker", help="Marker name when using --models")
    parser.add_argument("--seqs", required=True, help="Input marker FASTA sequences")
    parser.add_argument("--out", required=True, help="Output aligned FASTA file")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs to use")
    return parser.parse_args()


def _load_single_hmm(model_path: str):
    with plan7.HMMFile(model_path) as hmm_file:
        hmm_profile = next(iter(hmm_file), None)
    if hmm_profile is None:
        raise ValueError(f"No HMM profile found in file: {model_path}")
    return hmm_profile


def _load_marker_hmm(models_path: str, marker: str):
    with plan7.HMMFile(models_path) as hmm_file:
        for hmm_profile in hmm_file:
            hmm_name = hmm_profile.name
            if isinstance(hmm_name, bytes):
                hmm_name = hmm_name.decode("utf-8", errors="replace")
            if str(hmm_name) == marker:
                return hmm_profile
    raise ValueError(f"Marker '{marker}' not found in marker-set HMM file: {models_path}")


def main() -> None:
    args = parse_args()

    if args.model:
        hmm_profile = _load_single_hmm(args.model)
    else:
        if not args.marker:
            raise ValueError("--marker is required when using --models")
        hmm_profile = _load_marker_hmm(args.models, args.marker)
    with easel.SequenceFile(args.seqs, digital=True, alphabet=hmm_profile.alphabet) as seq_file:
        try:
            msa = hmmer.hmmalign(
                hmm_profile,
                seq_file,
                cpus=max(1, args.cpus),
                trim=True,
            )
        except (PermissionError, OSError):
            # Restricted runtimes can block pyhmmer multiprocessing primitives.
            seq_file.rewind()
            msa = hmmer.hmmalign(
                hmm_profile,
                seq_file,
                cpus=1,
                trim=True,
            )

    with open(args.out, "wb") as out_handle:
        msa.write(out_handle, format="afa")
    normalize_fasta(args.out)


def normalize_fasta(fasta_path: str) -> None:
    records = []
    header = None
    seq_chunks = []
    with open(fasta_path) as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = "|".join(line.split("|")[0:])
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append((header, "".join(seq_chunks)))

    records.sort(key=lambda item: item[0])
    with open(fasta_path, "w") as f:
        for header, seq in records:
            f.write(f"{header}\n{seq}\n")


if __name__ == "__main__":
    main()
