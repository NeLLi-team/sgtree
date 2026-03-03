#!/usr/bin/env python
"""Run pyhmmer hmmalign for one marker and write aligned FASTA."""
import argparse

from pyhmmer import easel, hmmer, plan7


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pyhmmer hmmalign for a marker")
    parser.add_argument("--models", required=True, help="Combined marker-set HMM file")
    parser.add_argument("--marker", required=True, help="Marker name to align against")
    parser.add_argument("--seqs", required=True, help="Input marker FASTA sequences")
    parser.add_argument("--out", required=True, help="Output aligned FASTA file")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs to use")
    return parser.parse_args()


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

    hmm_profile = _load_marker_hmm(args.models, args.marker)
    with easel.SequenceFile(args.seqs, digital=True, alphabet=hmm_profile.alphabet) as seq_file:
        msa = hmmer.hmmalign(
            hmm_profile,
            seq_file,
            cpus=max(1, args.cpus),
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
