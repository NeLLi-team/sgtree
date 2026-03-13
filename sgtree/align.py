import os
import glob
import subprocess
import shutil

from pyhmmer import easel, hmmer, plan7

from sgtree.config import Config
from sgtree.parallel import map_processed, map_threaded


def _run_mafft(args):
    extracted_seqs_dir, aligned_dir, filename, threads = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    aligned_dest = os.path.join(aligned_dir, filename)
    cmd = ["mafft", "--auto", "--thread", str(threads), "--quiet", filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    with open(aligned_dest, "w") as f:
        f.write(result.stdout.decode("utf-8") + "\n")


def _run_mafft_linsi(args):
    extracted_seqs_dir, aligned_dir, filename, threads = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    aligned_dest = os.path.join(aligned_dir, filename)
    cmd = ["mafft-linsi", "--auto", "--thread", str(threads), "--quiet", filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    with open(aligned_dest, "w") as f:
        f.write(result.stdout.decode("utf-8") + "\n")


def _run_hmmalign(args):
    extracted_seqs_dir, aligned_dir, model_hmm_path, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    model = filename.split(".")[0]
    faa_path = os.path.join(aligned_dir, model + ".faa")
    with plan7.HMMFile(model_hmm_path) as hmm_file:
        hmm_profile = next(iter(hmm_file), None)
    if hmm_profile is None:
        raise ValueError(f"Could not read HMM profile from {model_hmm_path}")

    with easel.SequenceFile(filepath, digital=True, alphabet=hmm_profile.alphabet) as seq_file:
        msa = hmmer.hmmalign(hmm_profile, seq_file, cpus=1, trim=True)
    with open(faa_path, "wb") as out_handle:
        msa.write(out_handle, format="afa")
    _normalize_fasta(faa_path)


def _normalize_fasta(fasta_path: str) -> None:
    """Normalize FASTA headers and sequence order for deterministic downstream trees."""
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

def _split_models(models_path: str, split_dir: str) -> None:
    if os.path.isdir(split_dir):
        shutil.rmtree(split_dir)
    os.makedirs(split_dir, exist_ok=True)
    with plan7.HMMFile(models_path) as hmm_file:
        for hmm_profile in hmm_file:
            marker = hmm_profile.name
            if isinstance(marker, bytes):
                marker = marker.decode("utf-8", errors="replace")
            out_path = os.path.join(split_dir, f"{marker}.hmm")
            with open(out_path, "wb") as out_handle:
                hmm_profile.write(out_handle)


def run_alignment(
    cfg: Config,
    *,
    extracted_seqs_dir: str | None = None,
    aligned_dir: str | None = None,
):
    """Run sequence alignment using the configured method.

    By default this uses the standard extracted/aligned directories from the
    config, but marker-selection cleanup can provide alternate directories when
    alignments need to be rebuilt from cleaned sequence sets.
    """
    extracted_seqs_dir = extracted_seqs_dir or cfg.extracted_seqs_dir
    aligned_dir = aligned_dir or cfg.aligned_dir
    os.makedirs(aligned_dir, exist_ok=True)

    files = [
        os.path.basename(f)
        for f in glob.glob(os.path.join(extracted_seqs_dir, "*"))
    ]

    print(f"- ...running {cfg.aln_method}")
    n_jobs = max(1, min(cfg.num_cpus, len(files) if files else 1))
    threads_per_job = max(1, cfg.num_cpus // n_jobs)

    if cfg.aln_method == "mafft":
        map_threaded(_run_mafft, [
            (extracted_seqs_dir, aligned_dir, f, threads_per_job) for f in files
        ], n_jobs)
    elif cfg.aln_method == "mafft-linsi":
        map_threaded(_run_mafft_linsi, [
            (extracted_seqs_dir, aligned_dir, f, threads_per_job) for f in files
        ], n_jobs)
    else:
        split_dir = os.path.join(cfg.outdir, "models_split")
        _split_models(cfg.models_path, split_dir)
        map_processed(_run_hmmalign, [
            (extracted_seqs_dir, aligned_dir, os.path.join(split_dir, f"{f.split('.')[0]}.hmm"), f)
            for f in files
        ], n_jobs)
