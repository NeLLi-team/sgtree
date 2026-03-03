import os
import glob
import subprocess
import multiprocessing as mp

from pyhmmer import easel, hmmer, plan7

from sgtree.config import Config


def _run_mafft(args):
    extracted_seqs_dir, aligned_dir, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    aligned_dest = os.path.join(aligned_dir, filename)
    cmd = ["mafft", "--auto", "--thread", "4", "--quiet", filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    with open(aligned_dest, "w") as f:
        f.write(result.stdout.decode("utf-8") + "\n")


def _run_mafft_linsi(args):
    extracted_seqs_dir, aligned_dir, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    aligned_dest = os.path.join(aligned_dir, filename)
    cmd = ["mafft-linsi", "--auto", "--thread", "4", "--quiet", filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    with open(aligned_dest, "w") as f:
        f.write(result.stdout.decode("utf-8") + "\n")


def _run_hmmalign(args):
    extracted_seqs_dir, aligned_dir, modelset_path, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    model = filename.split(".")[0]
    faa_path = os.path.join(aligned_dir, model + ".faa")
    hmm_profile = None
    with plan7.HMMFile(modelset_path) as hmm_file:
        for hmm_profile_candidate in hmm_file:
            hmm_name = hmm_profile_candidate.name
            if isinstance(hmm_name, bytes):
                hmm_name = hmm_name.decode("utf-8", errors="replace")
            if str(hmm_name) == model:
                hmm_profile = hmm_profile_candidate
                break
    if hmm_profile is None:
        raise ValueError(f"Could not find marker '{model}' in marker-set HMM file: {modelset_path}")

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


def _map_with_fallback(func, args, workers: int):
    if not args:
        return
    n_workers = max(1, min(workers, len(args)))
    if n_workers == 1:
        for item in args:
            func(item)
        return
    try:
        with mp.Pool(n_workers) as pool:
            pool.map(func, args)
    except (PermissionError, OSError) as e:
        print(f"warning: multiprocessing unavailable ({e}); falling back to serial execution")
        for item in args:
            func(item)


def run_alignment(cfg: Config):
    """Run sequence alignment using the configured method (mafft, mafft-linsi, or hmmalign)."""
    os.makedirs(cfg.aligned_dir, exist_ok=True)

    files = [
        os.path.basename(f)
        for f in glob.glob(os.path.join(cfg.extracted_seqs_dir, "*"))
    ]

    print(f"- ...running {cfg.aln_method}")

    if cfg.aln_method == "mafft":
        _map_with_fallback(_run_mafft, [
            (cfg.extracted_seqs_dir, cfg.aligned_dir, f) for f in files
        ], cfg.num_cpus)
    elif cfg.aln_method == "mafft-linsi":
        _map_with_fallback(_run_mafft_linsi, [
            (cfg.extracted_seqs_dir, cfg.aligned_dir, f) for f in files
        ], cfg.num_cpus)
    else:
        _map_with_fallback(_run_hmmalign, [
            (cfg.extracted_seqs_dir, cfg.aligned_dir, cfg.models_path, f) for f in files
        ], cfg.num_cpus)
