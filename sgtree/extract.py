import os
import glob
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO

from sgtree.config import Config


def extract_hits(cfg: Config, df: pd.DataFrame):
    """Extract hit identifiers per marker model, write per-model files to extracted/."""
    os.makedirs(cfg.extracted_dir, exist_ok=True)

    # build list of (sequence_id, model_name) pairs
    ls_seq_model = [
        (df.iloc[i]["savedname"].replace("/", "|"), df.iloc[i]["namemodel"].split("/")[1])
        for i in range(len(df))
    ]

    # create empty files for each model
    unique_models = set(pair[1] for pair in ls_seq_model)
    for model in unique_models:
        open(os.path.join(cfg.extracted_dir, model), "w").close()

    # write sequence IDs to their model files (grouped by model)
    model_seqs = {}
    for seq_id, model in ls_seq_model:
        if model not in model_seqs:
            model_seqs[model] = []
        model_seqs[model].append(seq_id)

    for model, seqs in model_seqs.items():
        with open(os.path.join(cfg.extracted_dir, model), "w") as f:
            f.write("\n".join(seqs) + "\n")


def _write_seqs_worker(args):
    """Worker: extract sequences for one model file (must be top-level for multiprocessing)."""
    filepath, ref_proteomes_path, extracted_seqs_dir = args
    try:
        with open(filepath) as f:
            ls_of_seq = f.read().strip().split("\n")

        fasta_parser = SeqIO.parse(ref_proteomes_path, "fasta")
        wanted = [rec for rec in fasta_parser if rec.id in ls_of_seq]

        hmm_match = os.path.basename(filepath)
        SeqIO.write(wanted, os.path.join(extracted_seqs_dir, hmm_match + ".faa"), "fasta")
    except Exception:
        import traceback
        print("Error for writing extracted sequences", traceback.format_exc())
        raise


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


def write_extracted_sequences(cfg: Config):
    """Retrieve actual protein sequences from proteome FASTA files, write per-model FASTAs."""
    os.makedirs(cfg.extracted_seqs_dir, exist_ok=True)

    # build combined proteomes file (query + reference)
    with open(cfg.proteomes_path) as fp:
        data = fp.read()

    if cfg.ref is not None:
        ref_proteomes = os.path.join(cfg.ref_dir_path(), "proteomes")
        with open(ref_proteomes) as fp:
            data += "\n" + fp.read()

    with open(cfg.ref_proteomes_path, "w") as fp:
        fp.write(data)

    # extract sequences for each model
    ls_of_files = glob.glob(os.path.join(cfg.extracted_dir, "*"))

    args = [
        (f, cfg.ref_proteomes_path, cfg.extracted_seqs_dir)
        for f in ls_of_files
    ]
    _map_with_fallback(_write_seqs_worker, args, cfg.num_cpus)
