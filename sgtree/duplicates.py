import os
import glob
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO

from sgtree.config import Config


def _is_duplicate(seq_id, all_ids):
    """Check if a genome (first part of ID before |) appears more than once."""
    genome = seq_id.split("|")[0]
    return sum(1 for x in all_ids if x.split("|")[0] == genome) > 1


def _get_score(identifier, df_fordups):
    """Look up the hmmsearch score for an identifier."""
    key = identifier.replace("|", "/")
    row = df_fordups.loc[key]
    return f"{row['savedname']}:{row.iloc[7]}"


def _process_file_worker(args):
    """Worker: eliminate duplicates for one aligned marker file."""
    filepath, aln_spectree_dir, df_fordups = args
    try:
        record_dict = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))
        all_ids = list(record_dict.keys())

        # find duplicates (same genome, multiple hits)
        dups = {}
        for key in all_ids:
            genome = key.split("|")[0]
            other_genomes = [k.split("|")[0] for k in all_ids if k != key]
            if genome in other_genomes:
                if genome in dups:
                    dups[genome] = dups[genome] + "," + key
                else:
                    dups[genome] = key

        # split comma-separated values and get scores
        for key in dups:
            dups[key] = dups[key].split(",")
            dups[key] = [_get_score(v, df_fordups) for v in dups[key]]

        # for each set of duplicates, remove the best score from the "to-remove" list
        ids_to_remove = set()
        for key, scored_ids in dups.items():
            scores = [float(s.split(":")[1]) for s in scored_ids]
            best_idx = scores.index(max(scores))
            for i, scored_id in enumerate(scored_ids):
                if i != best_idx:
                    raw_id = scored_id.split(":")[0].replace("/", "|")
                    ids_to_remove.add(raw_id)

        # keep only non-removed records
        kept_ids = [k for k in all_ids if k not in ids_to_remove]

        out_path = os.path.join(aln_spectree_dir, os.path.basename(filepath))
        with open(out_path, "w") as out:
            for seq_id in kept_ids:
                SeqIO.write(record_dict[seq_id], out, "fasta")

    except Exception as e:
        import sys
        print("elimination of duplicates exception", sys.exc_info())
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


def eliminate_duplicates(cfg: Config, df_fordups: pd.DataFrame):
    """For each aligned marker, keep only the highest-scoring hit per genome."""
    os.makedirs(cfg.aln_spectree_dir, exist_ok=True)

    aligned_files = glob.glob(os.path.join(cfg.aligned_dir, "*.faa"))
    args = [(f, cfg.aln_spectree_dir, df_fordups) for f in aligned_files]

    _map_with_fallback(_process_file_worker, args, cfg.num_cpus)
