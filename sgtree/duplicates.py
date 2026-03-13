import os
import glob

import pandas as pd
from Bio import SeqIO

from sgtree.config import Config
from sgtree.parallel import map_threaded


SCORE_COLUMNS = ("score_bits", "7")


def _resolve_score_column(df_fordups) -> str:
    score_col = next((col for col in SCORE_COLUMNS if col in df_fordups.columns), None)
    if score_col is None:
        raise ValueError(
            f"Missing score column in duplicate table; expected one of: {', '.join(SCORE_COLUMNS)}"
        )
    return score_col


def _is_duplicate(seq_id, all_ids):
    """Check if a genome (first part of ID before |) appears more than once."""
    genome = seq_id.split("|")[0]
    return sum(1 for x in all_ids if x.split("|")[0] == genome) > 1


def _build_score_lookup(df_fordups, score_col: str) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in df_fordups.reset_index(drop=True).itertuples(index=False):
        key = str(row.savedname).replace("/", "|")
        lookup[key] = f"{row.savedname}:{float(getattr(row, score_col))}"
    return lookup


def _pick_best_scored_id(scored_ids: list[str]) -> str:
    """Pick highest-scoring identifier with deterministic tie-break on id."""
    if not scored_ids:
        raise ValueError("Expected at least one scored identifier")
    return min(
        scored_ids,
        key=lambda entry: (-float(entry.rsplit(":", 1)[1]), entry.rsplit(":", 1)[0]),
    )


def _process_file_worker(args):
    """Worker: eliminate duplicates for one aligned marker file."""
    filepath, aln_spectree_dir, score_lookup = args
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
            dups[key] = [score_lookup[v] for v in dups[key]]

        # for each set of duplicates, remove the best score from the "to-remove" list
        ids_to_remove = set()
        for key, scored_ids in dups.items():
            best_scored_id = _pick_best_scored_id(scored_ids)
            for scored_id in scored_ids:
                if scored_id != best_scored_id:
                    raw_id = scored_id.split(":", 1)[0].replace("/", "|")
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

def eliminate_duplicates(cfg: Config, df_fordups: pd.DataFrame):
    """For each aligned marker, keep only the highest-scoring hit per genome."""
    os.makedirs(cfg.aln_spectree_dir, exist_ok=True)
    score_col = _resolve_score_column(df_fordups)
    score_lookup = _build_score_lookup(df_fordups, score_col)

    aligned_files = glob.glob(os.path.join(cfg.aligned_dir, "*.faa"))
    args = [(f, cfg.aln_spectree_dir, score_lookup) for f in aligned_files]

    map_threaded(_process_file_worker, args, cfg.num_cpus)
