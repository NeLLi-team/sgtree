import os
import glob

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

    # Build id->models mapping from extracted marker ID lists.
    ls_of_files = glob.glob(os.path.join(cfg.extracted_dir, "*"))
    id_to_models: dict[str, list[str]] = {}
    models: list[str] = []
    for filepath in ls_of_files:
        model = os.path.basename(filepath)
        models.append(model)
        with open(filepath) as f:
            for line in f:
                seq_id = line.strip()
                if not seq_id:
                    continue
                id_to_models.setdefault(seq_id, []).append(model)

    # Single-pass sequence extraction across the combined proteomes file.
    handles = {
        model: open(os.path.join(cfg.extracted_seqs_dir, model + ".faa"), "w")
        for model in models
    }
    try:
        for rec in SeqIO.parse(cfg.ref_proteomes_path, "fasta"):
            for model in id_to_models.get(rec.id, []):
                SeqIO.write(rec, handles[model], "fasta")
    finally:
        for handle in handles.values():
            handle.close()
