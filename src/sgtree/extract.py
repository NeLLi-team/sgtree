"""Extract selected marker sequences from normalized proteome FASTA files."""

from __future__ import annotations

import shutil
from contextlib import ExitStack
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from sgtree.config import Config

_COPY_CHUNK = 1 << 20


def _concat_proteomes(source_paths: list[str], dest_path: str) -> None:
    """Concatenate proteome FASTA files into dest_path by streaming.

    A newline byte is written between consecutive sources so that a source
    which omits its trailing newline does not glue its last sequence to the
    next file's first header.
    """
    with Path(dest_path).open("wb") as destination:
        for source_index, source_path in enumerate(source_paths):
            if source_index > 0:
                destination.write(b"\n")
            with Path(source_path).open("rb") as source:
                shutil.copyfileobj(source, destination, _COPY_CHUNK)


def extract_hits(cfg: Config, df: pd.DataFrame) -> None:
    """Extract hit identifiers per marker model, write per-model files to extracted/."""
    extracted_dir = Path(cfg.extracted_dir)
    extracted_dir.mkdir(parents=True, exist_ok=True)

    # build list of (sequence_id, model_name) pairs
    ls_seq_model = [
        (
            df.iloc[i]["savedname"].replace("/", "|"),
            df.iloc[i]["namemodel"].split("/")[1],
        )
        for i in range(len(df))
    ]

    # create empty files for each model
    unique_models = {pair[1] for pair in ls_seq_model}
    for model in unique_models:
        (extracted_dir / model).touch()

    # write sequence IDs to their model files (grouped by model)
    model_seqs: dict[str, list[str]] = {}
    for seq_id, model in ls_seq_model:
        model_seqs.setdefault(model, []).append(seq_id)

    for model, seqs in model_seqs.items():
        with (extracted_dir / model).open("w") as handle:
            handle.write("\n".join(seqs) + "\n")


def write_extracted_sequences(cfg: Config) -> None:
    """Write selected protein sequences into one FASTA per marker model."""
    extracted_seqs_dir = Path(cfg.extracted_seqs_dir)
    extracted_seqs_dir.mkdir(parents=True, exist_ok=True)

    sources = [cfg.proteomes_path]
    if cfg.ref is not None:
        ref_dir = cfg.ref_dir_path()
        assert ref_dir is not None
        sources.append(str(Path(ref_dir) / "proteomes"))
    _concat_proteomes(sources, cfg.ref_proteomes_path)

    # Build id->models mapping from extracted marker ID lists.
    extracted_files = list(Path(cfg.extracted_dir).glob("*"))
    id_to_models: dict[str, list[str]] = {}
    models: list[str] = []
    for filepath in extracted_files:
        model = filepath.name
        models.append(model)
        with filepath.open() as handle:
            for line in handle:
                seq_id = line.strip()
                if not seq_id:
                    continue
                id_to_models.setdefault(seq_id, []).append(model)

    # Single-pass sequence extraction across the combined proteomes file.
    with ExitStack() as stack:
        handles = {
            model: stack.enter_context((extracted_seqs_dir / f"{model}.faa").open("w"))
            for model in models
        }
        for rec in SeqIO.parse(cfg.ref_proteomes_path, "fasta"):
            for model in id_to_models.get(rec.id, []):
                SeqIO.write(rec, handles[model], "fasta")
