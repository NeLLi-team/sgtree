"""Small FASTA utility helpers shared across the pipeline stages."""

from __future__ import annotations

import os
from pathlib import Path

from Bio import SeqIO


def fasta_contig_bases_stats(path: str | os.PathLike[str]) -> tuple[int, int]:
    """Return ``(contig_count, total_bases)`` over all records in a FASTA file.

    Accepts either ``str`` or ``os.PathLike``. Uses BioPython's ``SeqIO.parse``
    to remain consistent with the previous per-module implementations in
    ``sgtree/input_stage.py`` and ``sgtree/ani_clustering.py``.
    """
    path = Path(path)
    contigs = 0
    total_bases = 0
    with path.open() as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contigs += 1
            total_bases += len(record.seq)
    return contigs, total_bases
