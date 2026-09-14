"""Create and parse stable SGTree genome, contig, and sequence identifiers."""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

GENERIC_CONTIG_PREFIXES = {
    "protein",
    "proteins",
    "prot",
    "gene",
    "genes",
    "cds",
}


def sanitize_token(token: str, fallback: str) -> str:
    """Return one delimiter-safe identifier token or the supplied fallback."""
    parts = (token or "").strip().split()
    if not parts:
        return fallback
    text = parts[0]
    text = text.replace("|", "_").replace("/", "_")
    text = re.sub(r"[^A-Za-z0-9._:-]+", "_", text).strip("_")
    return text or fallback


def assign_genome_ids(paths: Sequence[str]) -> list[tuple[str, str]]:
    """Assign unique genome IDs from input filenames."""
    assigned: list[tuple[str, str]] = []
    sources_by_genome: dict[str, str] = {}
    for file_index, path in enumerate(paths, start=1):
        stem = Path(path).stem
        genome_id = sanitize_token(stem, f"genome_{file_index:05d}")
        previous_source = sources_by_genome.get(genome_id)
        if previous_source is not None:
            raise ValueError(
                f"Input filenames normalize to duplicate genome ID {genome_id!r}: "
                f"{previous_source!r} and {path!r}"
            )
        sources_by_genome[genome_id] = path
        assigned.append((path, genome_id))
    return assigned


def infer_contig_id(token: str, *, fallback: str = "unknown_contig") -> tuple[str, str]:
    """Infer a contig ID from a protein token and report the inference source."""
    text = sanitize_token(token, fallback)
    if text == fallback:
        return fallback, "unknown"

    for pattern in (
        r"^(.+)_\d+$",
        r"^(.+)-\d+$",
        r"^(.+):\d+$",
    ):
        match = re.match(pattern, text)
        if match:
            candidate = sanitize_token(match.group(1), fallback)
            if candidate.lower() not in GENERIC_CONTIG_PREFIXES:
                return candidate, "suffix"

    return fallback, "unknown"


def build_sequence_id(genome_id: str, contig_id: str, gene_id: str) -> str:
    """Build a normalized ``genome|contig|gene`` identifier."""
    return "|".join(
        [
            sanitize_token(genome_id, "unknown_genome"),
            sanitize_token(contig_id, "unknown_contig"),
            sanitize_token(gene_id, "unknown_gene"),
        ]
    )


def parse_sequence_id(identifier: str) -> tuple[str, str, str]:
    """Parse an SGTree sequence identifier, filling absent legacy fields."""
    parts = (identifier or "").split("|")
    if len(parts) >= 3:
        return parts[0], parts[1], "|".join(parts[2:])
    if len(parts) == 2:
        return parts[0], "unknown_contig", parts[1]
    if len(parts) == 1:
        return "unknown_genome", "unknown_contig", parts[0]
    return "unknown_genome", "unknown_contig", "unknown_gene"


def parse_savedname(savedname: str) -> tuple[str, str, str]:
    """Parse the slash-delimited identifier stored in SGTree tables."""
    return parse_sequence_id(savedname.replace("/", "|"))
