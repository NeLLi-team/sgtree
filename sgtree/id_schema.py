from __future__ import annotations

import re


GENERIC_CONTIG_PREFIXES = {
    "protein",
    "proteins",
    "prot",
    "gene",
    "genes",
    "cds",
}


def sanitize_token(token: str, fallback: str) -> str:
    text = (token or "").strip().split()[0]
    text = text.replace("|", "_").replace("/", "_")
    text = re.sub(r"[^A-Za-z0-9._:-]+", "_", text).strip("_")
    return text or fallback


def infer_contig_id(token: str, *, fallback: str = "unknown_contig") -> tuple[str, str]:
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
    return "|".join(
        [
            sanitize_token(genome_id, "unknown_genome"),
            sanitize_token(contig_id, "unknown_contig"),
            sanitize_token(gene_id, "unknown_gene"),
        ]
    )


def parse_sequence_id(identifier: str) -> tuple[str, str, str]:
    parts = (identifier or "").split("|")
    if len(parts) >= 3:
        return parts[0], parts[1], "|".join(parts[2:])
    if len(parts) == 2:
        return parts[0], "unknown_contig", parts[1]
    if len(parts) == 1:
        return "unknown_genome", "unknown_contig", parts[0]
    return "unknown_genome", "unknown_contig", "unknown_gene"


def parse_savedname(savedname: str) -> tuple[str, str, str]:
    return parse_sequence_id(savedname.replace("/", "|"))
