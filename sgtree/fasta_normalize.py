"""Normalize proteome FASTA inputs to stable SGTree IDs."""

from __future__ import annotations

import glob
import io
import os
import re
from typing import Iterator

from Bio import SeqIO

VALID_AA = set("ABCDEFGHIKLMNPQRSTVWYBXZJUO")


def _sanitize_token(token: str, fallback: str) -> str:
    text = (token or "").strip().split()[0]
    text = text.replace("|", "_").replace("/", "_")
    text = re.sub(r"[^A-Za-z0-9._:-]+", "_", text).strip("_")
    return text or fallback


def _iter_fasta_records(path: str) -> Iterator[tuple[str, str]]:
    with open(path) as handle:
        text = handle.read()
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    # Repair malformed FASTA where headers were glued to sequence lines (for example '*>next_header').
    text = re.sub(r"(?<!^)(?<!\n)>", "\n>", text)
    if not text.endswith("\n"):
        text += "\n"
    for record in SeqIO.parse(io.StringIO(text), "fasta"):
        yield record.description, str(record.seq)


def _clean_sequence(sequence: str) -> tuple[str, int]:
    seq = re.sub(r"\s+", "", sequence).upper().replace("*", "")
    cleaned = []
    replaced = 0
    for aa in seq:
        if aa in VALID_AA:
            cleaned.append(aa)
        else:
            cleaned.append("X")
            replaced += 1
    return "".join(cleaned), replaced


def _iter_input_files(genomedir: str) -> list[str]:
    if os.path.isdir(genomedir):
        files = sorted(glob.glob(os.path.join(genomedir, "*.faa")))
        if not files:
            raise ValueError(f"No .faa files found in directory: {genomedir}")
        return files
    if not os.path.isfile(genomedir):
        raise FileNotFoundError(f"Input path does not exist: {genomedir}")
    return [genomedir]


def _normalize_ids(
    raw_header: str,
    file_genome: str,
    protein_index: int,
    infer_genome_from_header: bool,
) -> tuple[str, str]:
    token = raw_header.split()[0] if raw_header else ""
    if infer_genome_from_header and "|" in token:
        genome_part, protein_part = token.split("|", 1)
        genome_id = _sanitize_token(genome_part, file_genome)
    else:
        genome_id = file_genome
        protein_part = token.split("|", 1)[1] if "|" in token else token
    protein_id = _sanitize_token(protein_part, f"protein_{protein_index:06d}")
    return genome_id, protein_id


def normalize_and_concat_proteomes(
    genomedir: str,
    out_fasta: str,
    map_path: str | None = None,
) -> dict[str, int]:
    files = _iter_input_files(genomedir)
    infer_genome_from_header = os.path.isfile(genomedir)

    total_records = 0
    total_genomes = set()
    total_invalid = 0

    map_handle = open(map_path, "w") if map_path else None
    if map_handle:
        map_handle.write("source_file\toriginal_header\tnormalized_header\n")

    try:
        with open(out_fasta, "w") as out:
            for file_index, path in enumerate(files, start=1):
                stem = os.path.splitext(os.path.basename(path))[0]
                file_genome = _sanitize_token(stem, f"genome_{file_index:05d}")
                seen_proteins = set()

                for protein_index, (raw_header, raw_seq) in enumerate(_iter_fasta_records(path), start=1):
                    genome_id, protein_id = _normalize_ids(
                        raw_header, file_genome, protein_index, infer_genome_from_header
                    )
                    base_id = protein_id
                    suffix = 2
                    while protein_id in seen_proteins:
                        protein_id = f"{base_id}_{suffix}"
                        suffix += 1
                    seen_proteins.add(protein_id)

                    seq, replaced = _clean_sequence(raw_seq)
                    if not seq:
                        continue

                    normalized_id = f"{genome_id}|{protein_id}"
                    out.write(f">{normalized_id}\n{seq}\n")

                    total_records += 1
                    total_genomes.add(genome_id)
                    total_invalid += replaced

                    if map_handle:
                        header_clean = raw_header.replace("\t", " ").strip()
                        map_handle.write(f"{path}\t{header_clean}\t{normalized_id}\n")
    finally:
        if map_handle:
            map_handle.close()

    return {
        "genomes": len(total_genomes),
        "records": total_records,
        "invalid_chars_replaced": total_invalid,
    }
