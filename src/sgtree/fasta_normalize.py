"""Normalize proteome FASTA inputs to stable SGTree IDs."""

from __future__ import annotations

import io
import re
from collections.abc import Iterator
from contextlib import ExitStack
from pathlib import Path

from Bio import SeqIO

from sgtree.id_schema import (
    assign_genome_ids,
    build_sequence_id,
    infer_contig_id,
    sanitize_token,
)
from sgtree.input_stage import _list_files

VALID_AA = set("ABCDEFGHIKLMNPQRSTVWYBXZJUO")


def _iter_fasta_records(path: str) -> Iterator[tuple[str, str]]:
    with Path(path).open() as handle:
        text = handle.read()
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    # Repair malformed FASTA where a header is glued to the prior sequence.
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


def _normalize_ids(
    raw_header: str,
    file_genome: str,
    protein_index: int,
    *,
    infer_genome_from_header: bool,
) -> tuple[str, str, str, str]:
    token = raw_header.split(maxsplit=1)[0] if raw_header else ""
    if "|" in token:
        parts = token.split("|")
        if infer_genome_from_header and parts:
            genome_id = sanitize_token(parts[0], file_genome)
        else:
            genome_id = file_genome
        if len(parts) >= 3:
            contig_id = sanitize_token(parts[1], "unknown_contig")
            protein_part = "|".join(parts[2:])
            contig_inference = "header"
        elif len(parts) == 2:
            protein_part = parts[1]
            contig_id, contig_inference = infer_contig_id(protein_part)
        else:
            protein_part = token
            contig_id, contig_inference = infer_contig_id(protein_part)
    else:
        genome_id = file_genome
        protein_part = token
        contig_id, contig_inference = infer_contig_id(protein_part)
    protein_id = sanitize_token(protein_part, f"protein_{protein_index:06d}")
    return genome_id, contig_id, protein_id, contig_inference


def normalize_and_concat_proteomes(
    genomedir: str,
    out_fasta: str,
    map_path: str | None = None,
) -> dict[str, int]:
    """Normalize input proteomes into one FASTA and optional identifier map."""
    files = _list_files(genomedir)
    genome_files = assign_genome_ids(files)
    infer_genome_from_header = Path(genomedir).is_file()

    total_records = 0
    total_genomes = set()
    total_invalid = 0
    contigs_seen = set()

    with ExitStack() as stack:
        output_handle = stack.enter_context(Path(out_fasta).open("w"))
        map_handle = stack.enter_context(Path(map_path).open("w")) if map_path else None
        if map_handle:
            map_handle.write(
                "source_file\toriginal_header\tnormalized_header\tgenome_id\t"
                "contig_id\tgene_id\tcontig_inference\n"
            )

        for path, file_genome in genome_files:
            seen_proteins = set()

            for protein_index, (raw_header, raw_seq) in enumerate(
                _iter_fasta_records(path), start=1
            ):
                genome_id, contig_id, protein_id, contig_inference = _normalize_ids(
                    raw_header,
                    file_genome,
                    protein_index,
                    infer_genome_from_header=infer_genome_from_header,
                )
                base_id = protein_id
                suffix = 2
                while (contig_id, protein_id) in seen_proteins:
                    protein_id = f"{base_id}_{suffix}"
                    suffix += 1
                seen_proteins.add((contig_id, protein_id))

                seq, replaced = _clean_sequence(raw_seq)
                if not seq:
                    continue

                normalized_id = build_sequence_id(genome_id, contig_id, protein_id)
                output_handle.write(f">{normalized_id}\n{seq}\n")

                total_records += 1
                total_genomes.add(genome_id)
                contigs_seen.add((genome_id, contig_id))
                total_invalid += replaced

                if map_handle:
                    header_clean = raw_header.replace("\t", " ").strip()
                    map_handle.write(
                        f"{path}\t{header_clean}\t{normalized_id}\t{genome_id}\t"
                        f"{contig_id}\t{protein_id}\t{contig_inference}\n"
                    )

    return {
        "genomes": len(total_genomes),
        "contigs": len(contigs_seen),
        "records": total_records,
        "invalid_chars_replaced": total_invalid,
    }
