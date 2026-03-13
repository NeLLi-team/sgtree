from __future__ import annotations

import glob
import os
from dataclasses import dataclass

from Bio import SeqIO

from sgtree.id_schema import build_sequence_id, sanitize_token


NUCLEOTIDE_EXTENSIONS = (".fna", ".fa", ".fasta")
PROTEIN_EXTENSIONS = (".faa",)


@dataclass(frozen=True)
class InputStageStats:
    input_format: str
    staged_source: str
    staged_records: int
    staged_genomes: int
    contigs: int


def _fasta_size_stats(path: str) -> tuple[int, int]:
    contigs = 0
    total_bases = 0
    with open(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contigs += 1
            total_bases += len(record.seq)
    return contigs, total_bases


def _list_files(input_path: str) -> list[str]:
    if os.path.isdir(input_path):
        return sorted(
            path for path in glob.glob(os.path.join(input_path, "*"))
            if os.path.isfile(path)
        )
    if os.path.isfile(input_path):
        return [input_path]
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def detect_input_format(input_path: str) -> str:
    files = _list_files(input_path)
    exts = {os.path.splitext(path)[1].lower() for path in files}
    if exts and exts <= set(PROTEIN_EXTENSIONS):
        return "faa"
    if exts and exts <= set(NUCLEOTIDE_EXTENSIONS):
        return "fna"

    first = files[0]
    sequence = ""
    with open(first) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            sequence = line.strip().upper()
            if sequence:
                break
    if sequence and set(sequence) <= set("ACGTNWSMKRYBDHV"):
        return "fna"
    return "faa"


def gene_call_inputs(
    input_path: str,
    output_dir: str,
    map_path: str,
) -> InputStageStats:
    try:
        import pyrodigal
    except ImportError as exc:
        raise RuntimeError(
            "FNA input requires pyrodigal; add it to the environment first"
        ) from exc

    os.makedirs(output_dir, exist_ok=True)
    files = _list_files(input_path)
    gene_finder = pyrodigal.GeneFinder(meta=True)

    total_records = 0
    total_contigs = 0
    genomes = set()

    with open(map_path, "w") as map_handle:
        map_handle.write(
            "source_file\tcontig_header\tnormalized_header\tgenome_id\tcontig_id\tgene_id\tbegin\tend\tstrand\ttranslation_table\n"
        )
        for file_index, path in enumerate(files, start=1):
            stem = os.path.splitext(os.path.basename(path))[0]
            genome_id = sanitize_token(stem, f"genome_{file_index:05d}")
            genomes.add(genome_id)
            out_path = os.path.join(output_dir, genome_id + ".faa")
            with open(out_path, "w") as out_handle:
                with open(path) as handle:
                    for contig_index, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
                        total_contigs += 1
                        contig_token = record.id or record.description or f"contig_{contig_index:06d}"
                        contig_id = sanitize_token(contig_token, f"contig_{contig_index:06d}")
                        genes = gene_finder.find_genes(bytes(record.seq))
                        for gene_index, gene in enumerate(genes, start=1):
                            gene_id = f"gene_{gene_index:06d}"
                            normalized_id = build_sequence_id(genome_id, contig_id, gene_id)
                            protein = str(
                                gene.translate(
                                    include_stop=False,
                                    strict=False,
                                )
                            )
                            if not protein:
                                continue
                            out_handle.write(f">{normalized_id}\n{protein}\n")
                            total_records += 1
                            map_handle.write(
                                "\t".join(
                                    [
                                        path,
                                        record.description.replace("\t", " ").strip(),
                                        normalized_id,
                                        genome_id,
                                        contig_id,
                                        gene_id,
                                        str(gene.begin),
                                        str(gene.end),
                                        str(gene.strand),
                                        str(gene.translation_table),
                                    ]
                                )
                                + "\n"
                            )

    return InputStageStats(
        input_format="fna",
        staged_source=output_dir,
        staged_records=total_records,
        staged_genomes=len(genomes),
        contigs=total_contigs,
    )


def write_genome_manifest(
    input_path: str,
    *,
    input_format: str,
    manifest_path: str,
    staged_source: str | None = None,
) -> None:
    files = _list_files(input_path)
    with open(manifest_path, "w") as handle:
        handle.write(
            "genome_id\tinput_format\tsource_file\tassembly_path\tstaged_proteome_path\tcontigs\ttotal_bases\n"
        )
        for file_index, path in enumerate(files, start=1):
            stem = os.path.splitext(os.path.basename(path))[0]
            genome_id = sanitize_token(stem, f"genome_{file_index:05d}")
            contigs = 0
            total_bases = 0
            assembly_path = ""
            staged_proteome_path = ""
            if input_format == "fna":
                contigs, total_bases = _fasta_size_stats(path)
                assembly_path = os.path.abspath(path)
                if staged_source is not None:
                    staged_proteome_path = os.path.abspath(os.path.join(staged_source, genome_id + ".faa"))
            handle.write(
                "\t".join(
                    [
                        genome_id,
                        input_format,
                        os.path.abspath(path),
                        assembly_path,
                        staged_proteome_path,
                        str(contigs),
                        str(total_bases),
                    ]
                )
                + "\n"
            )
