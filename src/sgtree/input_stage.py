"""Detect input sequence formats and stage assemblies for marker searches."""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from sgtree._fasta_utils import fasta_contig_bases_stats
from sgtree.id_schema import assign_genome_ids, build_sequence_id, sanitize_token

NUCLEOTIDE_EXTENSIONS = (".fna", ".fa", ".fasta")
PROTEIN_EXTENSIONS = (".faa",)


@dataclass(frozen=True)
class InputStageStats:
    """Summarize records produced while staging an input collection."""

    input_format: str
    staged_source: str
    staged_records: int
    staged_genomes: int
    contigs: int


def _fasta_size_stats(path: str) -> tuple[int, int]:
    return fasta_contig_bases_stats(path)


def _list_files(input_path: str) -> list[str]:
    path = Path(input_path)
    if path.is_dir():
        files = sorted(
            str(candidate)
            for candidate in path.iterdir()
            if candidate.is_file() and not candidate.name.startswith(".")
        )
        if not files:
            raise ValueError(f"No input files found in directory: {input_path}")
        return files
    if path.is_file():
        return [str(path)]
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def _validate_unique_contig_ids(files: list[str]) -> None:
    """Reject assembly records that normalize to the same contig ID."""
    for path in files:
        headers_by_contig: dict[str, str] = {}
        with Path(path).open(encoding="utf-8") as handle:
            for contig_index, record in enumerate(
                SeqIO.parse(handle, "fasta"), start=1
            ):
                fallback = f"contig_{contig_index:06d}"
                contig_token = record.id or record.description or fallback
                contig_id = sanitize_token(contig_token, fallback)
                previous_header = headers_by_contig.get(contig_id)
                if previous_header is not None:
                    raise ValueError(
                        f"Duplicate normalized contig ID {contig_id!r} in {path!r}: "
                        f"{previous_header!r} and {record.description!r}"
                    )
                headers_by_contig[contig_id] = record.description


def _detect_file_format(path: str) -> str:
    with Path(path).open(encoding="utf-8") as handle:
        header = next((line.strip() for line in handle if line.strip()), "")
        if not header.startswith(">"):
            raise ValueError(f"Input file does not start with a FASTA header: {path}")
        sequence = next(
            (
                line.strip().upper()
                for line in handle
                if line.strip() and not line.startswith(">")
            ),
            "",
        )
    if not sequence:
        raise ValueError(f"Input FASTA has no sequence: {path}")
    extension = Path(path).suffix.lower()
    if extension in PROTEIN_EXTENSIONS:
        return "faa"
    if extension in NUCLEOTIDE_EXTENSIONS:
        return "fna"
    return "fna" if set(sequence) <= set("ACGTNWSMKRYBDHV") else "faa"


def detect_input_format(input_path: str) -> str:
    """Return the common FASTA format, rejecting mixed or non-FASTA inputs."""
    formats = {_detect_file_format(path) for path in _list_files(input_path)}
    if len(formats) != 1:
        raise ValueError(f"Mixed nucleotide and protein FASTA inputs: {input_path}")
    return formats.pop()


def gene_call_inputs(
    input_path: str,
    output_dir: str,
    map_path: str,
) -> InputStageStats:
    """Call genes in nucleotide FASTA files and write normalized proteins."""
    try:
        pyrodigal = importlib.import_module("pyrodigal")
    except ImportError as exc:
        raise RuntimeError(
            "FNA input requires pyrodigal; add it to the environment first"
        ) from exc

    files = _list_files(input_path)
    genome_files = assign_genome_ids(files)
    _validate_unique_contig_ids(files)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    gene_finder = pyrodigal.GeneFinder(meta=True)

    total_records = 0
    total_contigs = 0
    genomes = set()
    total_files = len(files)

    print(f"-... gene-calling {total_files} genome files with pyrodigal", flush=True)

    with Path(map_path).open("w", encoding="utf-8") as map_handle:
        map_handle.write(
            "source_file\tcontig_header\tnormalized_header\tgenome_id\tcontig_id\tgene_id\tbegin\tend\tstrand\ttranslation_table\n"
        )
        for file_index, (path, genome_id) in enumerate(genome_files, start=1):
            print(
                f"- ...gene-calling genome {file_index}/{total_files}: {genome_id}",
                flush=True,
            )
            genomes.add(genome_id)
            output_path = Path(output_dir) / f"{genome_id}.faa"
            with (
                output_path.open("w", encoding="utf-8") as output_handle,
                Path(path).open(encoding="utf-8") as input_handle,
            ):
                for contig_index, record in enumerate(
                    SeqIO.parse(input_handle, "fasta"), start=1
                ):
                    total_contigs += 1
                    contig_token = (
                        record.id or record.description or f"contig_{contig_index:06d}"
                    )
                    contig_id = sanitize_token(
                        contig_token, f"contig_{contig_index:06d}"
                    )
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
                        output_handle.write(f">{normalized_id}\n{protein}\n")
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
    """Write input provenance and assembly size statistics for each genome."""
    files = _list_files(input_path)
    genome_files = assign_genome_ids(files)
    with Path(manifest_path).open("w", encoding="utf-8") as handle:
        handle.write(
            "genome_id\tinput_format\tsource_file\tassembly_path\tstaged_proteome_path\tcontigs\ttotal_bases\n"
        )
        for path, genome_id in genome_files:
            contigs = 0
            total_bases = 0
            assembly_path = ""
            staged_proteome_path = ""
            if input_format == "fna":
                contigs, total_bases = _fasta_size_stats(path)
                assembly_path = str(Path(path).absolute())
                if staged_source is not None:
                    staged_proteome_path = str(
                        (Path(staged_source) / f"{genome_id}.faa").absolute()
                    )
            handle.write(
                "\t".join(
                    [
                        genome_id,
                        input_format,
                        str(Path(path).absolute()),
                        assembly_path,
                        staged_proteome_path,
                        str(contigs),
                        str(total_bases),
                    ]
                )
                + "\n"
            )
