from __future__ import annotations

import gzip
import shutil
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


DEFAULT_BURKHOLDERIACEAE_SPECIES_PLAN: tuple[dict[str, object], ...] = (
    {"genus": "Burkholderia", "species": "Burkholderia thailandensis", "target_count": 1, "role": "singleton"},
    {"genus": "Ralstonia", "species": "Ralstonia nicotianae", "target_count": 1, "role": "singleton"},
    {"genus": "Comamonas", "species": "Comamonas acidovorans", "target_count": 1, "role": "singleton"},
    {"genus": "Taylorella", "species": "Taylorella equigenitalis", "target_count": 1, "role": "singleton"},
    {"genus": "Sutterella", "species": "Sutterella wadsworthensis", "target_count": 1, "role": "singleton"},
    {"genus": "Diaphorobacter", "species": "Diaphorobacter nitroreducens", "target_count": 1, "role": "singleton"},
    {"genus": "Alcaligenes", "species": "Alcaligenes phenolicus", "target_count": 1, "role": "singleton"},
    {"genus": "Parasutterella", "species": "Parasutterella excrementihominis", "target_count": 1, "role": "singleton"},
    {"genus": "Caballeronia", "species": "Caballeronia zhejiangensis", "target_count": 1, "role": "singleton"},
    {"genus": "Paracidovorax", "species": "Paracidovorax avenae", "target_count": 1, "role": "singleton"},
    {"genus": "Cupriavidus", "species": "Cupriavidus metallidurans", "target_count": 1, "role": "singleton"},
    {"genus": "Pandoraea", "species": "Pandoraea apista", "target_count": 1, "role": "singleton"},
    {"genus": "Kerstersia", "species": "Kerstersia gyiorum", "target_count": 1, "role": "singleton"},
    {"genus": "Acidovorax", "species": "Acidovorax facilis", "target_count": 1, "role": "singleton"},
    {"genus": "Herbaspirillum", "species": "Herbaspirillum huttiense", "target_count": 1, "role": "singleton"},
    {"genus": "Janthinobacterium", "species": "Janthinobacterium lividum", "target_count": 1, "role": "singleton"},
    {"genus": "Paraburkholderia", "species": "Paraburkholderia fungorum", "target_count": 1, "role": "singleton"},
    {"genus": "Thiomonas", "species": "Thiomonas arsenitoxydans", "target_count": 1, "role": "singleton"},
    {"genus": "Achromobacter", "species": "Achromobacter xylosoxidans", "target_count": 6, "role": "strain_set_6"},
    {"genus": "Bordetella", "species": "Bordetella pertussis", "target_count": 26, "role": "strain_set_26"},
)


def _species_plan_frame(
    species_plan: tuple[dict[str, object], ...] | list[dict[str, object]] | None = None,
) -> pd.DataFrame:
    frame = pd.DataFrame(species_plan or DEFAULT_BURKHOLDERIACEAE_SPECIES_PLAN).copy()
    required = {"genus", "species", "target_count", "role"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Species plan is missing required columns: {', '.join(sorted(missing))}")
    frame["target_count"] = frame["target_count"].astype(int)
    frame["plan_order"] = range(1, len(frame) + 1)
    return frame


def select_species_rows(
    candidates: pd.DataFrame,
    species_plan: tuple[dict[str, object], ...] | list[dict[str, object]] | None = None,
    *,
    prefix: str = "BURK__",
) -> pd.DataFrame:
    """Pick benchmark genomes from a candidate table using stable quality-aware ranking."""
    required = {
        "assembly_accession",
        "organism_name",
        "ftp_path",
        "genome_size_bp",
        "assembly_level",
        "refseq_category",
        "phylum",
        "class",
        "order_name",
        "family",
        "genus",
        "species",
        "taxonomy_source",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table is missing required columns: {', '.join(sorted(missing))}")

    plan = _species_plan_frame(species_plan)
    merged = candidates.merge(plan, on=["genus", "species"], how="inner", validate="many_to_one").copy()
    if merged.empty:
        raise ValueError("No candidate genomes matched the requested Burkholderiaceae species plan")

    merged["assembly_rank"] = merged["assembly_level"].map(
        {
            "Complete Genome": 0,
            "Chromosome": 1,
            "Scaffold": 2,
            "Contig": 3,
        }
    ).fillna(4).astype(int)
    merged["refseq_rank"] = merged["refseq_category"].map(
        {
            "reference genome": 0,
            "representative genome": 1,
        }
    ).fillna(2).astype(int)
    merged["size_delta_bp"] = (
        merged["genome_size_bp"]
        - merged.groupby(["genus", "species"])["genome_size_bp"].transform("median")
    ).abs()
    merged["pick_rank"] = (
        merged.sort_values(
            [
                "plan_order",
                "genus",
                "species",
                "assembly_rank",
                "refseq_rank",
                "size_delta_bp",
                "assembly_accession",
            ],
            ascending=[True, True, True, True, True, True, True],
        )
        .groupby(["genus", "species"], group_keys=False)
        .cumcount()
        + 1
    )

    selected = merged[merged["pick_rank"] <= merged["target_count"]].copy()
    counts = selected.groupby(["genus", "species"]).size()
    missing_species = []
    for row in plan.itertuples(index=False):
        picked = int(counts.get((row.genus, row.species), 0))
        if picked != int(row.target_count):
            missing_species.append(f"{row.genus} {row.species}: expected {row.target_count}, found {picked}")
    if missing_species:
        raise ValueError("Could not satisfy benchmark sampling plan: " + "; ".join(missing_species))

    selected = selected.sort_values(["plan_order", "pick_rank", "assembly_accession"]).reset_index(drop=True)
    selected["accession_mod"] = (
        selected["assembly_accession"].astype(str).str.replace("_", "-", regex=False).str.replace(".", "-", regex=False)
    )
    selected["genome_id"] = prefix + selected["accession_mod"]
    selected["taxonomy_lookup"] = selected[
        ["phylum", "class", "order_name", "family", "genus", "species"]
    ].fillna("unclassified").astype(str).agg("|".join, axis=1)
    selected["taxonomy_string"] = selected.apply(
        lambda row: ";".join(
            [
                f"p__{row['phylum']}",
                f"c__{row['class']}",
                f"o__{row['order_name']}",
                f"f__{row['family']}",
                f"g__{row['genus']}",
                f"s__{row['species']}",
            ]
        ),
        axis=1,
    )
    return selected


def build_burkholderiaceae_selection(
    taxonomy_db_path: Path,
    *,
    species_plan: tuple[dict[str, object], ...] | list[dict[str, object]] | None = None,
    prefix: str = "BURK__",
) -> pd.DataFrame:
    try:
        import duckdb
    except ImportError as exc:
        raise RuntimeError("Burkholderiaceae benchmark preparation requires duckdb in the SGTree environment") from exc

    plan = _species_plan_frame(species_plan)
    if not taxonomy_db_path.exists():
        raise FileNotFoundError(taxonomy_db_path)

    with duckdb.connect(str(taxonomy_db_path), read_only=True) as con:
        con.register("species_plan", plan[["genus", "species"]])
        candidates = con.execute(
            """
            with combined as (
                select
                    assembly_accession,
                    organism_name,
                    ftp_path,
                    genome_size_bp,
                    assembly_level,
                    refseq_category,
                    phylum,
                    class,
                    order_name,
                    family,
                    genus,
                    species,
                    'gtdb' as taxonomy_source
                from gtdb_genomes
                union all
                select
                    assembly_accession,
                    organism_name,
                    ftp_path,
                    genome_size_bp,
                    assembly_level,
                    refseq_category,
                    phylum,
                    class,
                    order_name,
                    family,
                    genus,
                    species,
                    'non_gtdb' as taxonomy_source
                from non_gtdb_genomes
            ),
            deduplicated as (
                select
                    *,
                    row_number() over (
                        partition by assembly_accession
                        order by case when taxonomy_source = 'gtdb' then 0 else 1 end
                    ) as rn
                from combined
            )
            select
                assembly_accession,
                organism_name,
                ftp_path,
                genome_size_bp,
                assembly_level,
                refseq_category,
                phylum,
                class,
                order_name,
                family,
                genus,
                species,
                taxonomy_source
            from deduplicated
            inner join species_plan using (genus, species)
            where rn = 1
              and ftp_path is not null
              and ftp_path <> ''
            """
        ).fetchdf()
    return select_species_rows(candidates, species_plan=species_plan, prefix=prefix)


def _download_and_rewrite_fasta(url: str, output_path: Path, genome_id: str, *, timeout: int, retries: int) -> None:
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp.gz")
    for attempt in range(1, retries + 1):
        try:
            try:
                response = urllib.request.urlopen(url, timeout=timeout)
            except urllib.error.URLError as exc:
                reason = getattr(exc, "reason", None)
                if isinstance(reason, ssl.SSLCertVerificationError):
                    insecure_context = ssl._create_unverified_context()
                    response = urllib.request.urlopen(url, timeout=timeout, context=insecure_context)
                else:
                    raise
            with response, open(temp_path, "wb") as handle:
                shutil.copyfileobj(response, handle)
            with gzip.open(temp_path, "rt") as src, open(output_path, "w") as dst:
                for line in src:
                    if line.startswith(">"):
                        header = line[1:].strip().split()[0]
                        dst.write(f">{genome_id}|{header}\n")
                    else:
                        dst.write(line)
            temp_path.unlink(missing_ok=True)
            return
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            temp_path.unlink(missing_ok=True)
            if attempt == retries:
                raise RuntimeError(f"Failed to download {url}: {exc}") from exc
            time.sleep(float(attempt))


def download_benchmark_fna(
    selection: pd.DataFrame,
    output_dir: Path,
    *,
    overwrite: bool = False,
    timeout: int = 300,
    retries: int = 3,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for row in selection.itertuples(index=False):
        output_path = output_dir / f"{row.genome_id}.fna"
        if output_path.exists() and not overwrite:
            continue
        base_name = str(row.ftp_path).rstrip("/").split("/")[-1]
        genome_url = f"{row.ftp_path}/{base_name}_genomic.fna.gz"
        _download_and_rewrite_fasta(
            genome_url,
            output_path,
            row.genome_id,
            timeout=timeout,
            retries=retries,
        )


def write_benchmark_metadata(
    selection: pd.DataFrame,
    *,
    lookup_path: Path,
    taxonomy_tsv_path: Path,
    selection_tsv_path: Path,
) -> None:
    lookup_path.parent.mkdir(parents=True, exist_ok=True)
    taxonomy_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    selection_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    selection[["genome_id", "taxonomy_lookup"]].to_csv(
        lookup_path,
        sep="\t",
        index=False,
        header=False,
    )
    selection[
        [
            "genome_id",
            "assembly_accession",
            "organism_name",
            "taxonomy_source",
            "taxonomy_string",
            "phylum",
            "class",
            "order_name",
            "family",
            "genus",
            "species",
        ]
    ].to_csv(
        taxonomy_tsv_path,
        sep="\t",
        index=False,
    )
    selection[
        [
            "plan_order",
            "role",
            "target_count",
            "pick_rank",
            "genome_id",
            "assembly_accession",
            "assembly_level",
            "refseq_category",
            "genome_size_bp",
            "organism_name",
            "phylum",
            "class",
            "order_name",
            "family",
            "genus",
            "species",
            "ftp_path",
        ]
    ].to_csv(
        selection_tsv_path,
        sep="\t",
        index=False,
    )


def prepare_burkholderiaceae_benchmark_dataset(
    taxonomy_db_path: Path,
    output_dir: Path,
    *,
    lookup_path: Path,
    taxonomy_tsv_path: Path,
    selection_tsv_path: Path,
    prefix: str = "BURK__",
    overwrite: bool = False,
    timeout: int = 300,
    retries: int = 3,
) -> pd.DataFrame:
    selection = build_burkholderiaceae_selection(taxonomy_db_path, prefix=prefix)
    download_benchmark_fna(
        selection,
        output_dir,
        overwrite=overwrite,
        timeout=timeout,
        retries=retries,
    )
    write_benchmark_metadata(
        selection,
        lookup_path=lookup_path,
        taxonomy_tsv_path=taxonomy_tsv_path,
        selection_tsv_path=selection_tsv_path,
    )
    return selection
