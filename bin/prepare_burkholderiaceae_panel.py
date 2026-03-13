#!/usr/bin/env python
"""Prepare Burkholderiaceae benchmark selection manifests and local test data."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import duckdb
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
NELLI_ROOT = (ROOT / ".." / ".." / "nelli-genomes-db").resolve()
DB_PATH = NELLI_ROOT / "resources" / "database" / "gtdb_genomes.duckdb"
OUT_DIR = ROOT / "testgenomes"
REQUESTED_SELECTION = OUT_DIR / "burkholderiaceae_requested_50.tsv"
REQUESTED_LOOKUP = OUT_DIR / "burkholderiaceae_requested_50.lookup.tsv"
REQUESTED_LOOKUP_SIMPLE = OUT_DIR / "burkholderiaceae_requested_50.lookup"
LOCAL_DIR = OUT_DIR / "Burkholderiaceae_local"
LOCAL_LOOKUP = OUT_DIR / "burkholderiaceae_local.lookup.tsv"
LOCAL_LOOKUP_SIMPLE = OUT_DIR / "burkholderiaceae_local.lookup"

MULTI_SPECIES = [
    ("Bordetella", "Bordetella avium", 26, "strain_cluster_26"),
    ("Lautropia", "Lautropia mirabilis", 6, "strain_cluster_6"),
]
SINGLETON_GENERA = [
    "Achromobacter",
    "Acidovorax",
    "Alcaligenes",
    "Caballeronia",
    "Comamonas",
    "Cupriavidus",
    "Diaphorobacter",
    "Janthinobacterium",
    "Kinetoplastibacterium",
    "Limnobacter",
    "Limnohabitans",
    "Parasutterella",
    "Profftella",
    "Ralstonia",
    "Rhodoferax",
    "Sutterella",
    "Taylorella",
    "Zinderia",
]


def _combined_view_sql() -> str:
    return """
        with combined as (
            select
                'gtdb' as source,
                assembly_accession,
                organism_name,
                ftp_path,
                genome_size_bp,
                gtdb_representative,
                phylum,
                class,
                order_name,
                family,
                genus,
                species,
                assembly_level
            from gtdb_genomes
            union all
            select
                'ncbi' as source,
                assembly_accession,
                organism_name,
                ftp_path,
                genome_size_bp,
                'f' as gtdb_representative,
                phylum,
                class,
                order_name,
                family,
                genus,
                species,
                assembly_level
            from non_gtdb_genomes
        )
    """


def _taxonomy_string(row: pd.Series) -> str:
    return ";".join(
        [
            f"p__{row['phylum']}",
            f"c__{row['class']}",
            f"o__{row['order_name']}",
            f"f__{row['family']}",
            f"g__{row['genus']}",
            f"s__{row['species']}",
        ]
    )


def _simple_lookup_string(row: pd.Series) -> str:
    return "|".join(
        [
            str(row["phylum"]),
            str(row["class"]),
            str(row["order_name"]),
            str(row["family"]),
            str(row["genus"]),
            str(row["species"]),
        ]
    )


def _genome_id(accession: str) -> str:
    return "GAMMA__" + accession.replace("_", "-").replace(".", "-")


def _ranked_singletons(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    genera_sql = ", ".join(f"'{genus}'" for genus in SINGLETON_GENERA)
    query = f"""
        {_combined_view_sql()},
        ranked as (
            select
                *,
                row_number() over (
                    partition by genus
                    order by
                        case when gtdb_representative = 't' then 0 else 1 end,
                        case when source = 'gtdb' then 0 else 1 end,
                        genome_size_bp asc,
                        assembly_accession
                ) as rn
            from combined
            where family = 'Burkholderiaceae'
              and genus in ({genera_sql})
              and ftp_path is not null
        )
        select * from ranked where rn = 1 order by genus
    """
    frame = con.execute(query).fetchdf()
    frame["selection_role"] = "singleton"
    frame["target_cluster_size"] = 1
    return frame


def _ranked_multi_species(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for genus, species, target_size, role in MULTI_SPECIES:
        query = f"""
            {_combined_view_sql()}
            select *
            from combined
            where family = 'Burkholderiaceae'
              and genus = '{genus}'
              and species = '{species}'
              and ftp_path is not null
            order by
                case when gtdb_representative = 't' then 0 else 1 end,
                case when source = 'gtdb' then 0 else 1 end,
                genome_size_bp asc,
                assembly_accession
            limit {target_size}
        """
        frame = con.execute(query).fetchdf()
        frame["selection_role"] = role
        frame["target_cluster_size"] = target_size
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _find_local_fna(genome_id: str) -> str:
    matches = sorted(NELLI_ROOT.glob(f"data/**/fna/{genome_id}.fna"))
    return str(matches[0].resolve()) if matches else ""


def _write_lookup(frame: pd.DataFrame, tsv_path: Path, simple_path: Path) -> None:
    lookup = frame[["genome_id", "taxonomy"]].copy()
    lookup.to_csv(tsv_path, sep="\t", index=False)
    with simple_path.open("w") as handle:
        for _, row in frame.iterrows():
            handle.write(f"{row['genome_id']}\t{_simple_lookup_string(row)}\n")


def _load_local_burkholderiaceae(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    rows = []
    for path in sorted(NELLI_ROOT.glob("data/**/fna/*.fna")):
        genome_id = path.stem
        accession = genome_id.split("__", 1)[-1].replace("-", ".")
        if genome_id.startswith("GAMMA__"):
            accession = accession.replace("GCA.", "GCA_").replace("GCF.", "GCF_")
        rows.append({"genome_id": genome_id, "assembly_accession": accession, "local_fna_path": str(path.resolve())})
    local_df = pd.DataFrame(rows)
    if local_df.empty:
        return local_df
    query = f"""
        {_combined_view_sql()}
        select
            assembly_accession,
            organism_name,
            phylum,
            class,
            order_name,
            family,
            genus,
            species
        from combined
        where assembly_accession in ({", ".join(f"'{value}'" for value in local_df['assembly_accession'])})
    """
    meta = con.execute(query).fetchdf()
    merged = local_df.merge(meta, on="assembly_accession", how="left")
    merged = merged[merged["family"] == "Burkholderiaceae"].copy()
    merged = merged.sort_values(["genome_id", "local_fna_path"]).drop_duplicates(subset=["genome_id"], keep="first")
    merged["taxonomy"] = merged.apply(_taxonomy_string, axis=1)
    return merged.sort_values(["genus", "species", "genome_id"]).reset_index(drop=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with duckdb.connect(str(DB_PATH), read_only=True) as con:
        requested = pd.concat(
            [
                _ranked_multi_species(con),
                _ranked_singletons(con),
            ],
            ignore_index=True,
        )
        requested["genome_id"] = requested["assembly_accession"].map(_genome_id)
        requested["taxonomy"] = requested.apply(_taxonomy_string, axis=1)
        requested["local_fna_path"] = requested["genome_id"].map(_find_local_fna)
        requested["local_status"] = requested["local_fna_path"].apply(lambda value: "present" if value else "missing")
        requested = requested[
            [
                "selection_role",
                "target_cluster_size",
                "genome_id",
                "assembly_accession",
                "organism_name",
                "source",
                "gtdb_representative",
                "assembly_level",
                "genome_size_bp",
                "phylum",
                "class",
                "order_name",
                "family",
                "genus",
                "species",
                "taxonomy",
                "ftp_path",
                "local_status",
                "local_fna_path",
            ]
        ].sort_values(["selection_role", "genus", "species", "genome_id"]).reset_index(drop=True)
        requested.to_csv(REQUESTED_SELECTION, sep="\t", index=False)
        _write_lookup(requested, REQUESTED_LOOKUP, REQUESTED_LOOKUP_SIMPLE)

        local = _load_local_burkholderiaceae(con)
        if LOCAL_DIR.exists():
            shutil.rmtree(LOCAL_DIR)
        LOCAL_DIR.mkdir(parents=True, exist_ok=True)
        for row in local.itertuples(index=False):
            target = LOCAL_DIR / Path(str(row.local_fna_path)).name
            shutil.copy2(str(row.local_fna_path), target)
        if not local.empty:
            _write_lookup(local, LOCAL_LOOKUP, LOCAL_LOOKUP_SIMPLE)

        summary = {
            "requested_rows": len(requested),
            "requested_local_present": int((requested["local_status"] == "present").sum()),
            "local_materialized": len(local),
        }
        print(pd.Series(summary).to_string())
        print(f"requested_selection={REQUESTED_SELECTION}")
        print(f"local_dir={LOCAL_DIR}")


if __name__ == "__main__":
    main()
