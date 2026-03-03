import os
import glob
import shutil
import subprocess
import time

import pandas as pd
from pyhmmer import easel, hmmer, plan7

from sgtree.config import Config
from sgtree.fasta_normalize import normalize_and_concat_proteomes


def _cap_namemodel_duplicates(df: pd.DataFrame, max_per_group: int = 5) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    if "namemodel" not in df.columns:
        raise ValueError("Missing required column 'namemodel' for duplicate capping")
    if "score_bits" not in df.columns:
        raise ValueError("Missing required column 'score_bits' for duplicate capping")
    return (
        df.sort_values(["namemodel", "score_bits", "savedname"], ascending=[True, False, True])
        .groupby("namemodel", group_keys=False)
        .head(max_per_group)
        .copy()
    )


def _count_models_in_hmm(models_path: str) -> int:
    count = 0
    with open(models_path, "rb") as handle:
        for line in handle:
            if line.startswith(b"NAME"):
                count += 1
    return count


def concat_inputs(cfg: Config) -> int:
    """Stage HMM models and concatenate genome FASTA files into single files.

    Returns the number of models.
    """
    if os.path.isdir(cfg.modeldir):
        with open(cfg.models_path, "w") as dest:
            for filename in sorted(glob.glob(os.path.join(cfg.modeldir, "*.hmm"))):
                with open(filename) as src:
                    shutil.copyfileobj(src, dest)
    else:
        shutil.copyfile(cfg.modeldir, cfg.models_path)

    model_count = _count_models_in_hmm(cfg.models_path)
    if model_count == 0:
        raise ValueError(f"No HMM models found in marker set: {cfg.modeldir}")

    os.makedirs(cfg.outdir, exist_ok=True)
    map_path = os.path.join(cfg.outdir, "proteomes_header_map.tsv")
    stats = normalize_and_concat_proteomes(cfg.genomedir, cfg.proteomes_path, map_path)
    print(
        "-... normalized proteomes "
        f"(genomes={stats['genomes']}, records={stats['records']}, "
        f"invalid_chars_replaced={stats['invalid_chars_replaced']})"
    )

    return model_count


def run_hmmsearch(cfg: Config):
    """Run pyhmmer search on models vs proteomes with configurable threshold mode."""
    print("-... running hmmsearch")
    start = time.time()

    with plan7.HMMFile(cfg.models_path) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise ValueError(f"Marker set does not contain valid HMM entries: {cfg.models_path}")

    base_opts = {}
    if cfg.hmmsearch_cutoff == "cut_ga":
        base_opts["bit_cutoffs"] = "gathering"
    elif cfg.hmmsearch_cutoff == "cut_tc":
        base_opts["bit_cutoffs"] = "trusted"
    elif cfg.hmmsearch_cutoff == "cut_nc":
        base_opts["bit_cutoffs"] = "noise"
    else:
        base_opts["E"] = cfg.hmmsearch_evalue
        base_opts["domE"] = cfg.hmmsearch_evalue

    requested_cpus = max(1, cfg.num_cpus)

    def _run_search(cpus: int):
        search_opts = dict(base_opts)
        search_opts["cpus"] = cpus
        with easel.SequenceFile(cfg.proteomes_path, digital=True, alphabet=hmms[0].alphabet) as seq_file:
            with open(cfg.hitsoutdir, "wb") as hits_out:
                for i, hits in enumerate(hmmer.hmmsearch(hmms, seq_file, **search_opts)):
                    hits.write(hits_out, format="domains", header=(i == 0))

    try:
        _run_search(requested_cpus)
    except PermissionError as e:
        if requested_cpus == 1:
            raise
        print(
            f"warning: pyhmmer multiprocessing unavailable ({e}); "
            "retrying hmmsearch with cpus=1"
        )
        _run_search(1)

    elapsed = time.time() - start
    print(f"\nmarker protein detection done - runtime: {elapsed:.1f} seconds")
    print("=" * 80)
    return elapsed


def parse_hmmsearch(cfg: Config) -> tuple[pd.DataFrame, dict]:
    """Parse hmmsearch domtblout, count markers per genome, filter incomplete genomes.

    Returns (finaldf, dict_counts).
    """
    os.makedirs(cfg.tables_dir, exist_ok=True)
    min_models = cfg.min_models_fraction
    dict_counts = {}

    print("\n - ...extracting best hits")
    print("MINIMUM MODELS", round(min_models * cfg.genome_count))

    # optional length filtering
    if cfg.lflt_fraction > 0:
        df = pd.read_csv(cfg.hitsoutdir, comment="#", sep=r'\s+', header=None)
        filtered_out = df[df[2] < df.groupby(3)[2].transform('median') * cfg.lflt_fraction][0]
        list_torm_path = cfg.hitsoutdir + ".del.ls"
        filtered_out.to_csv(list_torm_path, index=False, header=False)
        out_lengthfilter = cfg.hitsoutdir + ".lfilt"
        with open(out_lengthfilter, "w") as l_out:
            subprocess.run(["grep", "-vwFf", list_torm_path, cfg.hitsoutdir], stdout=l_out)
        shutil.move(out_lengthfilter, cfg.hitsoutdir)

    # read domtblout as whitespace-separated, all strings
    finaldf = pd.read_csv(cfg.hitsoutdir, comment="#", sep=r'\s+', header=None, dtype=str)

    # count markers per genome
    seen = set()
    for _, row in finaldf.iterrows():
        name = row[0].split("|")[0]
        model = row[3]
        proteinid = row[0]
        key = f"{proteinid}_{model}"

        if key not in seen:
            seen.add(key)
            if name not in dict_counts:
                dict_counts[name] = {}
            if model in dict_counts[name]:
                dict_counts[name][model] += 1
            else:
                dict_counts[name][model] = 1

    # filter incomplete genomes
    incomplete_genomes = {
        g for g in dict_counts
        if len(dict_counts[g]) < (cfg.model_count * min_models)
    }

    removed_reasons = {g: [f"minmarker:{len(dict_counts[g]) / cfg.model_count:.4f}"] for g in incomplete_genomes}
    removed_genomes = set(incomplete_genomes)

    if cfg.max_sdup >= 0:
        high_single_dup = {
            g for g, counts in dict_counts.items()
            if counts and max(counts.values()) > cfg.max_sdup
        }
        for g in high_single_dup:
            removed_reasons.setdefault(g, []).append(f"maxsdup:{max(dict_counts[g].values())}")
        removed_genomes.update(high_single_dup)

    if cfg.max_dupl >= 0:
        high_dup_fraction = set()
        for g, counts in dict_counts.items():
            dup_fraction = sum(1 for v in counts.values() if v > 1) / cfg.model_count
            if dup_fraction > cfg.max_dupl:
                high_dup_fraction.add(g)
                removed_reasons.setdefault(g, []).append(f"maxdupl:{dup_fraction:.4f}")
        removed_genomes.update(high_dup_fraction)

    rows_to_drop = []
    for idx, row in finaldf.iterrows():
        name = row[0].split("|")[0]
        if name in removed_genomes:
            rows_to_drop.append(idx)
    finaldf = finaldf.drop(rows_to_drop)

    with open(os.path.join(cfg.outdir, "log_genomes_removed.txt"), "w") as f:
        for genome in sorted(removed_genomes):
            reasons = ";".join(removed_reasons.get(genome, ["filtered"]))
            f.write(f"{genome}\t{reasons}\n")

    print(f"AFTER hmmout {finaldf.shape} # rows deleted {len(rows_to_drop)}")

    # write marker count matrix
    count_mat = pd.DataFrame.from_dict(dict_counts).fillna(0)
    count_mat.to_csv(os.path.join(cfg.outdir, "marker_count_matrix.csv"))

    return finaldf, dict_counts


def build_working_df(cfg: Config, finaldf: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the working dataframe with savedname and namemodel columns.

    Handles duplicate filtering and reference merging.
    Returns (df, df_fordups).
    """
    finaldf = finaldf.copy()
    finaldf["savedname"] = finaldf[0].apply(lambda c: c.replace("|", "/"))
    df = finaldf.map(lambda x: x.split("|")[0] if isinstance(x, str) else x)
    score_col = 7 if 7 in df.columns else ("7" if "7" in df.columns else None)
    if score_col is None:
        raise ValueError("Expected HMMER bitscore column '7' in parsed domtblout table")
    df["score_bits"] = pd.to_numeric(df[score_col], errors="coerce")
    if df["score_bits"].isna().any():
        raise ValueError("Failed to parse one or more HMMER bitscores from column '7'")
    df['namemodel'] = df[0] + "/" + df[3]
    df = df.drop_duplicates(subset='savedname', keep='first')
    df.to_csv(os.path.join(cfg.tables_dir, "before_drops_elim_incompletes"))

    # cap duplicate copy count per genome+marker group (keep best-scoring copies)
    capped = _cap_namemodel_duplicates(df, max_per_group=5)
    capped.to_csv(os.path.join(cfg.tables_dir, "duplicates_namemodel"))

    dropped = df.loc[~df.index.isin(capped.index)]
    dropped.to_csv(os.path.join(cfg.tables_dir, "dropped_namemodel"))

    df = capped
    df.to_csv(os.path.join(cfg.tables_dir, "merged_final"))

    # merge with reference data if available
    if cfg.ref is not None:
        ref_path = os.path.join(cfg.ref_dir_path(), "tables", "merged_final")
        df_ref = pd.read_csv(ref_path)
        df = pd.concat([df, df_ref])

    df_fordups = df.set_index(df["savedname"])
    df_fordups.to_csv(os.path.join(cfg.outdir, "table_elim_dups"))

    return df, df_fordups
