import os
import glob
import shutil
import subprocess
import time

import pandas as pd
from pyhmmer import easel, hmmer, plan7

from sgtree.config import Config


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

    with open(cfg.proteomes_path, "w") as dest:
        if os.path.isdir(cfg.genomedir):
            for filename in sorted(glob.glob(os.path.join(cfg.genomedir, "*.faa"))):
                with open(filename) as src:
                    shutil.copyfileobj(src, dest)
        else:
            with open(cfg.genomedir) as src:
                shutil.copyfileobj(src, dest)

    return model_count


def run_hmmsearch(cfg: Config):
    """Run pyhmmer search with gathering cutoffs on models vs proteomes."""
    print("-... running hmmsearch")
    start = time.time()

    with plan7.HMMFile(cfg.models_path) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise ValueError(f"Marker set does not contain valid HMM entries: {cfg.models_path}")

    with easel.SequenceFile(cfg.proteomes_path, digital=True, alphabet=hmms[0].alphabet) as seq_file:
        with open(cfg.hitsoutdir, "wb") as hits_out:
            for i, hits in enumerate(
                hmmer.hmmsearch(
                    hmms,
                    seq_file,
                    cpus=max(1, cfg.num_cpus),
                    bit_cutoffs="gathering",
                )
            ):
                hits.write(hits_out, format="domains", header=(i == 0))

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
    incomplete_genomes = [
        g for g in dict_counts
        if len(dict_counts[g]) < (cfg.model_count * min_models)
    ]

    rows_to_drop = []
    for idx, row in finaldf.iterrows():
        name = row[0].split("|")[0]
        if name in incomplete_genomes:
            rows_to_drop.append(idx)
    finaldf = finaldf.drop(rows_to_drop)

    with open(os.path.join(cfg.outdir, "log_genomes_removed.txt"), "w") as f:
        f.write("\n".join(incomplete_genomes))

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
    df['namemodel'] = df[0] + "/" + df[3]
    df = df.drop_duplicates(subset='savedname', keep='first')
    df.to_csv(os.path.join(cfg.tables_dir, "before_drops_elim_incompletes"))

    # keep max 5 duplicates per namemodel
    helperdf = pd.concat(g for _, g in df.groupby("namemodel") if 5 >= len(g) >= 1)
    helperdf.to_csv(os.path.join(cfg.tables_dir, "duplicates_namemodel"))

    newdf = df.drop_duplicates(subset='namemodel', keep=False)
    newdf.to_csv(os.path.join(cfg.tables_dir, "dropped_namemodel"))

    df = pd.merge(helperdf, newdf, how='outer')
    df.to_csv(os.path.join(cfg.tables_dir, "merged_final"))

    # merge with reference data if available
    if cfg.ref is not None:
        ref_path = os.path.join(cfg.ref_dir_path(), "tables", "merged_final")
        df_ref = pd.read_csv(ref_path)
        df = pd.concat([df, df_ref])

    df_fordups = df.set_index(df["savedname"])
    df_fordups.to_csv(os.path.join(cfg.outdir, "table_elim_dups"))

    return df, df_fordups
