import os
import glob
import subprocess
import multiprocessing as mp
import fileinput

import pandas as pd
from Bio import SeqIO

from sgtree.config import Config


def _map_with_fallback(func, args, workers: int):
    if not args:
        return
    n_workers = max(1, min(workers, len(args)))
    if n_workers == 1:
        for item in args:
            func(item)
        return
    try:
        with mp.Pool(n_workers) as pool:
            pool.map(func, args)
    except (PermissionError, OSError) as e:
        print(f"warning: multiprocessing unavailable ({e}); falling back to serial execution")
        for item in args:
            func(item)


def _run_trimal_worker(args):
    """Run trimal on a single alignment file."""
    input_file, output_file = args
    cmd = ["trimal", "-in", input_file, "-out", output_file, "-gt", "0.1"]
    subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    # clean up fasta headers in input
    for line in fileinput.input(input_file, inplace=True):
        line = line.rstrip()
        if not line:
            continue
        if ">" in line:
            print("|".join(line.split("|")[0:]))
        else:
            print(line)


def run_trimal(cfg: Config, input_dir: str, output_dir: str):
    """Run trimal -gt 0.1 on all .faa files in input_dir, writing to output_dir."""
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(input_dir, "*.faa")))
    args = [
        (f, os.path.join(output_dir, os.path.basename(f)))
        for f in files
    ]

    _map_with_fallback(_run_trimal_worker, args, cfg.num_cpus)


def _trimal_simple_worker(args):
    """Worker: run trimal without header cleanup."""
    input_file, output_file = args
    cmd = ["trimal", "-in", input_file, "-out", output_file, "-gt", "0.1"]
    subprocess.run(cmd, stdout=subprocess.PIPE, check=True)


def run_trimal_simple(cfg: Config, input_dir: str, output_dir: str):
    """Run trimal without header cleanup (for protein tree trimming)."""
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(input_dir, "*.faa")))
    args = [(f, os.path.join(output_dir, os.path.basename(f))) for f in files]

    _map_with_fallback(_trimal_simple_worker, args, cfg.num_cpus)


def build_supermatrix(trimmed_dir: str, output_dir: str, table_path: str, concat_path: str):
    """Build concatenated alignment (supermatrix) from trimmed per-marker alignments.

    Fills missing markers with 'X' gap characters.
    """
    os.makedirs(output_dir, exist_ok=True)

    # build dataframe with all trimmed alignments
    df_conc = pd.DataFrame(columns=["SeqID"])
    for filepath in sorted(glob.glob(os.path.join(trimmed_dir, "*.faa"))):
        with open(filepath) as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        record_dict = {k: v.format("fasta").split("\n", 1)[1] for k, v in record_dict.items()}
        new_dict = {}
        for key in record_dict:
            genome_id = key.split("|")[0]
            if genome_id in new_dict:
                raise ValueError(
                    f"Duplicate genome id '{genome_id}' remains in alignment {os.path.basename(filepath)}"
                )
            new_dict[genome_id] = record_dict[key]
        new_df = pd.DataFrame(
            list(new_dict.items()),
            columns=["SeqID", os.path.basename(filepath)],
        )
        df_conc = pd.merge(new_df, df_conc, how="outer")

    marker_cols = sorted(col for col in df_conc.columns if col != "SeqID")
    df_conc = df_conc[["SeqID"] + marker_cols].sort_values("SeqID")

    # fill NaN cells with X characters of appropriate length
    _fill_nan_gaps(df_conc)

    # save intermediate table
    df_conc.to_csv(table_path)

    # rebuild from saved CSV and write concatenated FASTA
    df_conc = pd.read_csv(table_path)
    df_conc = df_conc.set_index("SeqID").sort_index()
    record_dict = df_conc.T.to_dict("list")
    record_dict = {k: v[1:] for k, v in record_dict.items()}
    record_dict = {k: "".join(str(x) for x in v) for k, v in record_dict.items()}
    record_dict = {k: v.replace("\n", "") for k, v in record_dict.items()}

    with open(concat_path, "w") as f:
        for k in sorted(record_dict):
            v = record_dict[k]
            f.write(f">{k}\n{v}\n")


def _fill_nan_gaps(df_conc: pd.DataFrame):
    """Replace NaN cells with X characters matching the alignment width of that column."""
    # forward pass: fill NaN cells using known lengths from previous rows
    for j in range(1, df_conc.shape[1]):
        len_string = 0
        for i in range(df_conc.shape[0]):
            val = df_conc.iloc[i, j]
            if isinstance(val, float):  # NaN
                if len_string > 0:
                    df_conc.iloc[i, j] = "X" * len_string
            else:
                len_string = len(str(val).replace("\n", ""))

    # backward pass: fill leading NaN rows (where first rows had no data)
    for j in range(1, df_conc.shape[1]):
        if isinstance(df_conc.iloc[0, j], float):
            # find reference length from last non-NaN value
            ref_len = 0
            for k in range(df_conc.shape[0] - 1, -1, -1):
                val = df_conc.iloc[k, j]
                if not isinstance(val, float):
                    ref_len = len(str(val).replace("\n", ""))
                    break
            for i in range(df_conc.shape[0]):
                if isinstance(df_conc.iloc[i, j], float):
                    df_conc.iloc[i, j] = "X" * ref_len
                else:
                    break
