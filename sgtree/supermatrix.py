import os
import glob
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from sgtree.config import Config
from sgtree.parallel import map_threaded


def _trim_alignment_fallback(input_file: str, output_file: str, *, gap_threshold: float = 0.1) -> None:
    records = list(SeqIO.parse(input_file, "fasta"))
    if not records:
        Path(output_file).write_text("")
        return

    sequences = [str(record.seq) for record in records]
    width = min(len(sequence) for sequence in sequences)
    keep_columns = []
    for idx in range(width):
        nongap_fraction = sum(sequence[idx] != "-" for sequence in sequences) / len(sequences)
        if nongap_fraction >= gap_threshold:
            keep_columns.append(idx)
    if not keep_columns:
        keep_columns = list(range(width))

    with open(output_file, "w") as handle:
        for record, sequence in zip(records, sequences):
            trimmed = "".join(sequence[idx] for idx in keep_columns)
            handle.write(f">{record.id}\n{trimmed}\n")


def _run_trimal_or_fallback(input_file: str, output_file: str) -> None:
    cmd = ["trimal", "-in", input_file, "-out", output_file, "-gt", "0.1"]
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    except FileNotFoundError:
        _trim_alignment_fallback(input_file, output_file, gap_threshold=0.1)


def _run_trimal_worker(args):
    """Run trimal on a single alignment file."""
    input_file, output_file = args
    _run_trimal_or_fallback(input_file, output_file)
    # clean up fasta headers in input without using global fileinput state.
    normalized_lines: list[str] = []
    with open(input_file) as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                normalized_lines.append("|".join(line.split("|")[0:]))
            else:
                normalized_lines.append(line)
    with open(input_file, "w") as handle:
        handle.write("\n".join(normalized_lines) + ("\n" if normalized_lines else ""))


def run_trimal(cfg: Config, input_dir: str, output_dir: str):
    """Run trimal -gt 0.1 on all .faa files in input_dir, writing to output_dir."""
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(input_dir, "*.faa")))
    args = [
        (f, os.path.join(output_dir, os.path.basename(f)))
        for f in files
    ]

    map_threaded(_run_trimal_worker, args, cfg.num_cpus)


def _trimal_simple_worker(args):
    """Worker: run trimal without header cleanup."""
    input_file, output_file = args
    _run_trimal_or_fallback(input_file, output_file)


def run_trimal_simple(cfg: Config, input_dir: str, output_dir: str):
    """Run trimal without header cleanup (for protein tree trimming)."""
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(input_dir, "*.faa")))
    args = [(f, os.path.join(output_dir, os.path.basename(f))) for f in files]

    map_threaded(_trimal_simple_worker, args, cfg.num_cpus)


def build_supermatrix(trimmed_dir: str, output_dir: str, table_path: str, concat_path: str):
    """Build concatenated alignment (supermatrix) from trimmed per-marker alignments.

    Fills missing markers with 'X' gap characters.
    """
    os.makedirs(output_dir, exist_ok=True)

    marker_frames: list[pd.DataFrame] = []
    for filepath in sorted(glob.glob(os.path.join(trimmed_dir, "*.faa"))):
        marker_name = os.path.basename(filepath)
        with open(filepath) as handle:
            record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        seqs: dict[str, str] = {}
        for key, rec in record_dict.items():
            genome_id = key.split("|")[0]
            if genome_id in seqs:
                raise ValueError(
                    f"Duplicate genome id '{genome_id}' remains in alignment {marker_name}"
                )
            seqs[genome_id] = rec.format("fasta").split("\n", 1)[1]
        marker_frames.append(pd.DataFrame({marker_name: pd.Series(seqs)}))

    if not marker_frames:
        raise ValueError(f"build_supermatrix: no marker alignments found in {trimmed_dir}")

    # One outer join across all markers in a single pass (O(markers), not O(markers^2)).
    df_conc = pd.concat(marker_frames, axis=1, join="outer").sort_index(axis=1).sort_index()
    df_conc.index.name = "SeqID"
    df_conc = df_conc.reset_index()

    _fill_nan_gaps(df_conc)

    df_conc.to_csv(table_path)

    # Build the concatenated FASTA directly from the in-memory DataFrame.
    marker_cols = [c for c in df_conc.columns if c != "SeqID"]
    with open(concat_path, "w") as fp:
        for seq_id in sorted(df_conc["SeqID"]):
            row = df_conc.loc[df_conc["SeqID"] == seq_id, marker_cols].iloc[0]
            seq = "".join(str(v).replace("\n", "") for v in row.values)
            fp.write(f">{seq_id}\n{seq}\n")


def _fill_nan_gaps(df_conc: pd.DataFrame):
    """Replace NaN cells with synthetic gap strings matching each column's alignment width.

    Vectorized: for every marker column (all columns except the first), compute
    the width of each non-NaN cell, then fill NaN cells in place with ``"X" * w``
    where ``w`` is the per-row width from a forward-filled width series. Leading
    NaN rows (NaN runs that precede any non-NaN cell) are filled using the
    bottom-most non-NaN width in the column — a semantic quirk preserved from
    the legacy implementation and pinned by
    ``test_fill_nan_gaps_uses_last_non_nan_width_for_leading_gaps``.
    """
    for j in range(1, df_conc.shape[1]):
        col = df_conc.columns[j]
        series = df_conc[col]
        nan_mask = series.isna()
        if not nan_mask.any():
            continue

        known = series[~nan_mask].map(lambda v: len(str(v).replace("\n", "")))
        if known.empty:
            # Whole column is NaN — legacy backward pass fills with "X" * 0 = "".
            df_conc[col] = series.fillna("")
            continue

        widths = known.reindex(series.index).ffill().fillna(known.iloc[-1]).astype(int)
        df_conc.loc[nan_mask, col] = widths.loc[nan_mask].map(lambda n: "X" * n)
