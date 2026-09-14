"""Align extracted marker sequences with external or HMM-based aligners."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from pyhmmer import easel, hmmer, plan7

from sgtree._subprocess import run_capture
from sgtree.config import Config
from sgtree.parallel import map_processed, map_threaded


def _clean_subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    mafft_prefix = Path.cwd() / ".pixi" / "envs" / "default" / "libexec" / "mafft"
    if mafft_prefix.is_dir():
        env["MAFFT_BINARIES"] = str(mafft_prefix)
    else:
        env.pop("MAFFT_BINARIES", None)
    return env


def _alignment_taxa_count(filepath: str) -> int:
    count = 0
    with Path(filepath).open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def _run_mafft(args: tuple[str, str, str, str, int]) -> None:
    binary, extracted_seqs_dir, aligned_dir, filename, threads = args
    filepath = str(Path(extracted_seqs_dir) / filename)
    aligned_dest = Path(aligned_dir) / filename
    cmd = [binary, "--auto", "--thread", str(threads), "--quiet", filepath]
    result = run_capture(cmd, env=_clean_subprocess_env())
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )
    with aligned_dest.open("w") as handle:
        handle.write(result.stdout + "\n")


def _run_famsa(args: tuple[str, str, str, int]) -> None:
    extracted_seqs_dir, aligned_dir, filename, threads = args
    filepath = str(Path(extracted_seqs_dir) / filename)
    aligned_dest = str(Path(aligned_dir) / filename)
    cmd = [
        "famsa",
        "-t",
        str(threads),
        "-refine_mode",
        "on",
        filepath,
        aligned_dest,
    ]
    result = run_capture(cmd, env=_clean_subprocess_env())
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr
        )


def _run_hmmalign(args: tuple[str, str, str, str]) -> None:
    extracted_seqs_dir, aligned_dir, model_hmm_path, filename = args
    filepath = str(Path(extracted_seqs_dir) / filename)
    model = Path(filename).stem
    faa_path = Path(aligned_dir) / f"{model}.faa"
    with plan7.HMMFile(model_hmm_path) as hmm_file:
        hmm_profile = next(iter(hmm_file), None)
    if hmm_profile is None:
        raise ValueError(f"Could not read HMM profile from {model_hmm_path}")

    with easel.SequenceFile(
        filepath, digital=True, alphabet=hmm_profile.alphabet
    ) as seq_file:
        msa = hmmer.hmmalign(hmm_profile, seq_file, cpus=1, trim=True)
    with faa_path.open("wb") as out_handle:
        msa.write(out_handle, format="afa")
    _normalize_fasta(str(faa_path))


def _normalize_fasta(fasta_path: str) -> None:
    """Normalize FASTA headers and sequence order for deterministic downstream trees."""
    records = []
    header = None
    seq_chunks = []
    with Path(fasta_path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = "|".join(line.split("|")[0:])
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append((header, "".join(seq_chunks)))

    records.sort(key=lambda item: item[0])
    with Path(fasta_path).open("w") as handle:
        for header, seq in records:
            handle.write(f"{header}\n{seq}\n")


def _split_models(models_path: str, split_dir: str) -> None:
    split_path = Path(split_dir)
    if split_path.is_dir():
        shutil.rmtree(split_path)
    split_path.mkdir(parents=True, exist_ok=True)
    with plan7.HMMFile(models_path) as hmm_file:
        for hmm_profile in hmm_file:
            marker = hmm_profile.name
            if isinstance(marker, bytes):
                marker = marker.decode("utf-8", errors="replace")
            out_path = split_path / f"{marker}.hmm"
            with out_path.open("wb") as out_handle:
                hmm_profile.write(out_handle)


def run_alignment(
    cfg: Config,
    *,
    extracted_seqs_dir: str | None = None,
    aligned_dir: str | None = None,
) -> None:
    """Run sequence alignment using the configured method.

    By default this uses the standard extracted/aligned directories from the
    config, but marker-selection cleanup can provide alternate directories when
    alignments need to be rebuilt from cleaned sequence sets.
    """
    extracted_seqs_dir = extracted_seqs_dir or cfg.extracted_seqs_dir
    aligned_dir = aligned_dir or cfg.aligned_dir
    Path(aligned_dir).mkdir(parents=True, exist_ok=True)

    files = [path.name for path in Path(extracted_seqs_dir).glob("*")]

    print(f"- ...running {cfg.aln_method}")

    if cfg.aln_method in ("mafft", "mafft-linsi", "famsa"):
        if cfg.aln_method == "famsa":
            famsa_args = [
                (extracted_seqs_dir, aligned_dir, filename, 1) for filename in files
            ]
            map_threaded(_run_famsa, famsa_args, cfg.num_cpus)
            return

        small: list[tuple[str, str, str, str, int]] = []
        large: list[tuple[str, str, str, str, int]] = []
        large_threads = max(1, min(4, cfg.num_cpus))
        for filename in files:
            filepath = str(Path(extracted_seqs_dir) / filename)
            taxa = _alignment_taxa_count(filepath)
            if taxa < 100:
                small.append(
                    (
                        cfg.aln_method,
                        extracted_seqs_dir,
                        aligned_dir,
                        filename,
                        1,
                    )
                )
            else:
                large.append(
                    (
                        cfg.aln_method,
                        extracted_seqs_dir,
                        aligned_dir,
                        filename,
                        large_threads,
                    )
                )
        map_threaded(
            _run_mafft,
            large,
            max(1, min(len(large), cfg.num_cpus // large_threads if large else 1)),
        )
        map_threaded(_run_mafft, small, cfg.num_cpus)
    else:
        split_dir = Path(cfg.outdir) / "models_split"
        _split_models(cfg.models_path, str(split_dir))
        n_jobs = max(1, min(cfg.num_cpus, len(files) if files else 1))
        map_processed(
            _run_hmmalign,
            [
                (
                    extracted_seqs_dir,
                    aligned_dir,
                    str(Path(split_dir) / f"{Path(f).stem}.hmm"),
                    f,
                )
                for f in files
            ],
            n_jobs,
        )
