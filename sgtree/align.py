import os
import glob
import subprocess
import multiprocessing as mp
import fileinput

from Bio import AlignIO

from sgtree.config import Config


def _run_mafft(args):
    extracted_seqs_dir, aligned_dir, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    aligned_dest = os.path.join(aligned_dir, filename)
    cmd = ["mafft", "--auto", "--thread", "4", "--quiet", filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    with open(aligned_dest, "w") as f:
        f.write(result.stdout.decode("utf-8") + "\n")


def _run_mafft_linsi(args):
    extracted_seqs_dir, aligned_dir, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    aligned_dest = os.path.join(aligned_dir, filename)
    cmd = ["mafft-linsi", "--auto", "--thread", "4", "--quiet", filepath]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    with open(aligned_dest, "w") as f:
        f.write(result.stdout.decode("utf-8") + "\n")


def _run_hmmalign(args):
    extracted_seqs_dir, aligned_dir, modelset_path, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    model = filename.split(".")[0]
    sto_path = os.path.join(aligned_dir, model + ".sto")
    faa_path = os.path.join(aligned_dir, model + ".faa")

    fetch_proc = subprocess.Popen(
        ["hmmfetch", modelset_path, model],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    try:
        with fetch_proc.stdout as fetch_out:
            cmd = ["hmmalign", "--trim", "-o", sto_path, "-", filepath]
            result = subprocess.run(
                cmd,
                stdin=fetch_out,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
    finally:
        fetch_rc = fetch_proc.wait()
    if fetch_rc != 0:
        fetch_err = fetch_proc.stderr.read().decode("utf-8", errors="replace")
        raise subprocess.CalledProcessError(fetch_rc, ["hmmfetch", modelset_path, model], stderr=fetch_err)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, stderr=result.stderr)

    aln = AlignIO.read(sto_path, "stockholm")
    AlignIO.write(aln, faa_path, "fasta")

    # clean up fasta headers
    for line in fileinput.input(faa_path, inplace=True):
        line = line.rstrip()
        if not line:
            continue
        if ">" in line:
            print("|".join(line.split("|")[0:]))
        else:
            print(line)


def run_alignment(cfg: Config):
    """Run sequence alignment using the configured method (mafft, mafft-linsi, or hmmalign)."""
    os.makedirs(cfg.aligned_dir, exist_ok=True)

    files = [
        os.path.basename(f)
        for f in glob.glob(os.path.join(cfg.extracted_seqs_dir, "*"))
    ]

    print(f"- ...running {cfg.aln_method}")

    pool = mp.Pool(cfg.num_cpus)
    if cfg.aln_method == "mafft":
        pool.map(_run_mafft, [
            (cfg.extracted_seqs_dir, cfg.aligned_dir, f) for f in files
        ])
    elif cfg.aln_method == "mafft-linsi":
        pool.map(_run_mafft_linsi, [
            (cfg.extracted_seqs_dir, cfg.aligned_dir, f) for f in files
        ])
    else:
        pool.map(_run_hmmalign, [
            (cfg.extracted_seqs_dir, cfg.aligned_dir, cfg.models_path, f) for f in files
        ])
    pool.close()
    pool.join()
