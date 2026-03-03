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
    extracted_seqs_dir, aligned_dir, modeldir, filename = args
    filepath = os.path.join(extracted_seqs_dir, filename)
    model = filename.split(".")[0]
    modfile = os.path.join(modeldir, model + ".hmm")
    sto_path = os.path.join(aligned_dir, model + ".sto")
    faa_path = os.path.join(aligned_dir, model + ".faa")

    cmd = ["hmmalign", "--trim", "-o", sto_path, modfile, filepath]
    subprocess.run(cmd, stdout=subprocess.PIPE, check=True)

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
            (cfg.extracted_seqs_dir, cfg.aligned_dir, cfg.modeldir, f) for f in files
        ])
    pool.close()
    pool.join()
