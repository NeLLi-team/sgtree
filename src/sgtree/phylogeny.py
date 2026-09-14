"""Build species and per-marker phylogenetic trees."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import Phylo, SeqIO

from sgtree._subprocess import run_capture, run_check
from sgtree.config import Config
from sgtree.parallel import map_threaded

_NEWICK_SAFE_ID = re.compile(r"^[A-Za-z0-9_.|+\-]+$")


def _fasttree_executable() -> str:
    for candidate in ("VeryFastTree", "veryfasttree"):
        if shutil.which(candidate):
            return candidate
    raise FileNotFoundError("Could not find VeryFastTree in PATH")


# Per-executable cache of whether the binary accepts the `-threads` flag.
# Populated lazily from the first real run; subsequent calls skip the probe.
_THREADS_SUPPORT: dict[str, bool] = {}


def _run_fasttree_with_optional_threads(
    threaded_cmd: list[str], fallback_cmd: list[str]
) -> None:
    executable = threaded_cmd[0]
    if _THREADS_SUPPORT.get(executable) is False:
        run_check(fallback_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return

    threaded = run_capture(threaded_cmd)
    if threaded.returncode == 0:
        _THREADS_SUPPORT[executable] = True
        return

    message = (threaded.stderr or "") + "\n" + (threaded.stdout or "")
    if "-threads" not in message:
        raise subprocess.CalledProcessError(
            returncode=threaded.returncode,
            cmd=threaded_cmd,
            output=threaded.stdout,
            stderr=threaded.stderr,
        )
    _THREADS_SUPPORT[executable] = False
    run_check(fallback_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_fasttree(input_fasta: str, output_tree: str, threads: int) -> None:
    """Run VeryFastTree to build an approximate ML tree."""
    executable = _fasttree_executable()
    threaded_cmd = [
        executable,
        "-threads",
        str(max(1, threads)),
        "-quiet",
        "-out",
        output_tree,
        input_fasta,
    ]
    fallback_cmd = [executable, "-quiet", "-out", output_tree, input_fasta]
    _run_fasttree_with_optional_threads(threaded_cmd, fallback_cmd)


def _run_iqtree(
    input_fasta: str,
    output_tree: str,
    cpus: int,
    model: str,
    *,
    fast: bool,
) -> None:
    """Run IQ-TREE and copy resulting treefile to output_tree."""
    prefix = output_tree + ".iqtree"
    cmd = [
        "iqtree",
        "--quiet",
        "-redo",
        "--prefix",
        prefix,
        "-m",
        model,
        "-T",
        str(max(1, cpus)),
    ]
    if fast:
        cmd.append("-fast")
    cmd.extend(["-s", input_fasta])
    run_check(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    treefile = Path(f"{prefix}.treefile")
    if not treefile.exists():
        raise FileNotFoundError(
            f"IQ-TREE did not produce expected treefile: {treefile}"
        )
    shutil.copyfile(treefile, output_tree)


def run_species_tree(cfg: Config, input_fasta: str, output_tree: str) -> None:
    """Run selected tree method for the species tree."""
    if cfg.tree_method == "iqtree":
        _run_iqtree(
            input_fasta,
            output_tree,
            cfg.num_cpus,
            cfg.iqtree_model,
            fast=cfg.iqtree_fast,
        )
    else:
        threads = int(os.environ.get("SGTREE_SPECIES_TREE_THREADS", cfg.num_cpus))
        if threads < 1:
            raise ValueError("SGTREE_SPECIES_TREE_THREADS must be at least 1")
        run_fasttree(input_fasta, output_tree, threads)


def _marker_tree_cache_context(
    filepath: str, executable: str
) -> tuple[str, dict[str, str], str] | None:
    """Return canonical marker content, leaf names, and cache key when safe to cache."""
    with Path(filepath).open() as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    canonical_records: list[tuple[str, str]] = []
    leaf_names: dict[str, str] = {}
    for record in records:
        if "|" not in record.id:
            return None
        genome_id = record.id.split("|", 1)[0]
        if (
            not genome_id
            or genome_id in leaf_names
            or not _NEWICK_SAFE_ID.fullmatch(genome_id)
            or not _NEWICK_SAFE_ID.fullmatch(record.id)
        ):
            return None
        leaf_names[genome_id] = record.id
        canonical_records.append((genome_id, str(record.seq)))

    if not canonical_records:
        return None

    canonical_fasta = "".join(
        f">{genome_id}\n{sequence}\n"
        for genome_id, sequence in sorted(canonical_records)
    )
    executable_path = shutil.which(executable) or executable
    executable_file = Path(executable_path)
    executable_identity = str(executable_file.resolve())
    try:
        stat = executable_file.stat()
        executable_identity += f":{stat.st_size}:{stat.st_mtime_ns}"
    except OSError:
        pass
    key_input = (
        f"{executable_identity}\n"
        "VeryFastTree -threads 1 -quiet (with no-threads fallback)\n"
        f"{canonical_fasta}"
    ).encode()
    return canonical_fasta, leaf_names, hashlib.sha256(key_input).hexdigest()


def _materialize_cached_marker_tree(
    cached_tree: str, output_tree: str, leaf_names: dict[str, str]
) -> None:
    tree = Phylo.read(cached_tree, "newick")
    cached_leaf_names = [terminal.name for terminal in tree.get_terminals()]
    if len(cached_leaf_names) != len(set(cached_leaf_names)) or set(
        cached_leaf_names
    ) != set(leaf_names):
        raise ValueError(
            f"Cached marker tree leaf set does not match alignment: {cached_tree}"
        )
    tree_text = Path(cached_tree).read_text()
    for canonical_name, restored_name in leaf_names.items():
        pattern = re.compile(rf"(?<=[(,]){re.escape(canonical_name)}(?=[:),;])")
        tree_text, replacement_count = pattern.subn(restored_name, tree_text)
        if replacement_count != 1:
            raise ValueError(
                f"Could not restore cached leaf {canonical_name}: {cached_tree}"
            )
    Path(output_tree).write_text(tree_text)


def _build_cached_veryfasttree(
    filepath: str, tree_out: str, executable: str, cache_dir: str
) -> bool:
    context = _marker_tree_cache_context(filepath, executable)
    if context is None:
        return False

    canonical_fasta, leaf_names, cache_key = context
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    cached_tree = cache_path / f"{cache_key}.nwk"
    if not cached_tree.exists():
        with tempfile.NamedTemporaryFile(
            mode="w",
            dir=cache_path,
            prefix=f".{cache_key}.",
            suffix=".faa",
            delete=False,
        ) as handle:
            handle.write(canonical_fasta)
            canonical_input = Path(handle.name)
        temporary_tree = Path(f"{canonical_input}.nwk")
        try:
            threaded_cmd = [
                executable,
                "-threads",
                "1",
                "-quiet",
                "-out",
                str(temporary_tree),
                str(canonical_input),
            ]
            fallback_cmd = [
                executable,
                "-quiet",
                "-out",
                str(temporary_tree),
                str(canonical_input),
            ]
            _run_fasttree_with_optional_threads(threaded_cmd, fallback_cmd)
            temporary_tree.replace(cached_tree)
        finally:
            canonical_input.unlink(missing_ok=True)
            temporary_tree.unlink(missing_ok=True)

    _materialize_cached_marker_tree(str(cached_tree), tree_out, leaf_names)
    return True


def _build_tree_worker(args: tuple[str, str, str, str, bool]) -> None:
    """Worker: build a single protein tree with selected tree method."""
    filepath, treeout_dir, tree_method, iqtree_model, iqtree_fast = args
    tree_out = str(Path(treeout_dir) / f"{Path(filepath).name}_tree.out")
    if tree_method == "iqtree":
        _run_iqtree(filepath, tree_out, 1, iqtree_model, fast=iqtree_fast)
    else:
        executable = _fasttree_executable()
        cache_dir = os.environ.get("SGTREE_MARKER_TREE_CACHE_DIR")
        if (
            tree_method == "veryfasttree"
            and cache_dir
            and _build_cached_veryfasttree(filepath, tree_out, executable, cache_dir)
        ):
            return
        threaded_cmd = [
            executable,
            "-threads",
            "1",
            "-quiet",
            "-out",
            tree_out,
            filepath,
        ]
        fallback_cmd = [executable, "-quiet", "-out", tree_out, filepath]
        _run_fasttree_with_optional_threads(threaded_cmd, fallback_cmd)


def run_fasttree_per_marker(
    cfg: Config,
    trimmed_dir: str,
    treeout_dir: str,
) -> None:
    """Build individual protein trees for each trimmed marker alignment."""
    Path(treeout_dir).mkdir(parents=True, exist_ok=True)

    files = [str(path) for path in Path(trimmed_dir).glob("*")]
    args = [
        (f, treeout_dir, cfg.tree_method, cfg.iqtree_model, cfg.iqtree_fast)
        for f in files
    ]
    if not args:
        return

    workers = max(1, min(cfg.num_cpus, len(args)))
    map_threaded(_build_tree_worker, args, workers)
