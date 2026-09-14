"""Build, validate, and reuse reference-genome pipeline caches."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from sgtree._subprocess import run_check
from sgtree.cleanup import (
    RUN_MARKER_CONTENT,
    RUN_MARKER_NAME,
    _archive_directory,
    _finalize_temp_entries,
    _move_to_directory,
    _prepare_temp_directory,
    _zip_file_in_place,
    is_generated_entry,
)
from sgtree.config import Config
from sgtree.id_schema import sanitize_token

REF_CACHE_META_FILE = "cache_meta.json"
REF_CACHE_SCHEMA_VERSION = 2
REQUIRED_TABLE_COLUMNS = ("savedname", "namemodel", "score_bits")
_HASH_CHUNK_SIZE = 1 << 20


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_signature(path: str, *, pattern: str | None = None) -> dict[str, object]:
    source = Path(path)
    if source.is_file():
        root = source.parent
        files = [source]
        source_type = "file"
    elif source.is_dir():
        root = source
        files = sorted(
            candidate
            for candidate in source.glob(pattern or "*")
            if not candidate.name.startswith(".") and candidate.is_file()
        )
        source_type = "directory"
    else:
        raise FileNotFoundError(f"input path does not exist: {source}")

    return {
        "path": str(source.resolve()),
        "type": source_type,
        "files": [
            {
                "path": str(file.relative_to(root)),
                "sha256": _sha256(file),
            }
            for file in files
        ],
    }


def _cache_signature(cfg: Config) -> dict[str, object]:
    if cfg.ref is None:
        raise ValueError("reference cache signature requires cfg.ref")
    return {
        "schema_version": REF_CACHE_SCHEMA_VERSION,
        "reference_input": _source_signature(cfg.ref),
        "model_input": _source_signature(cfg.modeldir, pattern="*.hmm"),
        "percent_models": cfg.percent_models,
        "max_sdup": cfg.max_sdup,
        "max_dupl": cfg.max_dupl,
        "hmmsearch_cutoff": cfg.hmmsearch_cutoff,
        "hmmsearch_evalue": cfg.hmmsearch_evalue,
        "aln_method": "hmmalign",
        "ref_lflt": 0,
        "tree_method": cfg.tree_method,
        "iqtree_fast": cfg.iqtree_fast,
        "iqtree_model": cfg.iqtree_model,
    }


def _legacy_cache_signature(cfg: Config) -> dict[str, object]:
    """Return the exact metadata signature written by schema version 1."""
    return {
        "schema_version": 1,
        "modeldir": str(Path(cfg.modeldir).absolute()),
        "percent_models": cfg.percent_models,
        "max_sdup": cfg.max_sdup,
        "max_dupl": cfg.max_dupl,
        "hmmsearch_cutoff": cfg.hmmsearch_cutoff,
        "hmmsearch_evalue": cfg.hmmsearch_evalue,
        "aln_method": "hmmalign",
        "ref_lflt": 0,
        "tree_method": cfg.tree_method,
        "iqtree_fast": cfg.iqtree_fast,
        "iqtree_model": cfg.iqtree_model,
    }


def _cache_meta_path(ref_dir: str) -> Path:
    return Path(ref_dir) / REF_CACHE_META_FILE


def _write_cache_meta(cfg: Config, ref_dir: str) -> None:
    meta_path = _cache_meta_path(ref_dir)
    with tempfile.TemporaryDirectory(
        dir=meta_path.parent,
        prefix=f".{meta_path.name}.",
    ) as temp_dir:
        temporary_path = Path(temp_dir) / meta_path.name
        with temporary_path.open(
            "x",
            encoding="utf-8",
        ) as handle:
            json.dump(_cache_signature(cfg), handle, sort_keys=True, indent=2)
            handle.write("\n")
        temporary_path.replace(meta_path)


def _load_cache_meta(ref_dir: str) -> dict[str, object] | None:
    meta_path = _cache_meta_path(ref_dir)
    if not meta_path.exists():
        return None
    if meta_path.is_symlink():
        raise ValueError(f"linked cache metadata is not allowed: {meta_path}")
    with meta_path.open(encoding="utf-8") as handle:
        loaded = json.load(handle)
    if not isinstance(loaded, dict):
        raise TypeError("cache metadata root must be an object")
    return loaded


def _validate_table_schema(table_path: str | Path) -> tuple[bool, str]:
    path = Path(table_path)
    if path.is_symlink() or not path.is_file():
        return False, f"missing {path.name}"
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            header = next(csv.reader(handle), [])
    except OSError as exc:
        return False, f"failed to read table schema ({exc})"

    missing = [column for column in REQUIRED_TABLE_COLUMNS if column not in header]
    if missing:
        return False, f"{path.name} missing required columns: {', '.join(missing)}"
    return True, ""


def _validate_nonempty_file(path: Path, relative_path: str) -> tuple[bool, str]:
    if path.is_symlink() or not path.is_file():
        return False, f"missing {relative_path}"
    if path.stat().st_size == 0:
        return False, f"empty {relative_path}"
    return True, ""


def _validate_reference_cache(cfg: Config, ref_dir: str) -> tuple[bool, str]:
    cache_dir = Path(ref_dir)
    if cache_dir.is_symlink():
        return False, "reference cache path is a symbolic link"

    artifact_checks = [
        _validate_nonempty_file(cache_dir / relative_path, relative_path)
        for relative_path in ("hits.hmmout", "proteomes")
    ]
    artifact_checks.extend(
        _validate_table_schema(cache_dir / relative_path)
        for relative_path in ("table_elim_dups", "tables/merged_final")
    )
    for artifact_valid, artifact_reason in artifact_checks:
        if not artifact_valid:
            return False, artifact_reason

    try:
        metadata = _load_cache_meta(ref_dir)
    except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
        return False, f"invalid {REF_CACHE_META_FILE} ({exc})"

    try:
        expected = _cache_signature(cfg)
    except OSError as exc:
        return False, f"failed to fingerprint cache inputs ({exc})"
    if metadata != expected:
        return False, "cache metadata does not match inputs or settings"

    return True, ""


def _mark_legacy_cache_for_rebuild(cfg: Config, ref_dir: Path) -> None:
    """Adopt an intact schema-1 cache so the nested CLI can clean it safely."""
    try:
        metadata = _load_cache_meta(str(ref_dir))
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return
    legacy_keys = _legacy_cache_signature(cfg).keys()
    if (
        metadata is None
        or metadata.get("schema_version") != 1
        or metadata.keys() != legacy_keys
    ):
        return

    artifact_checks = [
        _validate_nonempty_file(ref_dir / relative_path, relative_path)
        for relative_path in ("hits.hmmout", "proteomes")
    ]
    artifact_checks.extend(
        _validate_table_schema(ref_dir / relative_path)
        for relative_path in ("table_elim_dups", "tables/merged_final")
    )
    if not all(valid for valid, _reason in artifact_checks):
        return

    marker = ref_dir / RUN_MARKER_NAME
    if marker.exists() or marker.is_symlink():
        return
    with marker.open("x", encoding="utf-8") as handle:
        handle.write(RUN_MARKER_CONTENT)


def _genome_ids_from_dir(dirpath: str) -> list[str]:
    """Derive genome IDs from proteome filenames, matching input staging."""
    paths = sorted(
        path
        for path in Path(dirpath).glob("*")
        if not path.name.startswith(".") and path.is_file()
    )
    return [
        sanitize_token(path.stem, f"genome_{index:05d}")
        for index, path in enumerate(paths, start=1)
    ]


def _genome_ids(input_path: str) -> list[str]:
    if Path(input_path).is_dir():
        return _genome_ids_from_dir(input_path)

    genome_ids: list[str] = []
    with Path(input_path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            genome_id = line[1:].split("|", maxsplit=1)[0]
            if genome_id not in genome_ids:
                genome_ids.append(genome_id)
    return genome_ids


def check_duplicate_proteomes(genomedir: str, refdir: str) -> list[str]:
    """Exit when query and reference inputs contain the same genome ID."""
    genome_ids = _genome_ids(genomedir)
    reference_ids = _genome_ids(refdir)
    common = [genome_id for genome_id in genome_ids if genome_id in reference_ids]
    if not common:
        return reference_ids

    print(
        "WARNING: Duplicate proteomes found in reference directory, "
        "this will cause errors for --marker_selection yes (Noperm)\n"
        "please delete the following from either directory/concat file and try again."
    )
    for genome_id in common:
        print("please delete", genome_id)
    raise SystemExit(1)


def _archive_reference_outputs(ref_dir: Path) -> None:
    archives: list[Path] = []
    keep_names = {
        ".sgtree-run",
        "tables",
        "concat",
        "extracted_seqs",
        "marker_count_matrix.csv",
        "proteomes",
        "hits.hmmout",
        "table_elim_dups",
        "genome_manifest.tsv",
    }
    for path in ref_dir.iterdir():
        if path.name in keep_names or not is_generated_entry(path.name):
            continue
        if path.is_dir():
            archive = _archive_directory(str(path))
            if archive is not None:
                archives.append(archive)
            continue
        _zip_file_in_place(str(path))

    temp_dir = _prepare_temp_directory(ref_dir)
    _finalize_temp_entries(
        ref_dir,
        temp_dir,
        archives,
        ("ref_and_query_proteomes", "models", "tree.nwk"),
    )
    for path in ref_dir.glob("*.txt"):
        if is_generated_entry(path.name):
            _move_to_directory(path, temp_dir)

    for path in temp_dir.iterdir():
        if zipfile.is_zipfile(path):
            continue
        _zip_file_in_place(str(path))


def prepare_reference(cfg: Config) -> list[str] | None:
    """Build or reuse the reference cache and return reference filenames."""
    if cfg.ref is None:
        print("no reference directory")
        return None

    reference_ids = check_duplicate_proteomes(cfg.genomedir, cfg.ref)
    reference_filenames = [f"{genome_id}.faa" for genome_id in reference_ids]
    print(reference_filenames)

    ref_dir_value = cfg.ref_dir_path()
    assert ref_dir_value is not None
    ref_dir = Path(ref_dir_value)
    cache_valid, cache_reason = _validate_reference_cache(cfg, str(ref_dir))
    if cache_valid:
        print(f"already have files for references at {ref_dir}/")
        return reference_filenames

    if ref_dir.is_symlink():
        raise FileExistsError(f"refusing to replace linked reference cache: {ref_dir}")
    if ref_dir.exists():
        if not ref_dir.is_dir():
            raise FileExistsError(f"reference cache path is not a directory: {ref_dir}")
        print(f"reference cache invalid ({cache_reason}); rebuilding {ref_dir}/")
        _mark_legacy_cache_for_rebuild(cfg, ref_dir)

    command = [
        sys.executable,
        "-m",
        "sgtree",
        cfg.ref,
        cfg.modeldir,
        "--num_cpus",
        str(cfg.num_cpus),
        "--percent_models",
        str(cfg.percent_models),
        "--lflt",
        "0",
        "--aln",
        "hmmalign",
        "--tree_method",
        cfg.tree_method,
        "--iqtree_fast",
        "yes" if cfg.iqtree_fast else "no",
        "--iqtree_model",
        cfg.iqtree_model,
        "--hmmsearch_cutoff",
        cfg.hmmsearch_cutoff,
        "--hmmsearch_evalue",
        str(cfg.hmmsearch_evalue),
        "--max_sdup",
        str(cfg.max_sdup),
        "--max_dupl",
        str(cfg.max_dupl),
        "--save_dir",
        str(ref_dir),
        "--is_ref",
        "yes",
    ]
    print("- ... Creating new reference directory\n", command)
    run_check(command, stdout=subprocess.PIPE)

    _archive_reference_outputs(ref_dir)
    _write_cache_meta(cfg, str(ref_dir))
    return reference_filenames
