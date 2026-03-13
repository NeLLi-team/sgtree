from __future__ import annotations

import csv
import multiprocessing as mp
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from Bio import SeqIO


FASTA_SUFFIXES = {".fna", ".fa", ".fasta"}


@dataclass(frozen=True)
class GenomeRecord:
    genome_id: str
    source_role: str
    input_format: str
    source_file: str
    assembly_path: str | None
    staged_proteome_path: str | None
    total_bases: int
    contigs: int

    @property
    def has_assembly(self) -> bool:
        return bool(self.assembly_path)


@dataclass(frozen=True)
class DirectionalAni:
    ani: float
    query_coverage: float
    target_coverage: float


@dataclass(frozen=True)
class SamProjection:
    coverage: dict[str, list[tuple[int, int]]]
    mismatches: dict[tuple[str, int], str]


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return []
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def _interval_length(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start for start, end in _merge_intervals(intervals))


def _position_is_covered(intervals: list[tuple[int, int]], pos: int) -> bool:
    lo = 0
    hi = len(intervals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        start, end = intervals[mid]
        if start <= pos < end:
            return True
        if pos < start:
            hi = mid - 1
        else:
            lo = mid + 1
    return False


def _parse_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in {None, ""} else default
    except ValueError:
        return default


def _parse_int(value: str | None, default: int = 0) -> int:
    try:
        return int(value) if value not in {None, ""} else default
    except ValueError:
        return default


def load_genome_manifest(manifest_path: str | Path, source_role: str) -> list[GenomeRecord]:
    manifest_path = Path(manifest_path)
    rows: list[GenomeRecord] = []
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            assembly_path = row.get("assembly_path") or None
            staged_proteome_path = row.get("staged_proteome_path") or None
            rows.append(
                GenomeRecord(
                    genome_id=str(row["genome_id"]),
                    source_role=source_role,
                    input_format=str(row.get("input_format", "unknown")),
                    source_file=str(row.get("source_file", "")),
                    assembly_path=assembly_path,
                    staged_proteome_path=staged_proteome_path,
                    total_bases=_parse_int(row.get("total_bases")),
                    contigs=_parse_int(row.get("contigs")),
                )
            )
    if not rows:
        raise ValueError(f"No genomes found in manifest: {manifest_path}")
    return rows


def load_cluster_records(
    query_manifest: str | Path,
    ref_manifest: str | Path | None = None,
) -> list[GenomeRecord]:
    records = load_genome_manifest(query_manifest, "query")
    if ref_manifest is not None:
        records.extend(load_genome_manifest(ref_manifest, "ref"))
    genome_ids = [record.genome_id for record in records]
    if len(set(genome_ids)) != len(genome_ids):
        duplicates = sorted({gid for gid in genome_ids if genome_ids.count(gid) > 1})
        raise ValueError(f"Duplicate genome ids across query/ref manifests: {', '.join(duplicates)}")
    return records


def resolve_ani_backend(requested: str) -> str:
    backend = requested.strip().lower()
    if backend == "auto":
        if shutil.which("skani"):
            return "skani"
        if shutil.which("minimap2"):
            return "minimap2"
        raise RuntimeError("ANI clustering requires either skani or minimap2 in PATH")
    if backend == "skani":
        if not shutil.which("skani"):
            raise RuntimeError("ANI backend 'skani' requested but skani is not available")
        return "skani"
    if backend == "minimap2":
        if not shutil.which("minimap2"):
            raise RuntimeError("ANI backend 'minimap2' requested but minimap2 is not available")
        return "minimap2"
    raise ValueError(f"Unsupported ANI backend: {requested}")


def _run_cmd(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def _minimap2_direction(target: GenomeRecord, query: GenomeRecord) -> DirectionalAni:
    if not target.assembly_path or not query.assembly_path:
        raise RuntimeError("minimap2 ANI requires assembly paths for both genomes")
    result = _run_cmd(
        [
            "minimap2",
            "-x",
            "asm5",
            "-c",
            "-t",
            "1",
            target.assembly_path,
            query.assembly_path,
        ]
    )
    matches = 0
    aligned = 0
    query_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    target_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)

    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        query_intervals[fields[0]].append((_parse_int(fields[2]), _parse_int(fields[3])))
        target_intervals[fields[5]].append((_parse_int(fields[7]), _parse_int(fields[8])))
        matches += _parse_int(fields[9])
        aligned += _parse_int(fields[10])

    ani = matches / aligned if aligned else 0.0
    query_cov = (
        sum(_interval_length(intervals) for intervals in query_intervals.values()) / query.total_bases
        if query.total_bases
        else 0.0
    )
    target_cov = (
        sum(_interval_length(intervals) for intervals in target_intervals.values()) / target.total_bases
        if target.total_bases
        else 0.0
    )
    return DirectionalAni(ani=ani, query_coverage=query_cov, target_coverage=target_cov)


def _symmetrize_ani(a_to_b: DirectionalAni, b_to_a: DirectionalAni) -> tuple[float, float]:
    ani_values = [metric.ani for metric in (a_to_b, b_to_a) if metric.ani > 0]
    ani = sum(ani_values) / len(ani_values) if ani_values else 0.0
    aligned_fraction = min(
        a_to_b.query_coverage,
        a_to_b.target_coverage,
        b_to_a.query_coverage,
        b_to_a.target_coverage,
    )
    return ani, aligned_fraction


def _pairwise_minimap2(args: tuple[GenomeRecord, GenomeRecord]) -> dict[str, object]:
    genome_a, genome_b = args
    a_to_b = _minimap2_direction(genome_a, genome_b)
    b_to_a = _minimap2_direction(genome_b, genome_a)
    ani, aligned_fraction = _symmetrize_ani(a_to_b, b_to_a)
    return {
        "genome_a": genome_a.genome_id,
        "genome_b": genome_b.genome_id,
        "backend": "minimap2",
        "ani": ani,
        "aligned_fraction": aligned_fraction,
        "ani_a_to_b": a_to_b.ani,
        "ani_b_to_a": b_to_a.ani,
        "coverage_a_to_b": min(a_to_b.query_coverage, a_to_b.target_coverage),
        "coverage_b_to_a": min(b_to_a.query_coverage, b_to_a.target_coverage),
    }


def _map_with_fallback(func, args: list[tuple], workers: int) -> list[dict[str, object]]:
    if not args:
        return []
    n_workers = max(1, min(workers, len(args)))
    if n_workers == 1:
        return [func(item) for item in args]
    try:
        with mp.Pool(n_workers) as pool:
            return pool.map(func, args)
    except (PermissionError, OSError) as exc:
        print(f"warning: multiprocessing unavailable ({exc}); falling back to serial ANI execution")
        return [func(item) for item in args]


def _compute_skani_pairs(records: list[GenomeRecord], cpus: int) -> list[dict[str, object]]:
    cmd = ["skani", "triangle", "--medium", "--min-af", "0", "-E", "-t", str(max(1, cpus))]
    cmd.extend(record.assembly_path for record in records if record.assembly_path)
    result = _run_cmd(cmd)
    path_to_id = {
        str(Path(record.assembly_path).resolve()): record.genome_id
        for record in records
        if record.assembly_path
    }
    rows: list[dict[str, object]] = []
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith("Ref_file"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        ref_file = str(Path(fields[0]).resolve())
        query_file = str(Path(fields[1]).resolve())
        ani = _parse_float(fields[2])
        af_ref = _parse_float(fields[3])
        af_query = _parse_float(fields[4])
        if ani > 1.0:
            ani /= 100.0
        if af_ref > 1.0:
            af_ref /= 100.0
        if af_query > 1.0:
            af_query /= 100.0
        rows.append(
            {
                "genome_a": path_to_id[ref_file],
                "genome_b": path_to_id[query_file],
                "backend": "skani",
                "ani": ani,
                "aligned_fraction": min(af_ref, af_query),
                "ani_a_to_b": ani,
                "ani_b_to_a": ani,
                "coverage_a_to_b": af_ref,
                "coverage_b_to_a": af_query,
            }
        )
    return rows


def _run_mcl_binary(
    labels: list[str],
    rows: list[dict[str, object]],
    *,
    ani_threshold: float,
    inflation: float,
    graph_path: Path,
    out_path: Path,
) -> list[list[str]]:
    if not shutil.which("mcl"):
        return _run_python_mcl(
            labels,
            rows,
            ani_threshold=ani_threshold,
            inflation=inflation,
        )

    graph_path.parent.mkdir(parents=True, exist_ok=True)
    retained = 0
    with graph_path.open("w") as handle:
        for row in rows:
            ani = float(row["ani"])
            if ani < ani_threshold:
                continue
            handle.write(f"{row['genome_a']}\t{row['genome_b']}\t{ani:.6f}\n")
            retained += 1

    if retained == 0:
        out_path.write_text("")
        return [[label] for label in labels]

    subprocess.run(
        [
            "mcl",
            str(graph_path),
            "--abc",
            "-I",
            str(inflation),
            "-o",
            str(out_path),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    clusters: list[list[str]] = []
    seen: set[str] = set()
    for line in out_path.read_text().splitlines():
        members = sorted({token for token in line.split("\t") if token})
        if not members:
            continue
        clusters.append(members)
        seen.update(members)
    for label in labels:
        if label not in seen:
            clusters.append([label])
    clusters.sort(key=lambda members: (-len(members), members[0]))
    return clusters


def compute_pairwise_ani(
    records: list[GenomeRecord],
    *,
    backend: str,
    cpus: int,
) -> list[dict[str, object]]:
    if any(not record.has_assembly for record in records):
        missing = sorted(record.genome_id for record in records if not record.has_assembly)
        raise RuntimeError(
            "ANI clustering requires assembly input (*.fna/*.fa/*.fasta) for all genomes; missing assemblies for: "
            + ", ".join(missing)
        )
    if backend == "skani":
        return _compute_skani_pairs(records, cpus)
    pair_args = [
        (records[idx], records[jdx])
        for idx in range(len(records))
        for jdx in range(idx + 1, len(records))
    ]
    return _map_with_fallback(_pairwise_minimap2, pair_args, cpus)


def _run_python_mcl(
    labels: list[str],
    rows: list[dict[str, object]],
    *,
    ani_threshold: float,
    inflation: float,
    max_iterations: int = 100,
    expansion: int = 2,
    prune_threshold: float = 1e-5,
    tol: float = 1e-6,
) -> list[list[str]]:
    label_to_idx = {label: idx for idx, label in enumerate(labels)}
    matrix = np.zeros((len(labels), len(labels)), dtype=float)
    np.fill_diagonal(matrix, 1.0)

    for row in rows:
        ani = float(row["ani"])
        if ani < ani_threshold:
            continue
        idx = label_to_idx[str(row["genome_a"])]
        jdx = label_to_idx[str(row["genome_b"])]
        matrix[idx, jdx] = max(matrix[idx, jdx], ani)
        matrix[jdx, idx] = max(matrix[jdx, idx], ani)

    if not np.any(matrix):
        return [[label] for label in labels]

    matrix = matrix / np.where(matrix.sum(axis=0, keepdims=True) == 0, 1.0, matrix.sum(axis=0, keepdims=True))
    for _ in range(max_iterations):
        previous = matrix.copy()
        matrix = np.linalg.matrix_power(matrix, expansion)
        matrix = np.power(matrix, inflation)
        matrix[matrix < prune_threshold] = 0.0
        col_sums = matrix.sum(axis=0, keepdims=True)
        col_sums[col_sums == 0] = 1.0
        matrix = matrix / col_sums
        if np.allclose(matrix, previous, atol=tol):
            break

    adjacency = (matrix > prune_threshold).astype(int)
    clusters: list[list[str]] = []
    seen: set[int] = set()
    for start in range(len(labels)):
        if start in seen:
            continue
        stack = [start]
        component: list[str] = []
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            component.append(labels[current])
            neighbors = np.where(adjacency[current] | adjacency[:, current])[0]
            stack.extend(int(idx) for idx in neighbors if int(idx) not in seen)
        clusters.append(sorted(component))
    clusters.sort(key=lambda members: (-len(members), members[0]))
    return clusters


def choose_cluster_representative(cluster_members: list[GenomeRecord]) -> GenomeRecord:
    return min(
        cluster_members,
        key=lambda record: (
            0 if record.source_role == "ref" else 1,
            record.contigs if record.contigs > 0 else 10**9,
            -record.total_bases,
            record.genome_id,
        ),
    )


def run_ani_clustering(
    *,
    query_manifest: str | Path,
    ref_manifest: str | Path | None,
    outdir: str | Path,
    ani_threshold: float,
    inflation: float,
    backend: str,
    cpus: int,
) -> dict[str, str]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = load_cluster_records(query_manifest, ref_manifest)
    resolved_backend = resolve_ani_backend(backend)
    pair_rows = compute_pairwise_ani(records, backend=resolved_backend, cpus=cpus)
    pair_rows.sort(key=lambda row: (str(row["genome_a"]), str(row["genome_b"])))

    pairwise_path = outdir / "ani_pairwise.tsv"
    with pairwise_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "genome_a",
                "genome_b",
                "backend",
                "ani",
                "aligned_fraction",
                "ani_a_to_b",
                "ani_b_to_a",
                "coverage_a_to_b",
                "coverage_b_to_a",
            ],
        )
        writer.writeheader()
        writer.writerows(pair_rows)

    labels = sorted(record.genome_id for record in records)
    record_map = {record.genome_id: record for record in records}
    graph_path = outdir / "ani_graph.tsv"
    mcl_path = outdir / "ani_mcl_clusters.txt"
    clusters = _run_mcl_binary(
        labels,
        pair_rows,
        ani_threshold=ani_threshold,
        inflation=inflation,
        graph_path=graph_path,
        out_path=mcl_path,
    )

    clusters_path = outdir / "ani_clusters.tsv"
    representatives_path = outdir / "ani_representatives.tsv"
    kept_path = outdir / "ani_kept_genomes.txt"
    cluster_rows: list[dict[str, object]] = []
    representative_rows: list[dict[str, object]] = []
    kept_genomes: list[str] = []

    for idx, members in enumerate(clusters, start=1):
        member_records = [record_map[member] for member in members]
        representative = choose_cluster_representative(member_records)
        cluster_id = f"ani_cluster_{idx:03d}"
        kept_genomes.append(representative.genome_id)
        representative_rows.append(
            {
                "cluster_id": cluster_id,
                "representative_genome": representative.genome_id,
                "source_role": representative.source_role,
                "cluster_size": len(members),
                "total_bases": representative.total_bases,
                "contigs": representative.contigs,
            }
        )
        for member in member_records:
            cluster_rows.append(
                {
                    "cluster_id": cluster_id,
                    "representative_genome": representative.genome_id,
                    "genome_id": member.genome_id,
                    "source_role": member.source_role,
                    "cluster_size": len(members),
                    "total_bases": member.total_bases,
                    "contigs": member.contigs,
                    "kept_for_species_tree": "yes" if member.genome_id == representative.genome_id else "no",
                    "backend": resolved_backend,
                    "ani_threshold": f"{ani_threshold:.4f}",
                    "mcl_inflation": f"{inflation:.2f}",
                }
            )

    with clusters_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "cluster_id",
                "representative_genome",
                "genome_id",
                "source_role",
                "cluster_size",
                "total_bases",
                "contigs",
                "kept_for_species_tree",
                "backend",
                "ani_threshold",
                "mcl_inflation",
            ],
        )
        writer.writeheader()
        writer.writerows(cluster_rows)

    with representatives_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "cluster_id",
                "representative_genome",
                "source_role",
                "cluster_size",
                "total_bases",
                "contigs",
            ],
        )
        writer.writeheader()
        writer.writerows(representative_rows)

    kept_path.write_text("\n".join(sorted(kept_genomes)) + "\n")
    return {
        "clusters": str(clusters_path),
        "pairwise": str(pairwise_path),
        "representatives": str(representatives_path),
        "kept": str(kept_path),
        "backend": resolved_backend,
    }


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    value = ""
    parsed: list[tuple[int, str]] = []
    for char in cigar:
        if char.isdigit():
            value += char
            continue
        if not value:
            raise ValueError(f"Malformed CIGAR string: {cigar}")
        parsed.append((int(value), char))
        value = ""
    if value:
        raise ValueError(f"Malformed CIGAR string: {cigar}")
    return parsed


def _project_alignment_to_reference(
    reference_path: str,
    query_path: str,
) -> SamProjection:
    result = _run_cmd(
        [
            "minimap2",
            "-a",
            "--eqx",
            "-x",
            "asm5",
            "-t",
            "1",
            reference_path,
            query_path,
        ]
    )
    coverage: dict[str, list[tuple[int, int]]] = defaultdict(list)
    mismatches: dict[tuple[str, int], str] = {}

    for line in result.stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = _parse_int(fields[1])
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        rname = fields[2]
        if rname == "*" or fields[5] == "*":
            continue
        sequence = fields[9]
        if flag & 16:
            sequence = _reverse_complement(sequence)
        qpos = 0
        rpos = _parse_int(fields[3]) - 1
        for length, op in _parse_cigar(fields[5]):
            if op in {"=", "X", "M"}:
                coverage[rname].append((rpos, rpos + length))
                if op == "X":
                    for offset in range(length):
                        mismatches[(rname, rpos + offset)] = sequence[qpos + offset].upper()
                qpos += length
                rpos += length
            elif op in {"I", "S"}:
                qpos += length
            elif op in {"D", "N"}:
                rpos += length
            elif op in {"H", "P"}:
                continue
            else:
                raise ValueError(f"Unsupported CIGAR op '{op}' in {fields[5]}")

    merged = {
        contig: _merge_intervals(intervals)
        for contig, intervals in coverage.items()
    }
    return SamProjection(coverage=merged, mismatches=mismatches)


def _write_star_tree(genome_ids: list[str], tree_path: Path) -> None:
    tree_path.write_text("(" + ",".join(genome_ids) + ");\n")


def _run_fasttree_nt(alignment_path: Path, tree_path: Path) -> None:
    _run_cmd(["FastTree", "-nt", "-gtr", "-quiet", "-out", str(tree_path), str(alignment_path)])


def _load_cluster_members(clusters_path: str | Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with Path(clusters_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            grouped[str(row["cluster_id"])].append({key: str(value) for key, value in row.items()})
    return grouped


def build_snp_trees(
    *,
    clusters_path: str | Path,
    query_manifest: str | Path,
    ref_manifest: str | Path | None,
    outdir: str | Path,
    min_cluster_size: int,
    cpus: int,
) -> str:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = load_cluster_records(query_manifest, ref_manifest)
    record_map = {record.genome_id: record for record in records}
    cluster_members = _load_cluster_members(clusters_path)
    summary_rows: list[dict[str, object]] = []

    for cluster_id, members in sorted(cluster_members.items()):
        representative = str(members[0]["representative_genome"])
        genome_ids = [str(row["genome_id"]) for row in sorted(members, key=lambda row: row["genome_id"])]
        if len(genome_ids) < min_cluster_size:
            summary_rows.append(
                {
                    "cluster_id": cluster_id,
                    "representative_genome": representative,
                    "cluster_size": len(genome_ids),
                    "status": "skipped_cluster_too_small",
                    "snp_sites": 0,
                    "alignment_path": "",
                    "tree_path": "",
                    "members_path": "",
                }
            )
            continue

        representative_record = record_map[representative]
        if not representative_record.assembly_path:
            raise RuntimeError(f"Representative genome {representative} has no assembly path")

        cluster_dir = outdir / cluster_id
        cluster_dir.mkdir(parents=True, exist_ok=True)
        members_path = cluster_dir / "members.tsv"
        with members_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(members[0].keys()))
            writer.writeheader()
            writer.writerows(members)

        with open(representative_record.assembly_path) as handle:
            reference_sequences = {
                record.id: str(record.seq).upper()
                for record in SeqIO.parse(handle, "fasta")
            }
        projections: dict[str, SamProjection] = {
            representative: SamProjection(
                coverage={contig: [(0, len(seq))] for contig, seq in reference_sequences.items()},
                mismatches={},
            )
        }
        for genome_id in genome_ids:
            if genome_id == representative:
                continue
            member_record = record_map[genome_id]
            if not member_record.assembly_path:
                raise RuntimeError(f"Cluster member genome {genome_id} has no assembly path")
            projections[genome_id] = _project_alignment_to_reference(
                representative_record.assembly_path,
                member_record.assembly_path,
            )

        candidate_sites = set()
        for genome_id in genome_ids:
            if genome_id == representative:
                continue
            candidate_sites.update(projections[genome_id].mismatches.keys())
        candidate_sites = sorted(candidate_sites)

        site_columns: list[tuple[str, int, dict[str, str]]] = []
        for contig, pos in candidate_sites:
            reference_base = reference_sequences.get(contig, "")[pos : pos + 1].upper()
            if reference_base not in {"A", "C", "G", "T"}:
                continue
            column: dict[str, str] = {}
            valid = True
            for genome_id in genome_ids:
                projection = projections[genome_id]
                intervals = projection.coverage.get(contig, [])
                if not _position_is_covered(intervals, pos):
                    valid = False
                    break
                base = projection.mismatches.get((contig, pos), reference_base)
                if base not in {"A", "C", "G", "T"}:
                    valid = False
                    break
                column[genome_id] = base
            if valid and len(set(column.values())) > 1:
                site_columns.append((contig, pos, column))

        tree_path = cluster_dir / "tree.nwk"
        alignment_path = cluster_dir / "core_snps.fna"
        alignment_path_str = str(alignment_path)
        if site_columns:
            with alignment_path.open("w") as handle:
                for genome_id in genome_ids:
                    handle.write(f">{genome_id}\n")
                    handle.write("".join(column[genome_id] for _, _, column in site_columns) + "\n")
            _run_fasttree_nt(alignment_path, tree_path)
            status = "built"
        else:
            _write_star_tree(genome_ids, tree_path)
            status = "no_snp_sites"
            alignment_path_str = ""

        summary_rows.append(
            {
                "cluster_id": cluster_id,
                "representative_genome": representative,
                "cluster_size": len(genome_ids),
                "status": status,
                "snp_sites": len(site_columns),
                "alignment_path": alignment_path_str,
                "tree_path": str(tree_path),
                "members_path": str(members_path),
            }
        )

    summary_path = outdir / "snp_tree_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "cluster_id",
                "representative_genome",
                "cluster_size",
                "status",
                "snp_sites",
                "alignment_path",
                "tree_path",
                "members_path",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    return str(summary_path)
