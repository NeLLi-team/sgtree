import os
import glob
import csv
import logging
import math
import shutil
import statistics
import numpy as np

import pandas as pd
from Bio import SeqIO
from ete3 import Tree

from sgtree.config import Config
from sgtree.id_schema import parse_savedname, parse_sequence_id
from sgtree.marker_selection.loo_profile import score_loo_profiles
from sgtree.parallel import map_processed, map_threaded


logger = logging.getLogger("sgtree")

SCORE_COLUMNS = ("score_bits", "7")
LEGACY_SINGLETON_MODE_ALIASES = {
    "outlier": "topoknn",
}
COMPOSITE_SCORE_THRESHOLD = 2.0
COMPOSITE_SCORE_MARGIN = 0.35
RECIPIENT_CONSENSUS_Z_THRESHOLD = 3.0
RECIPIENT_CONSENSUS_MIN_SCORE = 1.5
RECIPIENT_CONSENSUS_RANK_MARGIN = 0.15
CONTIG_CONSENSUS_HIGH_OVERLAP = 0.6
CONTIG_CONSENSUS_LOW_OVERLAP = 0.2
CONTIG_CONSENSUS_MIN_SUPPORT = 2
NEIGHBOR_CLADE_MIN_PRESENT = 3
NEIGHBOR_CLADE_MIN_PRESENT_FRACTION = 0.6
NEIGHBOR_CLADE_MIN_SPECIES_ANCHOR = 0.65
NEIGHBOR_CLADE_MIN_NEIGHBOR_PURITY = 0.75
NEIGHBOR_CLADE_MIN_PURITY_DROP = 0.35
NEIGHBOR_CLADE_MAX_KNN_AGREEMENT = 0.4
NEIGHBOR_CLADE_MIN_SCORE = 2.0
NEIGHBOR_CLADE_MIN_RECIPIENT_SUPPORT = 1.0
NEIGHBOR_ML_HIGH_CONF_MIN_PRESENT = 3
NEIGHBOR_ML_HIGH_CONF_MIN_FRACTION = 0.6
NEIGHBOR_ML_HIGH_CONF_MIN_ANCHOR = 0.35
UNKNOWN_CONTIG_DELTA_FLOOR = 0.05
UNKNOWN_CONTIG_TOPOKNN_FLOOR = 1.5
UNKNOWN_CONTIG_RECIPIENT_FLOOR = 1.5

# GCP (Genome Consistency Profiling) constants
GCP_KEY_FEATURES = [
    "neighbor_overlap",
    "delta_rf",
    "attachment_gap",
    "recipient_consensus_score",
    "topoknn_score",
    "anchor_knn_agreement",
    "neighbor_clade_score",
]
GCP_COMBINED_THRESHOLD = 0.80
GCP_OUTLIER_COUNT_THRESHOLD = 3
GCP_IFOREST_THRESHOLD = 0.5
GCP_Z_THRESHOLD = 1.5
GCP_MIN_GENOMES = 3
GCP_MIN_MARKERS = 5


def choose_best_candidate(
    candidates: list[dict],
    preferred_protein_id: str | None = None,
) -> dict:
    """Rank RF-selection candidates deterministically.

    Preference order:
    1. lower RF distance
    2. more informative splits in the RF comparison
    3. stronger same-contig marker support
    4. higher HMM bitscore
    5. lexicographically smaller protein id for stability

    When a previously selected protein is supplied and remains tied on the
    RF-based criteria, keep that assignment instead of flipping on bitscore.
    """
    if not candidates:
        raise ValueError("Expected at least one RF candidate")
    best_rf_key = min(
        (float(item["rf_distance"]), -int(item.get("informative_splits", 0)))
        for item in candidates
    )
    if preferred_protein_id is not None:
        for item in candidates:
            rf_key = (
                float(item["rf_distance"]),
                -int(item.get("informative_splits", 0)),
            )
            if rf_key == best_rf_key and str(item["protein_id"]) == str(preferred_protein_id):
                return item
    return min(
        candidates,
        key=lambda item: (
            float(item["rf_distance"]),
            -int(item.get("informative_splits", 0)),
            -int(item.get("contig_marker_support", 0)),
            -float(item["bitscore"]),
            str(item["protein_id"]),
        ),
    )


def _load_score_table(table_path: str) -> tuple[pd.DataFrame, str]:
    dfa = pd.read_csv(table_path)
    if "savedname" not in dfa.columns:
        raise ValueError(f"Missing required column 'savedname' in {table_path}")
    score_col = next((col for col in SCORE_COLUMNS if col in dfa.columns), None)
    if score_col is None:
        raise ValueError(
            f"Missing score column in {table_path}; expected one of: {', '.join(SCORE_COLUMNS)}"
        )
    dfa = dfa.set_index("savedname")
    return dfa, score_col


_SCORE_TABLE_CACHE: dict[str, tuple[pd.DataFrame, str]] = {}
_SPECIES_TREE_CACHE: dict[str, Tree] = {}


def _cached_score_table(table_path: str) -> tuple[pd.DataFrame, str]:
    """Return ``(score_table, score_col)`` for ``table_path``.

    The parsed DataFrame is loaded once per process and reused; callers
    treat the DataFrame as read-only. Workers dispatched via
    ``map_processed`` populate the cache lazily on first hit and reuse it
    across subsequent marker tasks in the same worker lifetime.
    """
    cached = _SCORE_TABLE_CACHE.get(table_path)
    if cached is None:
        cached = _load_score_table(table_path)
        _SCORE_TABLE_CACHE[table_path] = cached
    return cached


def _cached_species_tree(species_tree_path: str) -> Tree:
    """Return a *copy* of the parsed species tree at ``species_tree_path``.

    The underlying ``ete3.Tree`` is parsed once per process. Each caller
    gets a defensive copy because several callers mutate their local tree
    (e.g. ``_propose_singleton_prune_worker`` and
    ``build_singleton_output_tree`` call ``.prune(...)``). Tree.copy() is
    O(nodes) but orders of magnitude cheaper than Newick parsing.
    """
    cached = _SPECIES_TREE_CACHE.get(species_tree_path)
    if cached is None:
        cached = Tree(species_tree_path)
        _SPECIES_TREE_CACHE[species_tree_path] = cached
    return cached.copy()


def _get_ascore(identifier: str, score_table: pd.DataFrame, score_col: str) -> str:
    row = score_table.loc[identifier.replace("|", "/")]
    return row.name + ":" + str(float(row[score_col]))


def _best_score(scored_list):
    """Return the entry with the highest score from a list of 'name:score' strings."""
    if not scored_list:
        return None
    return min(
        scored_list,
        key=lambda entry: (-float(entry.rsplit(":", 1)[1]), entry.rsplit(":", 1)[0]),
    )


def _removekey(d, key):
    """Return a copy of dict d without the given key."""
    return {k: v for k, v in d.items() if k != key}


def _split_scored_entry(scored_entry: str) -> tuple[str, float]:
    protein_id, bitscore = scored_entry.rsplit(":", 1)
    return protein_id, float(bitscore)


def _build_duplicate_map(lst_nodes: list[str], score_table: pd.DataFrame, score_col: str) -> dict[str, list[str]]:
    dups: dict[str, list[str]] = {}
    for node_name in lst_nodes:
        genome = node_name.split("|")[0]
        dups.setdefault(genome, []).append(_get_ascore(node_name, score_table, score_col))
    return dups


def _build_contig_support_map(score_table: pd.DataFrame) -> dict[str, int]:
    df = score_table.reset_index().copy()
    if "namemodel" not in df.columns:
        return {str(row.savedname): 1 for row in df.itertuples(index=False)}
    parsed = df["savedname"].apply(parse_savedname)
    df["genome_id"] = parsed.apply(lambda item: item[0])
    df["contig_id"] = parsed.apply(lambda item: item[1])
    df["marker_id"] = df["namemodel"].astype(str).str.split("/").str[-1]
    grouped = (
        df.groupby(["genome_id", "contig_id"])["marker_id"]
        .nunique()
        .to_dict()
    )
    support: dict[str, int] = {}
    for row in df.itertuples(index=False):
        key = (str(row.genome_id), str(row.contig_id))
        support[str(row.savedname)] = int(grouped.get(key, 1))
    return support


def _evaluate_candidate(
    marker_tree: Tree,
    species_tree: Tree,
    dups: dict[str, list[str]],
    selected: dict[str, str],
    genome: str,
    scored_entry: str,
    contig_support_map: dict[str, int],
) -> dict:
    alt_scored = []
    for group, entries in dups.items():
        chosen = scored_entry if group == genome else selected.get(group)
        if chosen is None and entries:
            chosen = _best_score(entries)
        if chosen is not None:
            alt_scored.append(chosen)

    alt = [entry.rsplit(":", 1)[0].replace("/", "|") for entry in alt_scored]
    species_tree_copy = species_tree.copy()
    species_tree_copy.prune([name.split("|")[0] for name in alt])

    prot_tree = marker_tree.copy()
    prot_tree.prune(alt)
    prot_tree_copy = prot_tree.copy()
    for leaf in prot_tree_copy.iter_leaves():
        leaf.name = leaf.name.split("|")[0]

    rf, maxrf, *_ = species_tree_copy.robinson_foulds(prot_tree_copy, unrooted_trees=True)
    protein_id, bitscore = _split_scored_entry(scored_entry)
    genome_id, contig_id, _gene_id = parse_savedname(protein_id)
    rf_distance = rf / maxrf if maxrf else 0.0
    return {
        "genome": genome,
        "protein_id": protein_id,
        "contig_id": contig_id,
        "contig_marker_support": int(contig_support_map.get(protein_id, 1)),
        "bitscore": bitscore,
        "rf_distance": rf_distance,
        "informative_splits": maxrf,
        "scored_entry": scored_entry,
    }


def _initial_selected_entries(
    dups: dict[str, list[str]],
    preferred_proteins: dict[str, str] | None = None,
) -> dict[str, str]:
    preferred_proteins = preferred_proteins or {}
    selected: dict[str, str] = {}
    for genome, entries in dups.items():
        if not entries:
            continue
        preferred_id = preferred_proteins.get(genome)
        if preferred_id is not None:
            matching = [
                entry for entry in entries
                if entry.rsplit(":", 1)[0] == preferred_id
            ]
            if matching:
                selected[genome] = matching[0]
                continue
        selected[genome] = _best_score(entries)
    return selected


def _optimize_selected_entries(
    marker_tree: Tree,
    species_tree: Tree,
    dups: dict[str, list[str]],
    contig_support_map: dict[str, int],
    selection_mode: str,
    max_rounds: int,
    locked_genomes: set[str],
    preferred_proteins: dict[str, str] | None = None,
) -> dict[str, str]:
    selected = _initial_selected_entries(dups, preferred_proteins=preferred_proteins)

    if selection_mode == "legacy":
        fixed_selected = dict(selected)
        for genome in sorted(dups):
            entries = dups[genome]
            if len(entries) <= 1 or genome in locked_genomes:
                continue
            candidates = [
                _evaluate_candidate(marker_tree, species_tree, dups, fixed_selected, genome, entry, contig_support_map)
                for entry in entries
            ]
            selected[genome] = choose_best_candidate(
                candidates,
                preferred_protein_id=fixed_selected.get(genome, "").rsplit(":", 1)[0] if fixed_selected.get(genome) else None,
            )["scored_entry"]
        return selected

    if selection_mode != "coordinate":
        raise ValueError(f"Unknown selection mode: {selection_mode}")

    rounds = max(1, max_rounds)
    for _ in range(rounds):
        changed = False
        for genome in sorted(dups):
            entries = dups[genome]
            if len(entries) <= 1 or genome in locked_genomes:
                continue
            candidates = [
                _evaluate_candidate(marker_tree, species_tree, dups, selected, genome, entry, contig_support_map)
                for entry in entries
            ]
            current_selected = selected.get(genome)
            best = choose_best_candidate(
                candidates,
                preferred_protein_id=current_selected.rsplit(":", 1)[0] if current_selected else None,
            )["scored_entry"]
            if best != selected.get(genome):
                selected[genome] = best
                changed = True
        if not changed:
            break
    return selected


def resolve_marker_tree(
    marker_tree_path: str,
    species_tree_path: str,
    table_path: str,
    marker_name: str,
    ls_refs: list[str] | None,
    selection_mode: str,
    max_rounds: int,
    lock_references: bool,
    initial_kept: dict[tuple[str, str], str] | None = None,
) -> tuple[list[str], list[dict]]:
    score_table, score_col = _cached_score_table(table_path)
    marker_tree = Tree(marker_tree_path)
    lst_nodes = [leaf.name for leaf in marker_tree.iter_leaves()]
    dups = _build_duplicate_map(lst_nodes, score_table, score_col)
    contig_support_map = _build_contig_support_map(score_table)
    species_tree = _cached_species_tree(species_tree_path)

    locked_genomes = set()
    if lock_references and ls_refs is not None:
        locked_genomes = {
            genome for genome in dups if f"{genome}.faa" in set(ls_refs)
        }

    preferred_proteins = {}
    if initial_kept is not None:
        preferred_proteins = {
            genome: protein_id
            for (marker, genome), protein_id in initial_kept.items()
            if marker == marker_name
        }

    selected = _optimize_selected_entries(
        marker_tree=marker_tree,
        species_tree=species_tree,
        dups=dups,
        contig_support_map=contig_support_map,
        selection_mode=selection_mode,
        max_rounds=max_rounds,
        locked_genomes=locked_genomes,
        preferred_proteins=preferred_proteins,
    )

    records: list[dict] = []
    bad_nodes: list[str] = []
    for genome in sorted(dups):
        entries = dups[genome]
        if len(entries) <= 1:
            continue
        candidates = [
            _evaluate_candidate(marker_tree, species_tree, dups, selected, genome, entry, contig_support_map)
            for entry in entries
        ]
        kept = selected[genome]
        for candidate in candidates:
            records.append(
                {
                    "genome": genome,
                    "protein_id": candidate["protein_id"],
                    "marker": marker_name,
                    "rf_distance": candidate["rf_distance"],
                    "informative_splits": candidate["informative_splits"],
                    "bitscore": candidate["bitscore"],
                    "status": "Kept" if candidate["scored_entry"] == kept else "Removed",
                }
            )
        bad_nodes.extend(
            entry.rsplit(":", 1)[0].replace("/", "|")
            for entry in entries
            if entry != kept
        )

    cleaned_nodes = [node for node in lst_nodes if node not in set(bad_nodes)]
    return cleaned_nodes, records


def _load_kept_assignments(rf_outfile: str) -> dict[tuple[str, str], str]:
    kept: dict[tuple[str, str], str] = {}
    if not os.path.exists(rf_outfile):
        return kept
    with open(rf_outfile) as handle:
        next(handle, None)
        for line in handle:
            parts = line.strip().split()
            if len(parts) != 4:
                continue
            protein_id, marker, _rf_distance, status = parts
            if status == "Kept":
                kept[(marker, protein_id.split("/")[0])] = protein_id
    return kept


def _process_tree_worker(args):
    """Worker: RF-distance based duplicate resolution for one marker tree.

    Returns ``records`` for the marker. The parent process aggregates and
    writes ``marker_selection_rf_values.txt`` once to avoid the concurrent
    append race the legacy code had when many workers appended to the
    shared file in parallel (nondeterministic ordering; potential write
    interleaving for payloads > ``PIPE_BUF`` on Linux).
    """
    (
        filepath,
        table_path,
        species_tree_path,
        outdir,
        ls_refs,
        selection_mode,
        max_rounds,
        lock_references,
        initial_kept,
    ) = args
    marker_name = os.path.basename(filepath).split(".")[0]
    cleaned_nodes, records = resolve_marker_tree(
        marker_tree_path=filepath,
        species_tree_path=species_tree_path,
        table_path=table_path,
        marker_name=marker_name,
        ls_refs=ls_refs,
        selection_mode=selection_mode,
        max_rounds=max_rounds,
        lock_references=lock_references,
        initial_kept=initial_kept,
    )

    with open(os.path.join(outdir, "removed", marker_name), "w") as f:
        removed = []
        for record in records:
            if record["status"] == "Removed":
                removed.append(record["protein_id"].replace("/", "|"))
        for item in removed:
            f.write(f"{item}\n")
        f.write(f"{len(removed)} {len(cleaned_nodes) + len(removed)}\n{'*' * 80}\n")

    t = Tree(filepath, format=1)
    t_final = t.copy()
    t_final.prune(cleaned_nodes)
    t_final.write(
        format=1,
        outfile=os.path.join(
            outdir, "protTrees", "no_duplicates", "out",
            f"_no_dups_{marker_name}_.nw",
        ),
    )
    return records


def run_noperm(
    cfg: Config,
    ls_refs: list[str] | None,
    species_tree_path: str | None = None,
    initial_kept: dict[tuple[str, str], str] | None = None,
):
    """RF-distance based marker selection.

    For each genome with duplicate hits in a marker, prune tree with each alternative,
    compare RF distance to species tree, keep the copy producing lowest RF distance.
    """
    treeout_dir = os.path.join(cfg.outdir, "treeouts_protTrees")
    table_path = os.path.join(cfg.outdir, "table_elim_dups")

    # create output directories
    for d in [
        os.path.join(cfg.outdir, "protTrees"),
        os.path.join(cfg.outdir, "protTrees", "no_duplicates"),
        os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out"),
        os.path.join(cfg.outdir, "protTrees", "no_singles"),
        os.path.join(cfg.outdir, "removed"),
    ]:
        os.makedirs(d, exist_ok=True)

    rf_outfile = os.path.join(cfg.outdir, "marker_selection_rf_values.txt")

    if species_tree_path is None:
        species_tree_path = os.path.join(cfg.outdir, "tree.nwk")
    ls_of_files = glob.glob(os.path.join(treeout_dir, "*"))
    args = [
        (
            f,
            table_path,
            species_tree_path,
            cfg.outdir,
            ls_refs,
            getattr(cfg, "selection_mode", "coordinate"),
            getattr(cfg, "selection_max_rounds", 5),
            getattr(cfg, "lock_references", False),
            initial_kept,
        )
        for f in ls_of_files
    ]

    worker_results = map_processed(_process_tree_worker, args, cfg.num_cpus)

    _write_rf_values_file(rf_outfile, worker_results)
    return _load_kept_assignments(rf_outfile)


def _write_rf_values_file(path: str, worker_results: list[list[dict]] | None) -> None:
    """Write ``marker_selection_rf_values.txt`` deterministically.

    Flattens per-marker record lists, sorts by ``(marker, protein_id)``, and
    writes once from the parent process. Replaces the legacy concurrent
    append-per-worker pattern that produced nondeterministic ordering and
    could interleave writes on Linux when records exceeded ``PIPE_BUF``.
    """
    worker_results = worker_results or []
    all_records = [rec for records in worker_results if records for rec in records]
    all_records.sort(key=lambda r: (str(r["marker"]), str(r["protein_id"])))
    with open(path, "w") as f:
        f.write("ProteinID MarkerGene RFdistance Status\n")
        for record in all_records:
            f.write(
                f"{record['protein_id']} {record['marker']} "
                f"{record['rf_distance']:.6f} {record['status']}\n"
            )


def _tree_to_genome_level(tree: Tree) -> Tree:
    reduced = tree.copy()
    seen_genomes = set()
    keep_leaves = []
    for leaf in reduced.iter_leaves():
        genome = leaf.name.split("|")[0]
        if genome not in seen_genomes:
            seen_genomes.add(genome)
            keep_leaves.append(leaf.name)
    reduced.prune(keep_leaves)
    for leaf in reduced.iter_leaves():
        leaf.name = leaf.name.split("|")[0]
    return reduced


def _rf_distance_between(species_tree: Tree, marker_tree: Tree) -> float:
    species = species_tree.copy()
    marker = _tree_to_genome_level(marker_tree)
    species.prune([leaf.name for leaf in marker.iter_leaves()])
    rf, maxrf, *_ = species.robinson_foulds(marker, unrooted_trees=True)
    return rf / maxrf if maxrf else 0.0


def choose_tree_by_rf(species_tree: Tree, original_tree: Tree, candidate_tree: Tree) -> Tree:
    """Accept singleton pruning only when it improves marker/species RF."""
    before = _rf_distance_between(species_tree, original_tree)
    after = _rf_distance_between(species_tree, candidate_tree)
    return candidate_tree if after < before else original_tree


def _nearest_genome_neighbors(tree: Tree, leaf_name: str, k: int) -> list[str]:
    ordered, _distances = _nearest_genome_neighbor_distances(tree, leaf_name, k)
    return ordered


def _nearest_genome_neighbor_distances(
    tree: Tree,
    leaf_name: str,
    k: int,
) -> tuple[list[str], dict[str, float]]:
    leaf = next(leaf for leaf in tree.iter_leaves() if leaf.name == leaf_name)
    neighbors = []
    for other in tree.iter_leaves():
        if other.name == leaf_name:
            continue
        neighbors.append((float(leaf.get_distance(other)), other.name.split("|")[0]))
    neighbors.sort(key=lambda item: (item[0], item[1]))

    ordered: list[str] = []
    distances: dict[str, float] = {}
    for dist, genome in neighbors:
        if genome in distances:
            continue
        ordered.append(genome)
        distances[genome] = dist
        if len(ordered) >= k:
            break
    return ordered, distances


def _leaf_rf_delta(species_tree: Tree, working_tree: Tree, leaf_name: str) -> tuple[float, Tree]:
    before = _rf_distance_between(species_tree, working_tree)
    remaining = [leaf.name for leaf in working_tree.iter_leaves() if leaf.name != leaf_name]
    candidate = working_tree.copy()
    candidate.prune(remaining)
    after = _rf_distance_between(species_tree, candidate)
    return before - after, candidate


def _branch_length_outlier(tree: Tree, leaf_name: str) -> float:
    leaves = list(tree.iter_leaves())
    distances = [float(leaf.dist) for leaf in leaves]
    if len(distances) < 3:
        return 0.0
    mean = sum(distances) / len(distances)
    var = sum((dist - mean) ** 2 for dist in distances) / len(distances)
    stdev = var ** 0.5
    if stdev == 0:
        return 0.0
    target = next(leaf for leaf in leaves if leaf.name == leaf_name)
    z = (float(target.dist) - mean) / stdev
    return max(0.0, z)


def _leaf_bitscore_outlier(
    tree: Tree,
    leaf_name: str,
    score_table: pd.DataFrame | None,
    score_col: str | None,
) -> float:
    if score_table is None or score_col is None:
        return 0.0
    leaves = [leaf.name for leaf in tree.iter_leaves()]
    scores = []
    for name in leaves:
        row = score_table.loc[name.replace("|", "/")]
        scores.append(float(row[score_col]))
    if len(scores) < 3:
        return 0.0
    mean = sum(scores) / len(scores)
    var = sum((score - mean) ** 2 for score in scores) / len(scores)
    stdev = var ** 0.5
    if stdev == 0:
        return 0.0
    target = float(score_table.loc[leaf_name.replace("|", "/")][score_col])
    z = (mean - target) / stdev
    return max(0.0, z)


def _leaf_neighbor_overlap(species_tree: Tree, working_tree: Tree, leaf_name: str, k: int) -> float:
    species_neighbors = _nearest_genome_neighbors(species_tree, leaf_name.split("|")[0], k)
    gene_neighbors = _nearest_genome_neighbors(working_tree, leaf_name, k)
    if not species_neighbors:
        return 1.0
    overlap = len(set(species_neighbors) & set(gene_neighbors))
    return overlap / len(species_neighbors)


def _canonical_singleton_mode(mode: str) -> str:
    return LEGACY_SINGLETON_MODE_ALIASES.get(mode, mode)


def singleton_mode_uses_global_rf_gate(mode: str) -> bool:
    """Return whether automatic pruning uses the final RF safeguard.

    Retained for callers of the previous policy helper. Every supported mode
    now uses the safeguard.
    """
    return True


def _nearest_genome_neighbor_profile(
    tree: Tree,
    leaf_name: str,
    k: int,
) -> tuple[list[str], dict[str, float]]:
    ordered, distances = _nearest_genome_neighbor_distances(tree, leaf_name, k)
    max_dist = max(distances.values(), default=0.0)
    if max_dist > 0:
        distances = {
            genome: dist / max_dist
            for genome, dist in distances.items()
        }
    return ordered, distances


def _leaf_topoknn_score(species_tree: Tree, working_tree: Tree, leaf_name: str, k: int) -> float:
    species_neighbors, species_dist = _nearest_genome_neighbor_profile(
        species_tree,
        leaf_name.split("|")[0],
        k,
    )
    gene_neighbors, gene_dist = _nearest_genome_neighbor_profile(working_tree, leaf_name, k)
    if not species_neighbors:
        return 0.0

    overlap = len(set(species_neighbors) & set(gene_neighbors)) / len(species_neighbors)
    species_pos = {genome: idx for idx, genome in enumerate(species_neighbors)}
    gene_pos = {genome: idx for idx, genome in enumerate(gene_neighbors)}
    max_rank = max(1, len(species_neighbors) - 1)

    rank_penalty = 0.0
    distance_penalty = 0.0
    for genome in species_neighbors:
        if genome in gene_pos:
            rank_penalty += abs(species_pos[genome] - gene_pos[genome]) / max_rank
        else:
            rank_penalty += 1.0
        distance_penalty += abs(species_dist.get(genome, 1.0) - gene_dist.get(genome, 1.0))

    rank_penalty /= len(species_neighbors)
    distance_penalty /= len(species_neighbors)
    return (1.25 * (1.0 - overlap)) + rank_penalty + (0.5 * distance_penalty)


def _lca_genome_purity(tree: Tree, leaf_names: list[str]) -> float:
    unique_leaf_names = list(dict.fromkeys(str(name) for name in leaf_names if str(name)))
    if not unique_leaf_names:
        return 0.0
    if len(unique_leaf_names) == 1:
        return 1.0
    clade = tree.get_common_ancestor(unique_leaf_names)
    subtree_genomes = {
        leaf.name.split("|")[0]
        for leaf in clade.iter_leaves()
    }
    if not subtree_genomes:
        return 0.0
    target_genomes = {name.split("|")[0] for name in unique_leaf_names}
    return len(target_genomes & subtree_genomes) / len(subtree_genomes)


def _build_species_anchor_context(
    species_tree: Tree,
    k: int,
) -> dict[str, dict[str, float | int | list[str]]]:
    leaves = [leaf.name for leaf in species_tree.iter_leaves()]
    raw_neighbors: dict[str, tuple[list[str], dict[str, float]]] = {}
    nearest_neighbor_distances: list[float] = []
    for leaf_name in leaves:
        ordered, distances = _nearest_genome_neighbor_distances(species_tree, leaf_name, k)
        raw_neighbors[leaf_name] = (ordered, distances)
        if distances:
            nearest_neighbor_distances.append(min(distances.values()))

    baseline_distance = statistics.median(nearest_neighbor_distances) if nearest_neighbor_distances else 0.0
    context: dict[str, dict[str, float | int | list[str]]] = {}
    for leaf_name, (ordered, distances) in raw_neighbors.items():
        compactness = statistics.median(distances.values()) if distances else 0.0
        if compactness > 0.0 and baseline_distance > 0.0:
            compactness_score = 1.0 / (1.0 + (compactness / baseline_distance))
        elif compactness > 0.0:
            compactness_score = 0.5
        else:
            compactness_score = 1.0
        first_distance = min(distances.values(), default=0.0)
        long_branch_z = _robust_zscore(first_distance, nearest_neighbor_distances)
        long_branch_support = 1.0 / (1.0 + long_branch_z)
        species_purity = _lca_genome_purity(species_tree, [leaf_name, *ordered])
        count_support = min(1.0, len(ordered) / max(1, k))
        anchor_score = max(
            0.0,
            min(
                1.0,
                (0.5 * species_purity)
                + (0.25 * compactness_score)
                + (0.15 * long_branch_support)
                + (0.10 * count_support),
            ),
        )
        context[leaf_name] = {
            "neighbors": ordered,
            "target_neighbor_count": len(ordered),
            "species_anchor_purity": species_purity,
            "species_anchor_compactness": compactness,
            "species_anchor_compactness_score": compactness_score,
            "species_long_branch_z": long_branch_z,
            "species_long_branch_support": long_branch_support,
            "species_anchor_score": anchor_score,
        }
    return context


def _leaf_neighbor_clade_metrics(
    species_tree: Tree,
    working_tree: Tree,
    leaf_name: str,
    *,
    k: int,
    recipient_consensus_score: float,
    species_anchor_context: dict[str, dict[str, float | int | list[str]]],
    genome_to_leaf: dict[str, str],
    leaf_lookup: dict[str, object],
) -> dict[str, float]:
    recipient_genome = leaf_name.split("|")[0]
    anchor = species_anchor_context.get(recipient_genome, {})
    anchor_neighbors = [str(genome) for genome in anchor.get("neighbors", [])]
    target_neighbor_count = int(anchor.get("target_neighbor_count", len(anchor_neighbors)))
    present_neighbors = [
        genome
        for genome in anchor_neighbors
        if genome in genome_to_leaf and genome != recipient_genome
    ]
    present_neighbor_count = len(present_neighbors)
    present_neighbor_fraction = (
        present_neighbor_count / max(1, target_neighbor_count)
        if target_neighbor_count > 0
        else 0.0
    )

    neighbor_anchor_purity = 0.0
    join_purity = 0.0
    purity_drop = 0.0
    anchor_knn_agreement = 0.0
    attachment_gap = 0.0
    if present_neighbors:
        anchor_leaf_names = [genome_to_leaf[genome] for genome in present_neighbors]
        neighbor_anchor_purity = _lca_genome_purity(working_tree, anchor_leaf_names)
        join_purity = _lca_genome_purity(working_tree, anchor_leaf_names + [leaf_name])
        purity_drop = max(0.0, neighbor_anchor_purity - join_purity)
        gene_neighbors = _nearest_genome_neighbors(working_tree, leaf_name, max(1, present_neighbor_count))
        anchor_knn_agreement = len(set(gene_neighbors) & set(present_neighbors)) / present_neighbor_count

        anchor_clade = working_tree.get_common_ancestor(anchor_leaf_names)
        anchor_distances = [
            float(leaf_lookup[name].get_distance(anchor_clade))
            for name in anchor_leaf_names
        ]
        candidate_distance = float(leaf_lookup[leaf_name].get_distance(anchor_clade))
        attachment_gap = max(0.0, candidate_distance - statistics.median(anchor_distances))

    knn_disagreement = 1.0 - anchor_knn_agreement if present_neighbor_count else 0.0
    recipient_support = min(max(float(recipient_consensus_score), 0.0), 3.0) / 3.0
    neighbor_clade_score = (
        1.25 * float(anchor.get("species_anchor_score", 0.0))
        + 1.25 * neighbor_anchor_purity
        + 2.0 * purity_drop
        + 1.0 * knn_disagreement
        + 0.5 * min(attachment_gap, 2.0)
        + 0.25 * recipient_support
    )

    return {
        "target_neighbor_count": target_neighbor_count,
        "present_neighbor_count": present_neighbor_count,
        "present_neighbor_fraction": present_neighbor_fraction,
        "species_anchor_purity": float(anchor.get("species_anchor_purity", 0.0)),
        "species_anchor_compactness": float(anchor.get("species_anchor_compactness", 0.0)),
        "species_anchor_compactness_score": float(anchor.get("species_anchor_compactness_score", 0.0)),
        "species_long_branch_z": float(anchor.get("species_long_branch_z", 0.0)),
        "species_long_branch_support": float(anchor.get("species_long_branch_support", 0.0)),
        "species_anchor_score": float(anchor.get("species_anchor_score", 0.0)),
        "neighbor_anchor_purity": neighbor_anchor_purity,
        "join_purity": join_purity,
        "purity_drop": purity_drop,
        "anchor_knn_agreement": anchor_knn_agreement,
        "attachment_gap": attachment_gap,
        "neighbor_clade_score": neighbor_clade_score,
    }


def _load_alignment_sequence_map(
    alignment_path: str | None,
    *,
    keep_ids: set[str] | None = None,
) -> dict[str, str]:
    if alignment_path is None or not os.path.exists(alignment_path):
        return {}
    sequences: dict[str, str] = {}
    with open(alignment_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if keep_ids is not None and record.id not in keep_ids:
                continue
            sequences[str(record.id)] = str(record.seq).upper()
    return sequences


def _gap_aware_hamming_distance(sequence: str, consensus: str) -> float:
    if not sequence or not consensus:
        return 0.0
    width = min(len(sequence), len(consensus))
    informative = 0
    mismatches = 0
    for index in range(width):
        left = sequence[index]
        right = consensus[index]
        if left == "-" and right == "-":
            continue
        informative += 1
        if left != right:
            mismatches += 1
    return (mismatches / informative) if informative else 0.0


def _alignment_consensus(sequences: list[str]) -> str:
    if not sequences:
        return ""
    width = min(len(sequence) for sequence in sequences)
    columns: list[str] = []
    for index in range(width):
        counts: dict[str, int] = {}
        for sequence in sequences:
            residue = sequence[index]
            if residue == "-":
                continue
            counts[residue] = counts.get(residue, 0) + 1
        if not counts:
            columns.append("-")
            continue
        columns.append(
            min(
                counts.items(),
                key=lambda item: (-item[1], item[0]),
            )[0]
        )
    return "".join(columns)


def _robust_zscore(value: float, background: list[float]) -> float:
    if len(background) < 3:
        return 0.0
    median = statistics.median(background)
    deviations = [abs(item - median) for item in background]
    mad = statistics.median(deviations)
    if mad > 0:
        return max(0.0, 0.6745 * (value - median) / mad)
    stdev = statistics.pstdev(background)
    if stdev == 0:
        return 0.0
    return max(0.0, (value - median) / stdev)


def _leaf_recipient_consensus_score(
    species_tree: Tree,
    working_tree: Tree,
    leaf_name: str,
    k: int,
    alignment_sequences: dict[str, str],
) -> float:
    if not alignment_sequences or leaf_name not in alignment_sequences:
        return 0.0

    genome_to_leaf = {
        leaf.name.split("|")[0]: leaf.name
        for leaf in working_tree.iter_leaves()
        if leaf.name in alignment_sequences
    }
    recipient_genome = leaf_name.split("|")[0]
    neighbor_genomes = _nearest_genome_neighbors(species_tree, recipient_genome, k)
    neighbor_sequences = [
        alignment_sequences[genome_to_leaf[genome]]
        for genome in neighbor_genomes
        if genome in genome_to_leaf
    ]
    if len(neighbor_sequences) < 3:
        return 0.0

    consensus = _alignment_consensus(neighbor_sequences)
    candidate_distance = _gap_aware_hamming_distance(alignment_sequences[leaf_name], consensus)
    neighbor_distances = [
        _gap_aware_hamming_distance(sequence, consensus)
        for sequence in neighbor_sequences
    ]
    return _robust_zscore(candidate_distance, neighbor_distances)


def _score_singleton_candidates(
    species_tree: Tree,
    working_tree: Tree,
    *,
    k: int,
    score_table: pd.DataFrame | None = None,
    score_col: str | None = None,
    alignment_path: str | None = None,
) -> list[dict]:
    working_leaves = list(working_tree.iter_leaves())
    keep_ids = {leaf.name for leaf in working_leaves}
    alignment_sequences = _load_alignment_sequence_map(
        alignment_path,
        keep_ids=keep_ids,
    )
    species_anchor_context = _build_species_anchor_context(species_tree, k)
    genome_to_leaf = {
        leaf.name.split("|")[0]: leaf.name
        for leaf in working_leaves
    }
    leaf_lookup = {
        leaf.name: leaf
        for leaf in working_leaves
    }
    candidates: list[dict] = []
    for leaf in working_leaves:
        delta_rf, candidate = _leaf_rf_delta(species_tree, working_tree, leaf.name)
        overlap = _leaf_neighbor_overlap(species_tree, working_tree, leaf.name, k)
        branch_outlier = _branch_length_outlier(working_tree, leaf.name)
        bitscore_outlier = _leaf_bitscore_outlier(working_tree, leaf.name, score_table, score_col)
        recipient_consensus_score = _leaf_recipient_consensus_score(
            species_tree,
            working_tree,
            leaf.name,
            k,
            alignment_sequences,
        )
        genome_id, contig_id, gene_id = parse_sequence_id(leaf.name)
        topoknn_score = (
            _leaf_topoknn_score(species_tree, working_tree, leaf.name, k)
            + (0.25 * branch_outlier)
            + (0.25 * bitscore_outlier)
        )
        neighbor_clade_metrics = _leaf_neighbor_clade_metrics(
            species_tree,
            working_tree,
            leaf.name,
            k=k,
            recipient_consensus_score=recipient_consensus_score,
            species_anchor_context=species_anchor_context,
            genome_to_leaf=genome_to_leaf,
            leaf_lookup=leaf_lookup,
        )
        candidates.append(
            {
                "leaf_name": leaf.name,
                "taxa_count": len(keep_ids),
                "delta_rf": delta_rf,
                "neighbor_overlap": overlap,
                "recipient_neighbor_overlap": overlap,
                "topoknn_score": topoknn_score,
                "branch_outlier": branch_outlier,
                "bitscore_outlier": bitscore_outlier,
                "recipient_consensus_score": recipient_consensus_score,
                "genome": genome_id,
                "contig_id": contig_id,
                "gene_id": gene_id,
                **neighbor_clade_metrics,
                "candidate_tree": candidate,
            }
        )
    return candidates


def _normalize_candidate_metric(candidates: list[dict], key: str) -> dict[str, float]:
    values = [float(candidate[key]) for candidate in candidates]
    if not values:
        return {}
    low = min(values)
    high = max(values)
    if high == low:
        normalized = 1.0 if high > 0 else 0.0
        return {candidate["leaf_name"]: normalized for candidate in candidates}
    return {
        candidate["leaf_name"]: (float(candidate[key]) - low) / (high - low)
        for candidate in candidates
    }


def _normalize_candidate_values(
    candidates: list[dict],
    value_fn,
) -> dict[str, float]:
    values = {str(candidate["leaf_name"]): float(value_fn(candidate)) for candidate in candidates}
    if not values:
        return {}
    low = min(values.values())
    high = max(values.values())
    if high == low:
        normalized = 1.0 if high > 0 else 0.0
        return {leaf_name: normalized for leaf_name in values}
    return {
        leaf_name: (value - low) / (high - low)
        for leaf_name, value in values.items()
    }


def _finalize_singleton_choice(candidate: dict, *, score_key: str) -> dict:
    chosen = dict(candidate)
    chosen["genome"] = chosen["leaf_name"].split("|")[0]
    chosen["score"] = float(chosen.get(score_key, 0.0))
    return chosen


def _best_and_runner_up(candidates: list[dict], score_key: str) -> tuple[dict | None, dict | None]:
    ranked = sorted(
        candidates,
        key=lambda candidate: (
            -float(candidate.get(score_key, 0.0)),
            -float(candidate.get("delta_rf", 0.0)),
            -float(candidate.get("topoknn_score", 0.0)),
            str(candidate.get("leaf_name", "")),
        ),
    )
    best = ranked[0] if ranked else None
    runner_up = ranked[1] if len(ranked) > 1 else None
    return best, runner_up


def _choose_composite_candidate(candidates: list[dict]) -> dict | None:
    if not candidates:
        return None
    rf_norm = _normalize_candidate_metric(candidates, "delta_rf")
    topoknn_norm = _normalize_candidate_metric(candidates, "topoknn_score")
    branch_norm = _normalize_candidate_metric(candidates, "branch_outlier")
    bitscore_norm = _normalize_candidate_metric(candidates, "bitscore_outlier")

    for candidate in candidates:
        leaf_name = candidate["leaf_name"]
        candidate["composite_score"] = (
            1.5 * rf_norm.get(leaf_name, 0.0)
            + 1.0 * topoknn_norm.get(leaf_name, 0.0)
            + 0.5 * branch_norm.get(leaf_name, 0.0)
            + 0.5 * bitscore_norm.get(leaf_name, 0.0)
        )

    best, runner_up = _best_and_runner_up(candidates, "composite_score")
    if best is None:
        return None
    runner_up_score = float(runner_up.get("composite_score", 0.0)) if runner_up is not None else 0.0
    best["composite_score_gap"] = float(best.get("composite_score", 0.0)) - runner_up_score
    if float(best.get("composite_score", 0.0)) < COMPOSITE_SCORE_THRESHOLD:
        return None
    if runner_up is not None and float(best["composite_score_gap"]) < COMPOSITE_SCORE_MARGIN:
        return None
    if (
        float(best.get("delta_rf", 0.0)) < 0.05
        and float(best.get("topoknn_score", 0.0)) < 0.75
        and float(best.get("bitscore_outlier", 0.0)) < 1.0
        and float(best.get("branch_outlier", 0.0)) < 1.0
    ):
        return None
    return best


def _composite_candidate_is_eligible(candidate: dict) -> bool:
    return (
        float(candidate.get("composite_score", 0.0)) >= COMPOSITE_SCORE_THRESHOLD
        and (
            float(candidate.get("delta_rf", 0.0)) >= 0.05
            or float(candidate.get("topoknn_score", 0.0)) >= 0.75
            or float(candidate.get("bitscore_outlier", 0.0)) >= 1.0
            or float(candidate.get("branch_outlier", 0.0)) >= 1.0
        )
    )


def _recipient_ranked_candidates(candidates: list[dict]) -> list[dict]:
    if not candidates:
        return []
    calibration = _recipient_consensus_calibration(len(candidates))
    eligible = [
        candidate
        for candidate in candidates
        if _passes_recipient_raw_floor(candidate, calibration)
    ]
    if not eligible:
        return []

    rf_norm = _normalize_candidate_metric(candidates, "delta_rf")
    topoknn_norm = _normalize_candidate_values(
        eligible,
        lambda candidate: min(float(candidate.get("topoknn_score", 0.0)), float(calibration["topoknn_cap"])),
    )
    recipient_norm = _normalize_candidate_values(
        eligible,
        lambda candidate: min(
            float(candidate.get("recipient_consensus_score", 0.0)),
            float(calibration["recipient_cap"]),
        ),
    )

    rf_weight = float(calibration["delta_weight"])
    topoknn_weight = float(calibration["topoknn_weight"])
    recipient_weight = float(calibration["recipient_weight"])
    for candidate in eligible:
        leaf_name = candidate["leaf_name"]
        candidate["recipient_rank_score"] = (
            rf_weight * rf_norm.get(leaf_name, 0.0)
            + topoknn_weight * topoknn_norm.get(leaf_name, 0.0)
            + recipient_weight * recipient_norm.get(leaf_name, 0.0)
        )
    return eligible


def _choose_ranked_recipient_candidate(candidates: list[dict]) -> dict | None:
    eligible = _recipient_ranked_candidates(candidates)
    if not eligible:
        return None

    calibration = _recipient_consensus_calibration(len(candidates))
    best, runner_up = _best_and_runner_up(eligible, "recipient_rank_score")
    if best is None:
        return None
    runner_up_score = float(runner_up.get("recipient_rank_score", 0.0)) if runner_up is not None else 0.0
    best["recipient_rank_gap"] = float(best.get("recipient_rank_score", 0.0)) - runner_up_score
    if runner_up is not None and float(best["recipient_rank_gap"]) < float(calibration["rank_margin"]):
        return None
    return best


def _choose_neighbor_clade_candidate(candidates: list[dict]) -> dict | None:
    if not candidates:
        return None
    eligible = [
        candidate
        for candidate in candidates
        if int(candidate.get("present_neighbor_count", 0)) >= 2
        and float(candidate.get("species_anchor_score", 0.0)) >= 0.25
    ]
    if not eligible:
        return None
    best, runner_up = _best_and_runner_up(eligible, "neighbor_clade_score")
    if best is None:
        return None
    runner_up_score = float(runner_up.get("neighbor_clade_score", 0.0)) if runner_up is not None else 0.0
    best["neighbor_clade_score_gap"] = float(best.get("neighbor_clade_score", 0.0)) - runner_up_score
    if float(best.get("neighbor_clade_score", 0.0)) < NEIGHBOR_CLADE_MIN_SCORE:
        return None
    return best


def choose_singleton_prune(
    species_tree: Tree,
    working_tree: Tree,
    *,
    mode: str,
    k: int,
    score_table: pd.DataFrame | None = None,
    score_col: str | None = None,
    alignment_path: str | None = None,
    scored_candidates: list[dict] | None = None,
) -> dict | None:
    mode = _canonical_singleton_mode(mode)
    if len(list(working_tree.iter_leaves())) <= 4:
        return None
    if scored_candidates is None:
        candidates = _score_singleton_candidates(
            species_tree,
            working_tree,
            k=k,
            score_table=score_table,
            score_col=score_col,
            alignment_path=alignment_path,
        )
    else:
        candidates = [dict(candidate) for candidate in scored_candidates]
    if not candidates:
        return None

    positive_rf = [candidate for candidate in candidates if float(candidate["delta_rf"]) > 0]
    if mode in {"delta_rf", "topoknn", "outlier", "hybrid"} and not positive_rf:
        return None

    if mode == "delta_rf":
        best = max(
            positive_rf,
            key=lambda candidate: (
                float(candidate["delta_rf"]),
                float(candidate["topoknn_score"]),
                str(candidate["leaf_name"]),
            ),
        )
        return _finalize_singleton_choice(best, score_key="delta_rf")

    if mode in {"topoknn", "outlier"}:
        best = max(
            positive_rf,
            key=lambda candidate: (
                float(candidate["topoknn_score"]),
                float(candidate["delta_rf"]),
                str(candidate["leaf_name"]),
            ),
        )
        if float(best["topoknn_score"]) < 0.75:
            return None
        return _finalize_singleton_choice(best, score_key="topoknn_score")

    if mode == "hybrid":
        rf_norm = _normalize_candidate_metric(positive_rf, "delta_rf")
        topoknn_norm = _normalize_candidate_metric(positive_rf, "topoknn_score")
        branch_norm = _normalize_candidate_metric(positive_rf, "branch_outlier")
        bitscore_norm = _normalize_candidate_metric(positive_rf, "bitscore_outlier")

        for candidate in positive_rf:
            leaf_name = candidate["leaf_name"]
            candidate["hybrid_score"] = (
                1.5 * rf_norm.get(leaf_name, 0.0)
                + 1.5 * topoknn_norm.get(leaf_name, 0.0)
                + 0.25 * branch_norm.get(leaf_name, 0.0)
                + 0.25 * bitscore_norm.get(leaf_name, 0.0)
            )

        best_rf = max(
            positive_rf,
            key=lambda candidate: (
                float(candidate["delta_rf"]),
                float(candidate["topoknn_score"]),
                str(candidate["leaf_name"]),
            ),
        )
        best_topoknn = max(
            positive_rf,
            key=lambda candidate: (
                float(candidate["topoknn_score"]),
                float(candidate["delta_rf"]),
                str(candidate["leaf_name"]),
            ),
        )
        best = max(
            positive_rf,
            key=lambda candidate: (
                float(candidate["hybrid_score"]),
                float(candidate["delta_rf"]),
                float(candidate["topoknn_score"]),
                str(candidate["leaf_name"]),
            ),
        )
        if best["leaf_name"] == best_rf["leaf_name"] == best_topoknn["leaf_name"]:
            return _finalize_singleton_choice(best, score_key="hybrid_score")
        if float(best["delta_rf"]) >= 0.05 and float(best["topoknn_score"]) >= 1.0:
            return _finalize_singleton_choice(best, score_key="hybrid_score")
        return None

    if mode == "composite":
        best = _choose_composite_candidate(candidates)
        if best is None:
            return None
        return _finalize_singleton_choice(best, score_key="composite_score")

    if mode == "contig_consensus":
        best = _choose_composite_candidate(candidates)
        if best is None:
            return None
        return _finalize_singleton_choice(best, score_key="composite_score")

    if mode == "recipient_consensus":
        best = _choose_ranked_recipient_candidate(candidates)
        if best is None:
            return None
        return _finalize_singleton_choice(best, score_key="recipient_rank_score")

    if mode == "neighbor_clade":
        best = _choose_neighbor_clade_candidate(candidates)
        if best is None:
            return None
        return _finalize_singleton_choice(best, score_key="neighbor_clade_score")

    if mode == "gcp":
        # GCP does global cross-marker classification, not per-marker selection
        return None

    raise ValueError(f"Unknown singleton mode: {mode}")


def prune_singletons(
    species_tree: Tree,
    working_tree: Tree,
    *,
    mode: str,
    k: int,
    score_table: pd.DataFrame | None = None,
    score_col: str | None = None,
    alignment_path: str | None = None,
) -> Tree:
    chosen = choose_singleton_prune(
        species_tree,
        working_tree,
        mode=mode,
        k=k,
        score_table=score_table,
        score_col=score_col,
        alignment_path=alignment_path,
    )
    if chosen is None:
        return working_tree.copy()
    return chosen["candidate_tree"]


def _marker_name_from_tree_path(filepath: str) -> str:
    return os.path.basename(filepath).split(".")[0].split("_")[-2]


def _count_genome_marker_support(files: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for filepath in files:
        genomes = {
            leaf.name.split("|")[0]
            for leaf in Tree(filepath).iter_leaves()
        }
        for genome in genomes:
            counts[genome] = counts.get(genome, 0) + 1
    return counts


def _load_contig_marker_context(table_path: str) -> dict[tuple[str, str], set[str]]:
    if not os.path.exists(table_path):
        return {}
    df = pd.read_csv(table_path)
    if "savedname" not in df.columns or "namemodel" not in df.columns:
        return {}
    if "genome_id" not in df.columns or "contig_id" not in df.columns or "gene_id" not in df.columns:
        parsed = df["savedname"].apply(parse_savedname)
        df = df.copy()
        df["genome_id"] = parsed.apply(lambda item: item[0])
        df["contig_id"] = parsed.apply(lambda item: item[1])
        df["gene_id"] = parsed.apply(lambda item: item[2])
    df = df.copy()
    df["marker_id"] = df["namemodel"].astype(str).str.split("/").str[-1]
    context: dict[tuple[str, str], set[str]] = {}
    for row in df.itertuples(index=False):
        key = (str(row.genome_id), str(row.contig_id))
        context.setdefault(key, set()).add(str(row.marker_id))
    return context


def _load_contig_marker_hit_map(table_path: str) -> dict[tuple[str, str, str], str]:
    if not os.path.exists(table_path):
        return {}
    df = pd.read_csv(table_path)
    if "savedname" not in df.columns or "namemodel" not in df.columns:
        return {}
    parsed = df["savedname"].apply(parse_savedname)
    df = df.copy()
    df["genome_id"] = parsed.apply(lambda item: item[0])
    df["contig_id"] = parsed.apply(lambda item: item[1])
    df["marker_id"] = df["namemodel"].astype(str).str.split("/").str[-1]
    return {
        (str(row.genome_id), str(row.contig_id), str(row.marker_id)): str(row.savedname).replace("/", "|")
        for row in df.itertuples(index=False)
    }


def _build_marker_neighbor_context(
    files: list[str],
    *,
    species_tree_path: str,
    table_path: str,
    k: int,
    num_cpus: int = 1,
) -> dict[tuple[str, str, str], dict[str, float]]:
    hit_map = _load_contig_marker_hit_map(table_path)
    if not hit_map:
        return {}
    species_tree = _cached_species_tree(species_tree_path)
    marker_paths = {
        _marker_name_from_tree_path(filepath): filepath
        for filepath in files
    }
    hits_by_marker: dict[str, list[tuple[tuple[str, str, str], str]]] = {}
    for key, leaf_name in hit_map.items():
        _genome, _contig, marker = key
        hits_by_marker.setdefault(marker, []).append((key, leaf_name))

    def _process_marker(marker: str) -> list[tuple[tuple[str, str, str], dict[str, float]]]:
        filepath = marker_paths.get(marker)
        if filepath is None:
            return []
        tree = Tree(filepath)
        leaves = {leaf.name for leaf in tree.iter_leaves()}
        results: list[tuple[tuple[str, str, str], dict[str, float]]] = []
        for key, leaf_name in hits_by_marker[marker]:
            if leaf_name not in leaves:
                continue
            overlap = _leaf_neighbor_overlap(species_tree, tree, leaf_name, k)
            results.append((key, {"recipient_neighbor_overlap": overlap}))
        return results

    markers = sorted(hits_by_marker.keys())
    workers = max(1, min(int(num_cpus), len(markers)))
    marker_results = map_threaded(_process_marker, markers, workers)

    context: dict[tuple[str, str, str], dict[str, float]] = {}
    for per_marker in marker_results:
        for key, value in per_marker:
            context[key] = value
    return context


def _classify_singleton_proposals_legacy(
    proposals: list[dict],
    *,
    contig_marker_context: dict[tuple[str, str], set[str]],
) -> list[dict]:
    suspicious_by_contig: dict[tuple[str, str], set[str]] = {}
    for proposal in proposals:
        key = (str(proposal["genome"]), str(proposal.get("contig_id", "unknown_contig")))
        suspicious_by_contig.setdefault(key, set()).add(str(proposal["marker_name"]))

    classified: list[dict] = []
    for proposal in proposals:
        updated = dict(proposal)
        genome = str(updated["genome"])
        contig = str(updated.get("contig_id", "unknown_contig"))
        key = (genome, contig)
        markers_on_contig = set(contig_marker_context.get(key, set()))
        suspicious_markers = suspicious_by_contig.get(key, set())

        if contig == "unknown_contig":
            updated["singleton_class"] = "ambiguous"
        elif not markers_on_contig or markers_on_contig == {updated["marker_name"]}:
            updated["singleton_class"] = "contamination_candidate"
        elif any(marker not in suspicious_markers for marker in markers_on_contig - {updated["marker_name"]}):
            updated["singleton_class"] = "hgt_candidate"
        else:
            updated["singleton_class"] = "contamination_candidate"
        classified.append(updated)
    return classified


def _strong_unknown_contig_signal(proposal: dict) -> bool:
    calibration = _recipient_consensus_calibration(int(proposal.get("taxa_count", 0) or 0))
    return (
        float(proposal.get("delta_rf", 0.0)) >= float(calibration["unknown_delta_floor"])
        and float(proposal.get("topoknn_score", 0.0)) >= float(calibration["unknown_topoknn_floor"])
        and float(proposal.get("recipient_consensus_score", 0.0)) >= float(calibration["unknown_recipient_floor"])
    )


def _recipient_consensus_calibration(taxa_count: int) -> dict[str, float]:
    if taxa_count <= 60:
        return {
            "delta_floor": UNKNOWN_CONTIG_DELTA_FLOOR,
            "topoknn_floor": UNKNOWN_CONTIG_TOPOKNN_FLOOR,
            "recipient_floor": UNKNOWN_CONTIG_RECIPIENT_FLOOR,
            "delta_weight": 1.5,
            "topoknn_weight": 1.0,
            "recipient_weight": 0.75,
            "topoknn_cap": 3.5,
            "recipient_cap": 3.0,
            "rank_margin": RECIPIENT_CONSENSUS_RANK_MARGIN,
            "unknown_delta_floor": UNKNOWN_CONTIG_DELTA_FLOOR,
            "unknown_topoknn_floor": UNKNOWN_CONTIG_TOPOKNN_FLOOR,
            "unknown_recipient_floor": UNKNOWN_CONTIG_RECIPIENT_FLOOR,
        }
    return {
        "delta_floor": 0.01,
        "topoknn_floor": 0.8,
        "recipient_floor": 0.25,
        "delta_weight": 2.5,
        "topoknn_weight": 1.0,
        "recipient_weight": 0.25,
        "topoknn_cap": 3.0,
        "recipient_cap": 2.0,
        "rank_margin": 0.2,
        "unknown_delta_floor": 0.01,
        "unknown_topoknn_floor": 0.8,
        "unknown_recipient_floor": 0.25,
    }


def _passes_recipient_raw_floor(candidate: dict, calibration: dict[str, float]) -> bool:
    delta_rf = float(candidate.get("delta_rf", 0.0))
    topoknn = float(candidate.get("topoknn_score", 0.0))
    recipient = float(candidate.get("recipient_consensus_score", 0.0))
    return delta_rf >= float(calibration["delta_floor"]) and (
        topoknn >= float(calibration["topoknn_floor"])
        or recipient >= float(calibration["recipient_floor"])
    )


def _classify_neighbor_clade_proposals(proposals: list[dict]) -> list[dict]:
    classified: list[dict] = []
    for proposal in proposals:
        updated = dict(proposal)
        target_neighbor_count = int(updated.get("target_neighbor_count", 0))
        present_neighbor_count = int(updated.get("present_neighbor_count", 0))
        present_neighbor_fraction = float(updated.get("present_neighbor_fraction", 0.0))
        if target_neighbor_count > 0 and present_neighbor_fraction <= 0.0:
            present_neighbor_fraction = present_neighbor_count / target_neighbor_count

        required_neighbors = min(NEIGHBOR_CLADE_MIN_PRESENT, max(1, target_neighbor_count))
        support_ok = (
            present_neighbor_count >= required_neighbors
            and present_neighbor_fraction >= NEIGHBOR_CLADE_MIN_PRESENT_FRACTION
        )
        species_anchor_score = float(updated.get("species_anchor_score", 0.0))
        neighbor_anchor_purity = float(updated.get("neighbor_anchor_purity", 0.0))
        purity_drop = float(updated.get("purity_drop", 0.0))
        anchor_knn_agreement = float(updated.get("anchor_knn_agreement", 1.0))
        recipient_support = float(updated.get("recipient_consensus_score", 0.0))
        neighbor_clade_score = float(updated.get("neighbor_clade_score", 0.0))

        if (
            support_ok
            and species_anchor_score >= NEIGHBOR_CLADE_MIN_SPECIES_ANCHOR
            and neighbor_anchor_purity >= NEIGHBOR_CLADE_MIN_NEIGHBOR_PURITY
            and purity_drop >= NEIGHBOR_CLADE_MIN_PURITY_DROP
            and anchor_knn_agreement <= NEIGHBOR_CLADE_MAX_KNN_AGREEMENT
            and neighbor_clade_score >= NEIGHBOR_CLADE_MIN_SCORE
            and recipient_support >= NEIGHBOR_CLADE_MIN_RECIPIENT_SUPPORT
        ):
            updated["singleton_class"] = "contamination_candidate"
        else:
            updated["singleton_class"] = "ambiguous"
        updated.setdefault("contig_consensus_score", 0.0)
        classified.append(updated)
    return classified


def _classify_recipient_consensus_fallback_proposals(
    proposals: list[dict],
    *,
    contig_marker_context: dict[tuple[str, str], set[str]] | None = None,
) -> list[dict]:
    contig_marker_context = contig_marker_context or {}
    by_marker: dict[str, list[dict]] = {}
    for proposal in proposals:
        by_marker.setdefault(str(proposal.get("marker_name", "")), []).append(dict(proposal))

    selected: dict[tuple[str, str], dict] = {}
    scored: dict[tuple[str, str], dict] = {}
    for marker_name, marker_proposals in by_marker.items():
        candidates = [dict(proposal) for proposal in marker_proposals]
        chosen = _choose_ranked_recipient_candidate(candidates)
        for candidate in candidates:
            key = (marker_name, str(candidate.get("leaf_name", "")))
            if "recipient_rank_score" in candidate:
                scored[key] = {
                    "recipient_rank_score": float(candidate["recipient_rank_score"]),
                    "score": float(candidate["recipient_rank_score"]),
                }
            if "recipient_rank_gap" in candidate:
                scored.setdefault(key, {})["recipient_rank_gap"] = float(candidate["recipient_rank_gap"])
        if chosen is not None:
            finalized = _finalize_singleton_choice(chosen, score_key="recipient_rank_score")
            selected[(marker_name, str(finalized.get("leaf_name", "")))] = finalized

    selected_classified = _classify_singleton_proposals_legacy(
        list(selected.values()),
        contig_marker_context=contig_marker_context,
    )
    selected_lookup: dict[tuple[str, str], dict] = {}
    for updated in selected_classified:
        if str(updated.get("contig_id", "unknown_contig")) == "unknown_contig" and _strong_unknown_contig_signal(updated):
            updated["singleton_class"] = "contamination_candidate"
        updated.setdefault("contig_consensus_score", 0.0)
        selected_lookup[(str(updated.get("marker_name", "")), str(updated.get("leaf_name", "")))] = updated

    classified: list[dict] = []
    for proposal in proposals:
        key = (str(proposal.get("marker_name", "")), str(proposal.get("leaf_name", "")))
        updated = dict(proposal)
        updated.update(scored.get(key, {}))
        if key in selected_lookup:
            updated.update(selected_lookup[key])
        else:
            updated.setdefault("score", float(updated.get("recipient_rank_score", 0.0)))
            updated.setdefault("singleton_class", "clean")
            updated.setdefault("contig_consensus_score", 0.0)
        classified.append(updated)
    return classified


def _classify_gcp_proposals(
    proposals: list[dict],
    *,
    contig_marker_context: dict[tuple[str, str], set[str]] | None = None,
) -> list[dict]:
    """Classify singleton proposals using Genome Consistency Profiling.

    Computes per-genome z-scores across markers, then combines outlier count,
    IsolationForest global anomaly, and HDBSCAN per-genome clustering to score
    each marker. Only the single most deviant marker per genome can be flagged.

    Falls back to recipient_consensus classification if the panel is too small.
    """
    if not proposals:
        return []

    df = pd.DataFrame(proposals)
    genomes = df["genome"].nunique()
    markers = df["marker_name"].nunique()

    if genomes < GCP_MIN_GENOMES or markers < GCP_MIN_MARKERS:
        logger.warning(
            "GCP fallback: panel has %d genomes and %d markers (need >=%d and >=%d). "
            "Using recipient_consensus classification.",
            genomes, markers, GCP_MIN_GENOMES, GCP_MIN_MARKERS,
        )
        return _classify_recipient_consensus_fallback_proposals(
            proposals,
            contig_marker_context=contig_marker_context,
        )

    feature_cols = [c for c in GCP_KEY_FEATURES if c in df.columns]
    genome_z_cols = []
    for col in feature_cols:
        z_col = f"{col}_genome_z"
        genome_means = df.groupby("genome")[col].transform("mean")
        genome_stds = df.groupby("genome")[col].transform("std")
        genome_stds = genome_stds.where(genome_stds > 0, 1.0)
        df[z_col] = (df[col].astype(float) - genome_means) / genome_stds
        df[z_col] = df[z_col].fillna(0.0)
        genome_z_cols.append(z_col)

    z_abs = df[genome_z_cols].abs()
    df["gcp_outlier_count"] = (z_abs > GCP_Z_THRESHOLD).sum(axis=1)
    df["gcp_maxz"] = z_abs.max(axis=1)
    df["gcp_meanz"] = z_abs.mean(axis=1)

    # HDBSCAN per-genome penalty
    df["gcp_hdbscan_penalty"] = 0.0
    try:
        import hdbscan as _hdbscan_mod
        _has_hdbscan = True
    except ImportError:
        _has_hdbscan = False

    if _has_hdbscan:
        feat_cols_present = [c for c in GCP_KEY_FEATURES[:7] if c in df.columns]
        for genome, grp in df.groupby("genome"):
            if len(grp) < 4:
                continue
            X = df.loc[grp.index, feat_cols_present].values.astype(float)
            mu, std = X.mean(axis=0), X.std(axis=0)
            std = np.where(std > 0, std, 1.0)
            X_std = (X - mu) / std
            min_cluster = max(2, len(grp) // 3)
            clusterer = _hdbscan_mod.HDBSCAN(
                min_cluster_size=min_cluster, min_samples=1, metric="euclidean",
            )
            labels = clusterer.fit_predict(X_std)
            if not np.all(labels == -1):
                unique, counts = np.unique(labels[labels >= 0], return_counts=True)
                if len(unique) > 0:
                    main_cluster = unique[np.argmax(counts)]
                    penalty = np.where(labels == main_cluster, 0.0, 1.0)
                    df.loc[grp.index, "gcp_hdbscan_penalty"] = penalty

    # IsolationForest global anomaly
    df["gcp_iforest_score"] = 0.0
    if feature_cols and len(df) >= 10:
        from sklearn.ensemble import IsolationForest
        X_iso = df[feature_cols].values.astype(float)
        mu, std = X_iso.mean(axis=0), X_iso.std(axis=0)
        std = np.where(std > 0, std, 1.0)
        X_iso_std = (X_iso - mu) / std
        iforest = IsolationForest(n_estimators=200, contamination=0.05, random_state=42)
        raw_scores = iforest.fit(X_iso_std).score_samples(X_iso_std)
        score_range = raw_scores.max() - raw_scores.min() + 1e-10
        df["gcp_iforest_score"] = 1 - (raw_scores - raw_scores.min()) / score_range

    # Combine signals
    df["gcp_outlier_count_pct"] = df["gcp_outlier_count"].rank(pct=True)
    df["gcp_maxz_pct"] = df["gcp_maxz"].rank(pct=True)
    df["gcp_meanz_pct"] = df["gcp_meanz"].rank(pct=True)
    df["gcp_iforest_pct"] = df["gcp_iforest_score"].rank(pct=True)

    genome_signal = (
        0.35 * df["gcp_outlier_count_pct"]
        + 0.35 * df["gcp_maxz_pct"]
        + 0.30 * df["gcp_meanz_pct"]
    )
    global_signal = (
        0.5 * df["gcp_iforest_pct"]
        + 0.3 * df["gcp_hdbscan_penalty"]
        + 0.2 * df["gcp_maxz_pct"]
    )
    df["gcp_combined"] = np.sqrt(genome_signal * global_signal)
    df["gcp_score"] = df["gcp_combined"]

    # Classify: only top-per-genome if it passes multi-gate thresholds
    df["singleton_class"] = "clean"
    for genome, grp in df.groupby("genome"):
        if len(grp) <= 2:
            df.loc[grp.index, "singleton_class"] = "ambiguous"
            continue
        top_loc = grp["gcp_combined"].idxmax()
        top_score = df.loc[top_loc, "gcp_combined"]
        top_outlier_count = df.loc[top_loc, "gcp_outlier_count"]
        top_iforest = df.loc[top_loc, "gcp_iforest_score"]
        if (
            top_score >= GCP_COMBINED_THRESHOLD
            and top_outlier_count >= GCP_OUTLIER_COUNT_THRESHOLD
            and top_iforest >= GCP_IFOREST_THRESHOLD
        ):
            df.loc[top_loc, "singleton_class"] = "contamination_candidate"

    classified = []
    for idx, row in df.iterrows():
        updated = dict(proposals[idx])
        updated["singleton_class"] = row["singleton_class"]
        updated["gcp_score"] = float(row["gcp_score"])
        updated["gcp_combined"] = float(row["gcp_combined"])
        updated["score"] = float(row["gcp_score"])
        updated.setdefault("contig_consensus_score", 0.0)
        classified.append(updated)
    return classified


def classify_singleton_proposals(
    proposals: list[dict],
    *,
    contig_marker_context: dict[tuple[str, str], set[str]],
    marker_neighbor_context: dict[tuple[str, str, str], dict[str, float]] | None = None,
    mode: str = "delta_rf",
) -> list[dict]:
    mode = _canonical_singleton_mode(mode)
    if mode == "gcp":
        return _classify_gcp_proposals(
            proposals,
            contig_marker_context=contig_marker_context,
        )
    if mode == "neighbor_clade":
        return _classify_neighbor_clade_proposals(proposals)
    if mode == "neighbor_ml":
        scored, _accepted = select_neighbor_ml_proposals(proposals)
        return scored
    if mode != "contig_consensus":
        classified = _classify_singleton_proposals_legacy(
            proposals,
            contig_marker_context=contig_marker_context,
        )
        if mode in {"composite", "recipient_consensus"}:
            for updated in classified:
                if str(updated.get("contig_id", "unknown_contig")) == "unknown_contig" and _strong_unknown_contig_signal(updated):
                    updated["singleton_class"] = "contamination_candidate"
                updated.setdefault("contig_consensus_score", 0.0)
        return classified

    marker_neighbor_context = marker_neighbor_context or {}
    classified: list[dict] = []
    for proposal in proposals:
        updated = dict(proposal)
        genome = str(updated["genome"])
        contig = str(updated.get("contig_id", "unknown_contig"))
        marker_name = str(updated.get("marker_name", ""))
        key = (genome, contig)
        markers_on_contig = set(contig_marker_context.get(key, set()))
        updated["contig_consensus_score"] = 0.0

        if contig == "unknown_contig":
            updated["singleton_class"] = (
                "contamination_candidate"
                if _strong_unknown_contig_signal(updated)
                else "ambiguous"
            )
            classified.append(updated)
            continue
        if not markers_on_contig or markers_on_contig == {marker_name}:
            updated["singleton_class"] = "contamination_candidate"
            updated["contig_consensus_score"] = 1.0
            classified.append(updated)
            continue

        other_overlaps = [
            float(marker_neighbor_context.get((genome, contig, marker), {}).get("recipient_neighbor_overlap", math.nan))
            for marker in sorted(markers_on_contig - {marker_name})
        ]
        other_overlaps = [value for value in other_overlaps if not math.isnan(value)]
        if len(other_overlaps) < CONTIG_CONSENSUS_MIN_SUPPORT:
            classified.extend(
                _classify_singleton_proposals_legacy(
                    [updated],
                    contig_marker_context=contig_marker_context,
                )
            )
            continue

        median_other_overlap = statistics.median(other_overlaps)
        candidate_overlap = float(updated.get("recipient_neighbor_overlap", updated.get("neighbor_overlap", 1.0)))
        recipient_consensus = float(updated.get("recipient_consensus_score", 0.0))
        updated["contig_consensus_score"] = max(0.0, median_other_overlap - candidate_overlap)

        if median_other_overlap < 0.4:
            updated["singleton_class"] = "ambiguous"
        elif (
            median_other_overlap >= CONTIG_CONSENSUS_HIGH_OVERLAP
            and candidate_overlap <= CONTIG_CONSENSUS_LOW_OVERLAP
            and recipient_consensus >= RECIPIENT_CONSENSUS_Z_THRESHOLD
        ):
            updated["singleton_class"] = "contamination_candidate"
        elif median_other_overlap >= CONTIG_CONSENSUS_HIGH_OVERLAP:
            updated["singleton_class"] = "hgt_candidate"
        else:
            updated["singleton_class"] = "ambiguous"
        classified.append(updated)
    return classified


def select_singleton_proposals(
    proposals: list[dict],
    *,
    genome_marker_counts: dict[str, int],
    min_markers_per_genome: int = 1,
    max_prunes_per_genome: int = 1,
) -> list[dict]:
    budgets = {
        genome: min(
            max(0, int(count) - int(min_markers_per_genome)),
            int(max_prunes_per_genome),
        )
        for genome, count in genome_marker_counts.items()
    }
    accepted: list[dict] = []
    ranked = sorted(
        proposals,
        key=lambda proposal: (
            -float(proposal.get("score", 0.0)),
            -float(proposal.get("delta_rf", 0.0)),
            -float(proposal.get("topoknn_score", 0.0)),
            str(proposal.get("marker_name", "")),
            str(proposal.get("leaf_name", "")),
        ),
    )
    for proposal in ranked:
        genome = str(proposal["genome"])
        if budgets.get(genome, 0) <= 0:
            continue
        budgets[genome] -= 1
        accepted.append(proposal)
    return accepted


def _independent_query_singleton_proposals(
    candidates: list[dict],
    *,
    mode: str,
) -> list[dict]:
    mode = _canonical_singleton_mode(mode)
    if not candidates:
        return []

    scoped_candidates = [dict(candidate) for candidate in candidates]

    if mode == "delta_rf":
        return [
            _finalize_singleton_choice(candidate, score_key="delta_rf")
            for candidate in scoped_candidates
            if float(candidate.get("delta_rf", 0.0)) > 0.0
        ]

    if mode in {"topoknn", "outlier"}:
        return [
            _finalize_singleton_choice(candidate, score_key="topoknn_score")
            for candidate in scoped_candidates
            if float(candidate.get("delta_rf", 0.0)) > 0.0
            and float(candidate.get("topoknn_score", 0.0)) >= 0.75
        ]

    if mode == "hybrid":
        positive_rf = [
            candidate
            for candidate in scoped_candidates
            if float(candidate.get("delta_rf", 0.0)) > 0.0
        ]
        if not positive_rf:
            return []
        rf_norm = _normalize_candidate_metric(positive_rf, "delta_rf")
        topoknn_norm = _normalize_candidate_metric(positive_rf, "topoknn_score")
        branch_norm = _normalize_candidate_metric(positive_rf, "branch_outlier")
        bitscore_norm = _normalize_candidate_metric(positive_rf, "bitscore_outlier")
        for candidate in positive_rf:
            leaf_name = candidate["leaf_name"]
            candidate["hybrid_score"] = (
                1.5 * rf_norm.get(leaf_name, 0.0)
                + 1.5 * topoknn_norm.get(leaf_name, 0.0)
                + 0.25 * branch_norm.get(leaf_name, 0.0)
                + 0.25 * bitscore_norm.get(leaf_name, 0.0)
            )
        return [
            _finalize_singleton_choice(candidate, score_key="hybrid_score")
            for candidate in positive_rf
            if float(candidate.get("delta_rf", 0.0)) >= 0.05
            and float(candidate.get("topoknn_score", 0.0)) >= 1.0
        ]

    if mode in {"composite", "contig_consensus"}:
        _choose_composite_candidate(scoped_candidates)
        return [
            _finalize_singleton_choice(candidate, score_key="composite_score")
            for candidate in scoped_candidates
            if _composite_candidate_is_eligible(candidate)
        ]

    if mode == "recipient_consensus":
        return [
            _finalize_singleton_choice(candidate, score_key="recipient_rank_score")
            for candidate in _recipient_ranked_candidates(scoped_candidates)
        ]

    if mode == "gcp":
        # GCP needs all candidates; classification happens globally in _classify_gcp_proposals
        return [
            {**candidate, "score": 0.0}
            for candidate in scoped_candidates
        ]

    return []


def singleton_proposals_from_results(
    results: list[dict],
    *,
    mode: str,
    reference_genomes: set[str] | None = None,
) -> tuple[list[dict], set[tuple[str, str]]]:
    mode = _canonical_singleton_mode(mode)
    reference_genomes = reference_genomes or set()
    proposals: list[dict] = []
    proposal_keys: set[tuple[str, str]] = set()
    for result in results:
        marker_name = str(result.get("marker_name", ""))
        if mode in {"neighbor_clade", "neighbor_ml", "gcp"}:
            score_key = "neighbor_clade_score" if mode != "gcp" else "delta_rf"
            for candidate in result.get("candidates", []):
                leaf_name = str(candidate.get("leaf_name", ""))
                if not leaf_name:
                    continue
                genome = str(candidate.get("genome", ""))
                if reference_genomes and genome in reference_genomes:
                    continue
                proposals.append(
                    {
                        **candidate,
                        "marker_name": marker_name,
                        "score": float(candidate.get(score_key, 0.0)),
                    }
                )
                proposal_keys.add((marker_name, leaf_name))
            continue

        if reference_genomes:
            query_candidates = [
                candidate
                for candidate in result.get("candidates", [])
                if str(candidate.get("genome", "")) not in reference_genomes
            ]
            for candidate in _independent_query_singleton_proposals(
                query_candidates,
                mode=mode,
            ):
                leaf_name = str(candidate.get("leaf_name", ""))
                if not leaf_name:
                    continue
                proposals.append(
                    {
                        **candidate,
                        "marker_name": marker_name,
                    }
                )
                proposal_keys.add((marker_name, leaf_name))
            continue

        chosen = result.get("chosen")
        if chosen is None:
            continue
        leaf_name = str(chosen.get("leaf_name", ""))
        if not leaf_name:
            continue
        proposals.append(
            {
                **chosen,
                "marker_name": marker_name,
            }
        )
        proposal_keys.add((marker_name, leaf_name))
    return proposals, proposal_keys


def _reference_genomes_from_dir(ref_dir: str | None) -> set[str]:
    if not ref_dir or not os.path.isdir(ref_dir):
        return set()
    return {
        os.path.basename(path).split(".")[0]
        for path in glob.glob(os.path.join(ref_dir, "*"))
        if os.path.isfile(path)
    }


def _prune_tree_to_query_genomes(tree: Tree, reference_genomes: set[str]) -> Tree:
    pruned = tree.copy()
    keep_leaves = [
        leaf.name
        for leaf in pruned.iter_leaves()
        if leaf.name.split("|")[0] not in reference_genomes
    ]
    if keep_leaves and len(keep_leaves) < len(list(pruned.iter_leaves())):
        pruned.prune(keep_leaves)
    return pruned


def _filter_reference_singleton_proposals(
    proposals: list[dict],
    proposal_keys: set[tuple[str, str]],
    *,
    reference_genomes: set[str],
) -> tuple[list[dict], set[tuple[str, str]]]:
    if not reference_genomes:
        return proposals, proposal_keys
    filtered = [
        proposal
        for proposal in proposals
        if str(proposal.get("genome", "")) not in reference_genomes
    ]
    filtered_pairs = {
        (str(proposal.get("marker_name", "")), str(proposal.get("leaf_name", "")))
        for proposal in filtered
    }
    filtered_keys = {
        key
        for key in proposal_keys
        if key in filtered_pairs
    }
    return filtered, filtered_keys


def _rank01_series(series: pd.Series) -> pd.Series:
    if series.empty:
        return series
    return series.rank(method="average", pct=True).fillna(0.0)


def _group_robust_zscores(
    df: pd.DataFrame,
    *,
    group_cols: list[str],
    feature: str,
) -> pd.Series:
    def _transform(series: pd.Series) -> pd.Series:
        median = float(series.median())
        mad = float((series - median).abs().median())
        if mad > 0:
            return 0.6745 * (series - median) / mad
        stdev = float(series.std(ddof=0))
        if stdev > 0:
            return (series - median) / stdev
        return pd.Series(np.zeros(len(series)), index=series.index)

    return (
        df.groupby(group_cols, observed=True)[feature]
        .transform(_transform)
        .fillna(0.0)
    )


_ML_NUMERIC_COLS = [
    "delta_rf",
    "neighbor_overlap",
    "topoknn_score",
    "branch_outlier",
    "bitscore_outlier",
    "recipient_consensus_score",
    "species_anchor_score",
    "species_anchor_purity",
    "species_anchor_compactness_score",
    "species_long_branch_z",
    "species_long_branch_support",
    "present_neighbor_count",
    "present_neighbor_fraction",
    "neighbor_anchor_purity",
    "join_purity",
    "purity_drop",
    "anchor_knn_agreement",
    "attachment_gap",
]

_ML_EMPTY_TRAIN_SCORE_COLUMNS = [
    "iforest_anomaly",
    "mahalanobis",
    "genome_outlier_score",
    "marker_outlier_score",
    "ensemble_score_v1",
    "genome_shape_signal",
    "shape_support",
    "anchor_support",
    "genome_focus_score_v2",
    "marker_top1_frequency",
    "marker_rarity_weight",
    "marker_penalized_score_v3",
    "shape_penalized_score_v4",
    "directed_support_v7",
    "genome_gap_v8",
    "genome_first_score_v8",
]


def _prepare_ml_features(
    proposals: list[dict],
) -> tuple[pd.DataFrame, list[str], list[str]]:
    df = pd.DataFrame(proposals).copy()
    for col in _ML_NUMERIC_COLS:
        df[col] = pd.to_numeric(df.get(col, 0.0), errors="coerce").fillna(0.0)

    df["panel_id"] = "runtime_panel"
    df["knn_disagreement"] = 1.0 - df["anchor_knn_agreement"]
    df["join_impurity"] = 1.0 - df["join_purity"]
    df["neighbor_impurity"] = 1.0 - df["neighbor_anchor_purity"]
    df["present_neighbor_deficit"] = 1.0 - df["present_neighbor_fraction"]
    df["high_confidence"] = (
        (df["present_neighbor_count"] >= NEIGHBOR_ML_HIGH_CONF_MIN_PRESENT)
        & (df["present_neighbor_fraction"] >= NEIGHBOR_ML_HIGH_CONF_MIN_FRACTION)
        & (df["species_anchor_score"] >= NEIGHBOR_ML_HIGH_CONF_MIN_ANCHOR)
    )

    raw_features = _ML_NUMERIC_COLS + [
        "knn_disagreement",
        "join_impurity",
        "neighbor_impurity",
        "present_neighbor_deficit",
    ]
    for feature in raw_features:
        df[f"zg_{feature}"] = _group_robust_zscores(df, group_cols=["genome"], feature=feature)
        df[f"zm_{feature}"] = _group_robust_zscores(df, group_cols=["marker_name"], feature=feature)

    feature_cols = list(raw_features)
    feature_cols.extend(f"zg_{feature}" for feature in raw_features)
    feature_cols.extend(f"zm_{feature}" for feature in raw_features)

    return df, feature_cols, raw_features


def _fit_anomaly_models(df: pd.DataFrame, feature_cols: list[str]) -> bool:
    """Fit IsolationForest + MinCovDet on high-confidence rows.

    Mutates df with ``iforest_anomaly`` and ``mahalanobis``. Returns True
    if models were fit, False if the high-confidence slice was empty; in
    the empty case, also writes zeros for every downstream score column so
    the orchestrator can return early without calling the ranking stage.
    """
    from sklearn.covariance import MinCovDet
    from sklearn.ensemble import IsolationForest
    from sklearn.preprocessing import RobustScaler

    train = (
        df.loc[df["high_confidence"], feature_cols]
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0.0)
    )
    if train.empty:
        for col in _ML_EMPTY_TRAIN_SCORE_COLUMNS:
            df[col] = 0.0
        return False

    scaler = RobustScaler()
    x_train = scaler.fit_transform(train)
    iso = IsolationForest(
        n_estimators=400,
        contamination=0.01,
        random_state=42,
        n_jobs=-1,
    )
    iso.fit(x_train)
    df["iforest_anomaly"] = 0.0
    df.loc[df["high_confidence"], "iforest_anomaly"] = -iso.score_samples(x_train)
    try:
        mcd = MinCovDet(random_state=42, support_fraction=0.8)
        mcd.fit(x_train)
        df["mahalanobis"] = 0.0
        df.loc[df["high_confidence"], "mahalanobis"] = mcd.mahalanobis(x_train)
    except Exception as exc:
        logger.warning(
            "MinCovDet fit failed; falling back to mahalanobis=0.0 (%s)",
            exc,
            exc_info=True,
        )
        df["mahalanobis"] = 0.0
    return True


def _rank_genome_proposals(
    df: pd.DataFrame, raw_features: list[str]
) -> pd.DataFrame:
    """Derive ensemble score, shape signal, marker-rarity penalty, and
    per-genome top-row rankings from the fitted anomaly columns.

    Returns a (possibly new) DataFrame because the genome_first_score_v8
    back-merge replaces ``df``.
    """
    genome_score_cols = [f"zg_{feature}" for feature in raw_features]
    marker_score_cols = [f"zm_{feature}" for feature in raw_features]
    df["genome_outlier_score"] = df[genome_score_cols].abs().mean(axis=1)
    df["marker_outlier_score"] = df[marker_score_cols].abs().mean(axis=1)
    score_parts = [
        _rank01_series(df["iforest_anomaly"].fillna(0.0)),
        _rank01_series(df["mahalanobis"].fillna(0.0)),
        _rank01_series(df["genome_outlier_score"].fillna(0.0)),
        _rank01_series(df["marker_outlier_score"].fillna(0.0)),
    ]
    df["ensemble_score_v1"] = sum(score_parts) / len(score_parts)
    df.loc[~df["high_confidence"], "ensemble_score_v1"] = 0.0

    directed_features = [
        "purity_drop",
        "recipient_consensus_score",
        "attachment_gap",
        "delta_rf",
        "knn_disagreement",
        "join_impurity",
    ]
    positive_genome_z = [df[f"zg_{feature}"].clip(lower=0.0) for feature in directed_features]
    df["genome_shape_signal"] = sum(positive_genome_z) / len(positive_genome_z)
    df["shape_support"] = (
        _rank01_series(df["purity_drop"].fillna(0.0))
        + _rank01_series(df["knn_disagreement"].fillna(0.0))
        + _rank01_series(df["recipient_consensus_score"].fillna(0.0))
        + _rank01_series(df["join_impurity"].fillna(0.0))
    ) / 4.0
    df["anchor_support"] = df["species_anchor_score"].clip(lower=0.0) * df["present_neighbor_fraction"].clip(lower=0.0)
    df["genome_focus_score_v2"] = (
        0.6 * _rank01_series(df["genome_shape_signal"].fillna(0.0))
        + 0.25 * _rank01_series(df["iforest_anomaly"].fillna(0.0))
        + 0.15 * _rank01_series(df["mahalanobis"].fillna(0.0))
    ) * (0.5 + 0.5 * df["anchor_support"])
    df.loc[~df["high_confidence"], "genome_focus_score_v2"] = 0.0

    top_marker_rows = (
        df.loc[df["high_confidence"]]
        .sort_values("genome_shape_signal", ascending=False)
        .groupby(["genome"], observed=True, as_index=False)
        .head(1)
    )
    top_counts = (
        top_marker_rows.groupby(["marker_name"], observed=True)
        .size()
        .rename("marker_top1_count")
        .reset_index()
    )
    panel_genome_count = max(1, int(df.loc[df["high_confidence"], "genome"].nunique()))
    top_counts["marker_top1_frequency"] = top_counts["marker_top1_count"] / panel_genome_count
    df = df.merge(top_counts[["marker_name", "marker_top1_frequency"]], on=["marker_name"], how="left")
    df["marker_top1_frequency"] = df["marker_top1_frequency"].fillna(0.0)
    df["marker_rarity_weight"] = 1.0 - df["marker_top1_frequency"].clip(upper=0.95)
    df["marker_penalized_score_v3"] = df["genome_focus_score_v2"] * df["marker_rarity_weight"]
    df.loc[~df["high_confidence"], "marker_penalized_score_v3"] = 0.0
    df["shape_penalized_score_v4"] = df["marker_penalized_score_v3"] * (0.5 + 0.5 * df["shape_support"])
    df.loc[~df["high_confidence"], "shape_penalized_score_v4"] = 0.0

    df["directed_support_v7"] = (
        _rank01_series(df["purity_drop"].fillna(0.0))
        + _rank01_series(df["recipient_consensus_score"].fillna(0.0))
        + _rank01_series(df["knn_disagreement"].fillna(0.0))
        + _rank01_series(df["join_impurity"].fillna(0.0))
    ) / 4.0
    top_rows = (
        df.loc[df["high_confidence"]]
        .sort_values(["genome", "shape_penalized_score_v4"], ascending=[True, False])
        .groupby(["genome"], observed=True, as_index=False)
        .head(1)
        .copy()
    )
    second_rows = (
        df.loc[df["high_confidence"]]
        .sort_values(["genome", "shape_penalized_score_v4"], ascending=[True, False])
        .groupby(["genome"], observed=True)
        .nth(1)
        .reset_index()[["genome", "shape_penalized_score_v4"]]
        .rename(columns={"shape_penalized_score_v4": "genome_second_score_v8"})
    )
    top_rows = top_rows.merge(second_rows, on=["genome"], how="left")
    top_rows["genome_second_score_v8"] = top_rows["genome_second_score_v8"].fillna(0.0)
    top_rows["genome_gap_v8"] = top_rows["shape_penalized_score_v4"] - top_rows["genome_second_score_v8"]
    top_rows["genome_top_rank_v8"] = _rank01_series(top_rows["shape_penalized_score_v4"])
    top_rows["genome_gap_rank_v8"] = _rank01_series(top_rows["genome_gap_v8"])
    top_rows["genome_directed_rank_v8"] = _rank01_series(top_rows["directed_support_v7"])
    top_rows["genome_first_score_v8"] = (
        0.25 * top_rows["genome_top_rank_v8"]
        + 0.75 * top_rows["genome_gap_rank_v8"]
        + 0.25 * top_rows["genome_directed_rank_v8"]
    )
    df["genome_gap_v8"] = 0.0
    df["genome_first_score_v8"] = 0.0
    df = df.merge(
        top_rows[["genome", "leaf_name", "genome_gap_v8", "genome_first_score_v8"]],
        on=["genome", "leaf_name"],
        how="left",
        suffixes=("", "_top"),
    )
    df["genome_gap_v8"] = df["genome_gap_v8_top"].fillna(df["genome_gap_v8"])
    df["genome_first_score_v8"] = df["genome_first_score_v8_top"].fillna(df["genome_first_score_v8"])
    df = df.drop(columns=["genome_gap_v8_top", "genome_first_score_v8_top"])
    df.loc[~df["high_confidence"], "genome_first_score_v8"] = 0.0
    return df


def score_neighbor_ml_proposals(proposals: list[dict]) -> list[dict]:
    if not proposals:
        return []

    df, feature_cols, raw_features = _prepare_ml_features(proposals)
    if not _fit_anomaly_models(df, feature_cols):
        return df.to_dict("records")
    df = _rank_genome_proposals(df, raw_features)
    return df.to_dict("records")


def select_neighbor_ml_proposals(proposals: list[dict]) -> tuple[list[dict], list[dict]]:
    scored = score_neighbor_ml_proposals(proposals)
    if not scored:
        return [], []

    df = pd.DataFrame(scored)
    panel_taxa = int(df["taxa_count"].max()) if "taxa_count" in df.columns and not df.empty else 0
    policy_variant = "genome_first_score_v8" if panel_taxa <= 60 else "shape_penalized_score_v4"
    top_genomes = (
        df.loc[df["high_confidence"]]
        .sort_values(policy_variant, ascending=False)
        .groupby("genome", observed=True, as_index=False)
        .head(1)
        .sort_values(policy_variant, ascending=False)
    )
    selected_keys = {
        (str(row.marker_name), str(row.leaf_name))
        for row in top_genomes.itertuples(index=False)
    }

    updated: list[dict] = []
    accepted: list[dict] = []
    for record in scored:
        item = dict(record)
        item["policy_variant"] = policy_variant
        item["policy_score"] = float(item.get(policy_variant, 0.0))
        item["score"] = item["policy_score"]
        key = (str(item["marker_name"]), str(item["leaf_name"]))
        if key in selected_keys and bool(item.get("high_confidence", False)):
            item["singleton_class"] = "contamination_candidate"
            accepted.append(item)
        else:
            item["singleton_class"] = "ambiguous"
        updated.append(item)
    return updated, accepted


def build_singleton_output_tree(
    *,
    marker_tree_path: str,
    species_tree_path: str,
    accepted_leaf_names: list[str],
    mode: str | None = None,
) -> tuple[Tree, str]:
    if mode is not None and not singleton_mode_uses_global_rf_gate(mode):
        raise ValueError(f"Singleton mode {mode!r} bypasses the RF safeguard")
    tf = Tree(marker_tree_path)
    accepted_leaf_set = {str(name) for name in accepted_leaf_names if str(name)}
    if not accepted_leaf_set:
        return tf.copy(), "kept"

    remaining = [
        leaf.name
        for leaf in tf.iter_leaves()
        if leaf.name not in accepted_leaf_set
    ]
    if not remaining:
        return tf.copy(), "kept"

    candidate_tree = tf.copy()
    candidate_tree.prune(remaining)
    td = _tree_to_genome_level(tf)
    ti = _cached_species_tree(species_tree_path)
    ti.prune([leaf.name for leaf in td.iter_leaves()])
    chosen_tree = choose_tree_by_rf(ti, tf, candidate_tree)
    if sorted(leaf.name for leaf in chosen_tree.iter_leaves()) == sorted(leaf.name for leaf in tf.iter_leaves()):
        return chosen_tree, "kept_rf_guard"
    return chosen_tree, "pruned"


def effective_singleton_mode(
    mode: str,
    rdist: float,
    *,
    duplicate_resolution_present: bool,
) -> str:
    """Use the chosen singleton cleanup implementation for all runtime runs."""
    return _canonical_singleton_mode(mode)


def _propose_singleton_prune_worker(args):
    (
        filepath,
        species_tree_path,
        num_nei_override,
        singles_min_rfdist,
        singles_mode,
        table_path,
        duplicate_markers,
        alignment_path,
        ref_dir,
    ) = args
    tf = Tree(filepath)
    td = _tree_to_genome_level(tf)
    td_leaves = list(td.iter_leaves())
    ti = _cached_species_tree(species_tree_path)
    ti.prune([n.name for n in td_leaves])

    rf, maxrf, *_ = ti.robinson_foulds(td, unrooted_trees=True)
    maxrf = maxrf + 0.0001
    rdist = rf / maxrf
    effective_mode = effective_singleton_mode(
        singles_mode,
        rdist,
        duplicate_resolution_present=bool(duplicate_markers),
    )
    leaf_count = len(list(tf.iter_leaves()))
    if num_nei_override > 0:
        num_nei = min(num_nei_override, max(1, leaf_count - 1))
    else:
        num_nei = max(2, min(10, leaf_count - 1))

    if effective_mode == "loo_profile":
        candidates = []
        for leaf in sorted(tf.iter_leaves(), key=lambda item: str(item.name)):
            genome_id, contig_id, gene_id = parse_sequence_id(str(leaf.name))
            candidates.append(
                {
                    "leaf_name": str(leaf.name),
                    "taxa_count": leaf_count,
                    "genome": genome_id,
                    "contig_id": contig_id,
                    "gene_id": gene_id,
                }
            )
        return {
            "filepath": filepath,
            "marker_name": _marker_name_from_tree_path(filepath),
            "mode": effective_mode,
            "rdist": rdist,
            "num_nei": num_nei,
            "chosen": None,
            "candidates": candidates,
        }

    score_table = None
    score_col = None
    if os.path.exists(table_path):
        score_table, score_col = _cached_score_table(table_path)
    candidates = _score_singleton_candidates(
        species_tree=ti,
        working_tree=tf,
        k=num_nei,
        score_table=score_table,
        score_col=score_col,
        alignment_path=alignment_path,
    )
    chosen = None
    if effective_mode == "neighbor_ml":
        chosen = None
    elif rdist >= singles_min_rfdist or effective_mode in {"neighbor_clade", "gcp"}:
        chosen = choose_singleton_prune(
            species_tree=ti,
            working_tree=tf,
            mode=effective_mode,
            k=num_nei,
            score_table=score_table,
            score_col=score_col,
            alignment_path=alignment_path,
            scored_candidates=candidates,
        )
        if chosen is not None:
            chosen = {
                key: value
                for key, value in chosen.items()
                if key != "candidate_tree"
            }
            genome_id, contig_id, gene_id = parse_sequence_id(chosen["leaf_name"])
            chosen["genome"] = genome_id
            chosen["contig_id"] = contig_id
            chosen["gene_id"] = gene_id

    return {
        "filepath": filepath,
        "marker_name": _marker_name_from_tree_path(filepath),
        "mode": effective_mode,
        "rdist": rdist,
        "num_nei": num_nei,
        "chosen": chosen,
        "candidates": [
            {
                key: value
                for key, value in candidate.items()
                if key != "candidate_tree"
            }
            for candidate in candidates
        ],
    }


def _singleton_decision(
    chosen: dict | None,
    *,
    marker_name: str,
    accepted_keys: set[tuple[str, str]],
) -> str:
    if chosen is None:
        return "kept"
    leaf_name = str(chosen.get("leaf_name", ""))
    if (marker_name, leaf_name) in accepted_keys:
        return "pruned"
    singleton_class = chosen.get("singleton_class", "unclassified")
    if singleton_class == "hgt_candidate":
        return "kept_hgt_candidate"
    if singleton_class == "ambiguous":
        return "kept_ambiguous"
    return "blocked_by_genome_budget"


def _write_singleton_result(
    result: dict,
    *,
    species_tree_path: str,
    outdir: str,
    accepted_keys: set[tuple[str, str]],
) -> None:
    filepath = result["filepath"]
    marker_name = result["marker_name"]
    classified_candidates = result.get("classified_candidates")
    if classified_candidates is None:
        chosen = result.get("chosen")
        classified_candidates = [chosen] if chosen is not None else []
    accepted_leaf_names = [
        str(candidate.get("leaf_name", ""))
        for candidate in classified_candidates
        if (marker_name, str(candidate.get("leaf_name", ""))) in accepted_keys
    ]
    chosen_tree, overall_decision = build_singleton_output_tree(
        marker_tree_path=filepath,
        species_tree_path=species_tree_path,
        accepted_leaf_names=accepted_leaf_names,
        mode=str(result.get("mode", "")),
    )

    with open(os.path.join(outdir, "removed", marker_name), "a") as f:
        if not classified_candidates:
            f.write(
                f"\nno_singleton_prune\tmode={result['mode']}\trdist={result['rdist']:.3f}\tnum_nei={result['num_nei']}\n"
            )
        else:
            for candidate in classified_candidates:
                candidate_leaf = str(candidate.get("leaf_name", ""))
                decision = (
                    overall_decision
                    if (marker_name, candidate_leaf) in accepted_keys and overall_decision == "kept_rf_guard"
                    else _singleton_decision(
                        candidate,
                        marker_name=marker_name,
                        accepted_keys=accepted_keys,
                    )
                )
                f.write(
                    (
                        f"\nsin{candidate_leaf}\tmode={result['mode']}"
                        f"\tdecision={decision}"
                        f"\tsingleton_class={candidate.get('singleton_class', 'unclassified')}"
                        f"\tscore={float(candidate.get('score', candidate.get('neighbor_clade_score', 0.0))):.3f}"
                        f"\tdelta_rf={float(candidate['delta_rf']):.3f}"
                        f"\ttopoknn={float(candidate['topoknn_score']):.3f}"
                        f"\toverlap={float(candidate['neighbor_overlap']):.3f}"
                        f"\tanchor_score={float(candidate.get('species_anchor_score', 0.0)):.3f}"
                        f"\tpurity_drop={float(candidate.get('purity_drop', 0.0)):.3f}"
                        f"\tneighbor_clade={float(candidate.get('neighbor_clade_score', 0.0)):.3f}"
                        f"\tbranch_outlier={float(candidate['branch_outlier']):.3f}"
                        f"\tbitscore_outlier={float(candidate['bitscore_outlier']):.3f}\n"
                    )
                )

    chosen_tree.write(
        format=1,
        outfile=os.path.join(
            outdir,
            "protTrees",
            "no_singles",
            os.path.basename(filepath).split(".")[0] + ".nw",
        ),
    )


def _write_singleton_candidate_table(
    results: list[dict],
    *,
    outdir: str,
    accepted_keys: set[tuple[str, str]],
    proposal_lookup: dict[tuple[str, str], dict],
    proposal_keys: set[tuple[str, str]],
) -> None:
    path = os.path.join(outdir, "singleton_candidates.tsv")

    def format_float(value):
        return "" if value is None else f"{float(value):.6f}"

    fieldnames = [
        "marker_name",
        "leaf_name",
        "genome",
        "contig_id",
        "gene_id",
        "mode",
        "rdist",
        "num_nei",
        "delta_rf",
        "neighbor_overlap",
        "topoknn_score",
        "branch_outlier",
        "bitscore_outlier",
        "recipient_consensus_score",
        "contig_consensus_score",
        "species_anchor_score",
        "species_anchor_purity",
        "species_anchor_compactness",
        "species_anchor_compactness_score",
        "species_long_branch_z",
        "species_long_branch_support",
        "target_neighbor_count",
        "present_neighbor_count",
        "present_neighbor_fraction",
        "neighbor_anchor_purity",
        "join_purity",
        "purity_drop",
        "anchor_knn_agreement",
        "attachment_gap",
        "neighbor_clade_score",
        "loo_class",
        "loo_abstention_reason",
        "loo_target_support",
        "loo_target_support_sides",
        "loo_attachment_taxa",
        "loo_attachment_clade",
        "loo_voter_count",
        "loo_voter_markers",
        "loo_coordinate_count",
        "loo_coordinate_taxa",
        "loo_target_discordance",
        "loo_voter_center",
        "loo_voter_mad",
        "loo_voter_upper",
        "loo_conflict_margin",
        "loo_score",
        "loo_marker_rank",
        "loo_marker_margin",
        "loo_conflict_beyond_dispersion",
        "loo_decision",
        "selected_candidate",
        "decision",
        "singleton_class",
        "truth_label",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for result in sorted(results, key=lambda item: item["filepath"]):
            chosen = result.get("chosen")
            for candidate in result.get("candidates", []):
                marker_name = str(result["marker_name"])
                leaf_name = str(candidate.get("leaf_name", ""))
                selected_candidate = (marker_name, leaf_name) in proposal_keys
                enriched = proposal_lookup.get((marker_name, leaf_name), {})
                is_loo_profile = str(result.get("mode", "")) == "loo_profile"
                numeric_default = None if is_loo_profile else 0.0
                count_default = None if is_loo_profile else 0
                conflict_beyond_dispersion = enriched.get("loo_conflict_beyond_dispersion")
                writer.writerow(
                    {
                        "marker_name": marker_name,
                        "leaf_name": leaf_name,
                        "genome": candidate.get("genome", ""),
                        "contig_id": candidate.get("contig_id", ""),
                        "gene_id": candidate.get("gene_id", ""),
                        "mode": result.get("mode", ""),
                        "rdist": f"{float(result.get('rdist', 0.0)):.6f}",
                        "num_nei": int(result.get("num_nei", 0)),
                        "delta_rf": format_float(candidate.get("delta_rf", numeric_default)),
                        "neighbor_overlap": format_float(candidate.get("neighbor_overlap", numeric_default)),
                        "topoknn_score": format_float(candidate.get("topoknn_score", numeric_default)),
                        "branch_outlier": format_float(candidate.get("branch_outlier", numeric_default)),
                        "bitscore_outlier": format_float(candidate.get("bitscore_outlier", numeric_default)),
                        "recipient_consensus_score": format_float(candidate.get("recipient_consensus_score", numeric_default)),
                        "contig_consensus_score": format_float(enriched.get("contig_consensus_score", numeric_default)),
                        "species_anchor_score": format_float(candidate.get("species_anchor_score", numeric_default)),
                        "species_anchor_purity": format_float(candidate.get("species_anchor_purity", numeric_default)),
                        "species_anchor_compactness": format_float(candidate.get("species_anchor_compactness", numeric_default)),
                        "species_anchor_compactness_score": format_float(candidate.get("species_anchor_compactness_score", numeric_default)),
                        "species_long_branch_z": format_float(candidate.get("species_long_branch_z", numeric_default)),
                        "species_long_branch_support": format_float(candidate.get("species_long_branch_support", numeric_default)),
                        "target_neighbor_count": candidate.get("target_neighbor_count", count_default),
                        "present_neighbor_count": candidate.get("present_neighbor_count", count_default),
                        "present_neighbor_fraction": format_float(candidate.get("present_neighbor_fraction", numeric_default)),
                        "neighbor_anchor_purity": format_float(candidate.get("neighbor_anchor_purity", numeric_default)),
                        "join_purity": format_float(candidate.get("join_purity", numeric_default)),
                        "purity_drop": format_float(candidate.get("purity_drop", numeric_default)),
                        "anchor_knn_agreement": format_float(candidate.get("anchor_knn_agreement", numeric_default)),
                        "attachment_gap": format_float(candidate.get("attachment_gap", numeric_default)),
                        "neighbor_clade_score": format_float(candidate.get("neighbor_clade_score", numeric_default)),
                        "loo_class": enriched.get("loo_class", ""),
                        "loo_abstention_reason": enriched.get("loo_abstention_reason") or "",
                        "loo_target_support": format_float(enriched.get("loo_target_support")),
                        "loo_target_support_sides": ";".join(
                            f"{','.join(str(taxon) for taxon in side.get('taxa', []))}:"
                            f"{format_float(side.get('support'))}"
                            for side in enriched.get("loo_target_support_sides", [])
                        ),
                        "loo_attachment_taxa": ",".join(
                            sorted(str(item) for item in enriched.get("loo_attachment_taxa", []))
                        ),
                        "loo_attachment_clade": enriched.get("loo_attachment_clade") or "",
                        "loo_voter_count": enriched.get("loo_voter_count", ""),
                        "loo_voter_markers": ",".join(
                            sorted(str(item) for item in enriched.get("loo_voter_markers", []))
                        ),
                        "loo_coordinate_count": enriched.get("loo_coordinate_count", ""),
                        "loo_coordinate_taxa": ",".join(
                            sorted(str(item) for item in enriched.get("loo_coordinate_taxa", []))
                        ),
                        "loo_target_discordance": format_float(enriched.get("loo_target_discordance")),
                        "loo_voter_center": format_float(enriched.get("loo_voter_center")),
                        "loo_voter_mad": format_float(enriched.get("loo_voter_mad")),
                        "loo_voter_upper": format_float(enriched.get("loo_voter_upper")),
                        "loo_conflict_margin": format_float(enriched.get("loo_conflict_margin")),
                        "loo_score": format_float(enriched.get("loo_score")),
                        "loo_marker_rank": enriched.get("loo_marker_rank") or "",
                        "loo_marker_margin": format_float(enriched.get("loo_marker_margin")),
                        "loo_conflict_beyond_dispersion": (
                            ""
                            if conflict_beyond_dispersion is None
                            else "yes" if conflict_beyond_dispersion else "no"
                        ),
                        "loo_decision": enriched.get("loo_decision", ""),
                        "selected_candidate": "yes" if selected_candidate else "no",
                        "decision": (
                            enriched.get("loo_decision", "")
                            if is_loo_profile
                            else _singleton_decision(
                                enriched,
                                marker_name=marker_name,
                                accepted_keys=accepted_keys,
                            )
                            if selected_candidate
                            else "not_selected"
                        ),
                        "singleton_class": enriched.get("singleton_class", ""),
                        "truth_label": "",
                    }
                )


def remove_singles(cfg: Config, species_tree_path: str | None = None):
    """Remove singleton contaminants with the delta-RF proposal workflow."""
    files = glob.glob(os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out", "*"))
    if species_tree_path is None:
        species_tree_path = os.path.join(cfg.outdir, "tree.nwk")
    duplicate_markers = {
        marker
        for marker, _genome in _load_kept_assignments(
            os.path.join(cfg.outdir, "marker_selection_rf_values.txt")
        )
    }
    args = [
        (
            f,
            species_tree_path,
            cfg.num_nei,
            cfg.singles_min_rfdist,
            cfg.singles_mode,
            os.path.join(cfg.outdir, "table_elim_dups"),
            duplicate_markers,
            os.path.join(cfg.outdir, "trimmed_protTrees", f"{_marker_name_from_tree_path(f)}.faa"),
            cfg.ref,
        )
        for f in files
    ]
    results = map_processed(_propose_singleton_prune_worker, args, cfg.num_cpus)
    reference_genomes = _reference_genomes_from_dir(cfg.ref)
    singleton_mode = _canonical_singleton_mode(cfg.singles_mode)
    if singleton_mode == "loo_profile":
        marker_trees = {
            str(result["marker_name"]): Tree(str(result["filepath"]), format=1)
            for result in sorted(results, key=lambda item: str(item["filepath"]))
        }
        loo_rows = score_loo_profiles(
            marker_trees,
            excluded_target_genomes=reference_genomes,
        )
        proposal_lookup = {
            (str(row["marker_name"]), str(row["leaf_name"])): row
            for row in loo_rows
        }
        accepted_keys: set[tuple[str, str]] = set()
        proposal_keys: set[tuple[str, str]] = set()
        _write_singleton_candidate_table(
            results,
            outdir=cfg.outdir,
            accepted_keys=accepted_keys,
            proposal_lookup=proposal_lookup,
            proposal_keys=proposal_keys,
        )
        no_singles_dir = os.path.join(cfg.outdir, "protTrees", "no_singles")
        os.makedirs(no_singles_dir, exist_ok=True)
        for result in sorted(results, key=lambda item: str(item["filepath"])):
            source = str(result["filepath"])
            shutil.copyfile(source, os.path.join(no_singles_dir, os.path.basename(source)))
        return

    proposals, proposal_keys = singleton_proposals_from_results(
        results,
        mode=cfg.singles_mode,
        reference_genomes=reference_genomes,
    )
    proposals, proposal_keys = _filter_reference_singleton_proposals(
        proposals,
        proposal_keys,
        reference_genomes=reference_genomes,
    )
    if singleton_mode in {"neighbor_ml", "gcp"}:
        proposals = classify_singleton_proposals(
            proposals,
            contig_marker_context=(
                _load_contig_marker_context(os.path.join(cfg.outdir, "table_elim_dups"))
                if singleton_mode == "gcp"
                else {}
            ),
            marker_neighbor_context=None,
            mode=cfg.singles_mode,
        )
        accepted_proposals = [
            proposal
            for proposal in proposals
            if proposal.get("singleton_class") == "contamination_candidate"
        ]
        if singleton_mode == "gcp":
            accepted = select_singleton_proposals(
                accepted_proposals,
                genome_marker_counts=_count_genome_marker_support(files),
                min_markers_per_genome=1,
            )
        else:
            accepted = accepted_proposals
    else:
        proposals = classify_singleton_proposals(
            proposals,
            contig_marker_context=_load_contig_marker_context(os.path.join(cfg.outdir, "table_elim_dups")),
            marker_neighbor_context=_build_marker_neighbor_context(
                files,
                species_tree_path=species_tree_path,
                table_path=os.path.join(cfg.outdir, "table_elim_dups"),
                k=max(2, cfg.num_nei or 5),
                num_cpus=cfg.num_cpus,
            ),
            mode=cfg.singles_mode,
        )
        accepted = select_singleton_proposals(
            [proposal for proposal in proposals if proposal.get("singleton_class") == "contamination_candidate"],
            genome_marker_counts=_count_genome_marker_support(files),
            min_markers_per_genome=1,
        )
    accepted_keys = {
        (str(proposal["marker_name"]), str(proposal["leaf_name"]))
        for proposal in accepted
    }
    proposal_lookup = {
        (proposal["marker_name"], proposal["leaf_name"]): proposal
        for proposal in proposals
    }
    proposals_by_marker: dict[str, list[dict]] = {}
    for proposal in proposals:
        proposals_by_marker.setdefault(str(proposal["marker_name"]), []).append(proposal)
    _write_singleton_candidate_table(
        results,
        outdir=cfg.outdir,
        accepted_keys=accepted_keys,
        proposal_lookup=proposal_lookup,
        proposal_keys=proposal_keys,
    )
    for result in sorted(results, key=lambda item: item["filepath"]):
        result = dict(result)
        result["classified_candidates"] = proposals_by_marker.get(str(result["marker_name"]), [])
        _write_singleton_result(
            result,
            species_tree_path=species_tree_path,
            outdir=cfg.outdir,
            accepted_keys=accepted_keys,
        )


def _write_cleaned_seq_worker(args):
    """Worker: write cleaned sequence FASTA for one newick file."""
    filepath, source_dir, output_dir = args
    t = Tree(filepath)
    lst_nodes = [node.name for node in t.traverse("postorder")]
    marker = os.path.basename(filepath).split("_")[3]
    seq_path = os.path.join(source_dir, marker + ".faa")

    with open(seq_path) as f:
        record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    # keep only nodes present in the cleaned tree
    for key in list(record_dict.keys()):
        if key not in lst_nodes:
            del record_dict[key]

    out_path = os.path.join(output_dir, marker + ".faa")
    with open(out_path, "w") as out:
        for k in record_dict:
            SeqIO.write(record_dict[k], out, "fasta")


def write_cleaned_sequences(cfg: Config, use_singles: bool | None = None) -> str:
    """Write cleaned sequence FASTAs based on marker-selection results.

    These cleaned sequence files are then re-aligned so the final tree is built
    from alignments that never contained the removed contaminants.
    """
    if use_singles is None:
        use_singles = cfg.singles

    if use_singles:
        newick_dir = os.path.join(cfg.outdir, "protTrees", "no_singles", "*")
    else:
        newick_dir = os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out", "*")

    source_dir = cfg.extracted_seqs_dir
    output_dir = os.path.join(cfg.outdir, "extracted_final")
    os.makedirs(output_dir, exist_ok=True)

    ls_of_files = glob.glob(newick_dir)
    args = [(f, source_dir, output_dir) for f in ls_of_files]

    map_threaded(_write_cleaned_seq_worker, args, cfg.num_cpus)
    return output_dir
