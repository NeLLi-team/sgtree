import os
import glob
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO
from ete3 import Tree

from sgtree.config import Config


SCORE_COLUMNS = ("score_bits", "7")


def choose_best_candidate(
    candidates: list[dict],
    preferred_protein_id: str | None = None,
) -> dict:
    """Rank RF-selection candidates deterministically.

    Preference order:
    1. lower RF distance
    2. more informative splits in the RF comparison
    3. higher HMM bitscore
    4. lexicographically smaller protein id for stability

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


def _evaluate_candidate(
    marker_tree: Tree,
    species_tree: Tree,
    dups: dict[str, list[str]],
    selected: dict[str, str],
    genome: str,
    scored_entry: str,
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
    rf_distance = rf / maxrf if maxrf else 0.0
    return {
        "genome": genome,
        "protein_id": protein_id,
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
                _evaluate_candidate(marker_tree, species_tree, dups, fixed_selected, genome, entry)
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
                _evaluate_candidate(marker_tree, species_tree, dups, selected, genome, entry)
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
    score_table, score_col = _load_score_table(table_path)
    marker_tree = Tree(marker_tree_path)
    lst_nodes = [leaf.name for leaf in marker_tree.iter_leaves()]
    dups = _build_duplicate_map(lst_nodes, score_table, score_col)
    species_tree = Tree(species_tree_path)

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
            _evaluate_candidate(marker_tree, species_tree, dups, selected, genome, entry)
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
    """Worker: RF-distance based duplicate resolution for one marker tree."""
    (
        filepath,
        table_path,
        species_tree_path,
        outdir,
        ls_refs,
        rf_outfile,
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

    with open(rf_outfile, "a") as f:
        for record in records:
            f.write(
                f"{record['protein_id']} {record['marker']} "
                f"{record['rf_distance']:.6f} {record['status']}\n"
            )

    t = Tree(filepath)
    t_final = t.copy()
    t_final.prune(cleaned_nodes)
    t_final.write(
        format=1,
        outfile=os.path.join(
            outdir, "protTrees", "no_duplicates", "out",
            f"_no_dups_{marker_name}_.nw",
        ),
    )


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
    with open(rf_outfile, "w") as f:
        f.write("ProteinID MarkerGene RFdistance Status\n")

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
            rf_outfile,
            getattr(cfg, "selection_mode", "coordinate"),
            getattr(cfg, "selection_max_rounds", 5),
            getattr(cfg, "lock_references", False),
            initial_kept,
        )
        for f in ls_of_files
    ]

    _map_with_fallback(_process_tree_worker, args, cfg.num_cpus)
    return _load_kept_assignments(rf_outfile)


def _score_func(ls_mod, ls_in):
    """Score neighbor list similarity between protein tree and species tree."""
    score = 0
    for i, init in enumerate(ls_in):
        if init in ls_mod:
            score += len(ls_in) - abs(ls_in.index(init) - ls_mod.index(init)) - i
    return score


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
    leaf = next(leaf for leaf in tree.iter_leaves() if leaf.name == leaf_name)
    neighbors = []
    for other in tree.iter_leaves():
        if other.name == leaf_name:
            continue
        neighbors.append((leaf.get_distance(other), other.name.split("|")[0]))
    neighbors.sort(key=lambda item: (item[0], item[1]))
    seen = []
    for _dist, genome in neighbors:
        if genome not in seen:
            seen.append(genome)
        if len(seen) >= k:
            break
    return seen


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


def _score_singleton_leaf(
    species_tree: Tree,
    working_tree: Tree,
    leaf_name: str,
    mode: str,
    *,
    k: int,
    score_table: pd.DataFrame | None = None,
    score_col: str | None = None,
) -> tuple[float, Tree]:
    delta_rf, candidate = _leaf_rf_delta(species_tree, working_tree, leaf_name)
    overlap = _leaf_neighbor_overlap(species_tree, working_tree, leaf_name, k)
    branch_outlier = _branch_length_outlier(working_tree, leaf_name)
    bitscore_outlier = _leaf_bitscore_outlier(working_tree, leaf_name, score_table, score_col)

    if mode == "delta_rf":
        score = delta_rf
    elif mode == "backbone":
        score = delta_rf + max(0.0, 0.5 - overlap)
    elif mode == "ensemble":
        score = (3.0 * delta_rf) + max(0.0, 0.5 - overlap) + (0.25 * branch_outlier) + (0.25 * bitscore_outlier)
    else:
        raise ValueError(f"Unknown singleton mode: {mode}")
    return score, candidate


def prune_singletons(
    species_tree: Tree,
    working_tree: Tree,
    *,
    mode: str,
    k: int,
    score_table: pd.DataFrame | None = None,
    score_col: str | None = None,
) -> Tree:
    current = working_tree.copy()
    max_pruned = max(1, len(list(current.iter_leaves())) // 5)
    removed = 0
    while len(list(current.iter_leaves())) > 4 and removed < max_pruned:
        best_leaf = None
        best_score = 0.0
        best_candidate = None
        for leaf in current.iter_leaves():
            score, candidate = _score_singleton_leaf(
                species_tree,
                current,
                leaf.name,
                mode,
                k=k,
                score_table=score_table,
                score_col=score_col,
            )
            if score > best_score:
                best_leaf = leaf.name
                best_score = score
                best_candidate = candidate
        if best_leaf is None or best_candidate is None or best_score <= 0:
            break
        current = best_candidate
        removed += 1
    return current


def _remove_singles_worker(args):
    """Worker: remove singleton markers for one file."""
    filepath, species_tree_path, outdir, num_nei_override, singles_min_rfdist, singles_mode, table_path = args
    tf = Tree(filepath)
    td = _tree_to_genome_level(tf)
    td_leaves = list(td.iter_leaves())
    lst_nodes = [node for node in next(tf.copy().traverse())]
    ti = Tree(species_tree_path)
    ti.prune([n.name for n in td_leaves])

    rf, maxrf, *_ = ti.robinson_foulds(td, unrooted_trees=True)
    maxrf = maxrf + 0.0001
    rdist = rf / maxrf
    if rdist < singles_min_rfdist:
        tf.write(
            format=1,
            outfile=os.path.join(
                outdir, "protTrees", "no_singles",
                os.path.basename(filepath).split(".")[0] + ".nw",
            ),
        )
        return
    if singles_mode in {"delta_rf", "backbone", "ensemble"}:
        score_table = None
        score_col = None
        if singles_mode == "ensemble":
            score_table, score_col = _load_score_table(table_path)
        if num_nei_override > 0:
            num_nei = min(num_nei_override, max(1, len(lst_nodes) - 1))
        else:
            num_nei = max(2, min(5, len(lst_nodes) - 1))
        chosen_tree = prune_singletons(
            species_tree=ti,
            working_tree=tf,
            mode=singles_mode,
            k=num_nei,
            score_table=score_table,
            score_col=score_col,
        )
        choose_tree_by_rf(ti, tf, chosen_tree).write(
            format=1,
            outfile=os.path.join(
                outdir, "protTrees", "no_singles",
                os.path.basename(filepath).split(".")[0] + ".nw",
            ),
        )
        return
    max_neighbors = max(1, len(lst_nodes) - 1)
    if num_nei_override > 0:
        num_nei = min(num_nei_override, max_neighbors)
    else:
        num_nei = max(1, round(len(lst_nodes) * (1 - rdist)))
    total_score = num_nei ** 2
    cutoff = total_score / 15

    # get neighbors from protein tree
    dict_neighbors = {}
    for node in lst_nodes:
        ori_leaf = node
        while len(node.get_leaves()) - 1 < num_nei:
            node = node.up
        dict_neighbors[ori_leaf.name] = {}
        for leaf in node.get_leaves():
            if leaf.name != ori_leaf.name:
                dict_neighbors[ori_leaf.name][leaf.name] = ori_leaf.get_distance(leaf.name)
        dict_neighbors[ori_leaf.name] = sorted(
            dict_neighbors[ori_leaf.name],
            key=lambda leaf: (dict_neighbors[ori_leaf.name][leaf], leaf),
        )[:num_nei]

    # get neighbors from species tree
    lst_nodesi = [node for node in next(ti.copy().traverse())]
    dict_neighborsi = {}
    for nodei in lst_nodesi:
        ori_leafi = nodei
        while len(nodei.get_leaves()) - 1 < num_nei:
            nodei = nodei.up
        dict_neighborsi[ori_leafi.name] = {}
        for leaf in nodei.get_leaves():
            if leaf.name != ori_leafi.name:
                dict_neighborsi[ori_leafi.name][leaf.name] = ori_leafi.get_distance(leaf.name)
        dict_neighborsi[ori_leafi.name] = sorted(
            dict_neighborsi[ori_leafi.name],
            key=lambda leaf: (dict_neighborsi[ori_leafi.name][leaf], leaf),
        )[:num_nei]

    # score and filter
    flagged = []
    removed_count = 0
    borderline_count = 0
    marker_name = os.path.basename(filepath).split(".")[0].split("_")[-2]

    for node in lst_nodes:
        ls_model = [n.split("|")[0] for n in dict_neighbors[node.name]]
        ls_init = dict_neighborsi[node.name.split("|")[0]]
        score = _score_func(ls_model, ls_init)
        margin = score - cutoff
        if score > cutoff:
            flagged.append(node)
            if abs(margin) <= 1:
                borderline_count += 1
        else:
            removed_count += 1
            with open(os.path.join(outdir, "removed", marker_name), "a") as f:
                f.write(
                    f"\nsin{node.name}\tscore={score:.3f}\tcutoff={cutoff:.3f}\tmargin={margin:.3f}\n"
                )

    with open(os.path.join(outdir, "removed", marker_name), "a") as f:
        f.write(
            f"\n{removed_count}/{len(lst_nodes)} {rdist} "
            f"num_nei={num_nei} cutoff={cutoff:.3f} borderline_kept={borderline_count}"
        )

    flagged_names = [x.name for x in flagged]
    candidate_tree = tf.copy()
    candidate_tree.prune(flagged_names)
    chosen_tree = choose_tree_by_rf(ti, tf, candidate_tree)
    chosen_tree.write(
        format=1,
        outfile=os.path.join(
            outdir, "protTrees", "no_singles",
            os.path.basename(filepath).split(".")[0] + ".nw",
        ),
    )


def remove_singles(cfg: Config, species_tree_path: str | None = None):
    """Remove singleton markers poorly placed by comparing neighbor topology."""
    files = glob.glob(os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out", "*"))
    if species_tree_path is None:
        species_tree_path = os.path.join(cfg.outdir, "tree.nwk")
    args = [
        (
            f,
            species_tree_path,
            cfg.outdir,
            cfg.num_nei,
            cfg.singles_min_rfdist,
            cfg.singles_mode,
            os.path.join(cfg.outdir, "table_elim_dups"),
        )
        for f in files
    ]

    _map_with_fallback(_remove_singles_worker, args, cfg.num_cpus)


def _write_cleaned_worker(args):
    """Worker: write cleaned alignment for one newick file."""
    filepath, aligned_dir, aligned_final_dir = args
    t = Tree(filepath)
    lst_nodes = [node.name for node in t.traverse("postorder")]
    marker = os.path.basename(filepath).split("_")[3]
    aln_path = os.path.join(aligned_dir, marker + ".faa")

    with open(aln_path) as f:
        record_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    # keep only nodes present in the cleaned tree
    for key in list(record_dict.keys()):
        if key not in lst_nodes:
            del record_dict[key]

    out_path = os.path.join(aligned_final_dir, marker + ".faa")
    with open(out_path, "w") as out:
        for k in record_dict:
            SeqIO.write(record_dict[k], out, "fasta")


def write_cleaned_alignments(cfg: Config, use_singles: bool | None = None):
    """Write cleaned alignments to aligned_final/ based on marker selection results."""
    if use_singles is None:
        use_singles = cfg.singles

    if use_singles:
        newick_dir = os.path.join(cfg.outdir, "protTrees", "no_singles", "*")
    else:
        newick_dir = os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out", "*")

    aligned_dir = cfg.aligned_dir
    aligned_final_dir = os.path.join(cfg.outdir, "aligned_final")
    os.makedirs(aligned_final_dir, exist_ok=True)

    ls_of_files = glob.glob(newick_dir)
    args = [(f, aligned_dir, aligned_final_dir) for f in ls_of_files]

    _map_with_fallback(_write_cleaned_worker, args, cfg.num_cpus)
