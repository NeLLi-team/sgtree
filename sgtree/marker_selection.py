import os
import glob
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO
from ete3 import Tree

from sgtree.config import Config


def _get_ascore(identifier, table_path):
    """Get the score for a given sequence from the elim_dups table."""
    dfa = pd.read_csv(table_path)
    dfa = dfa.set_index(dfa["savedname"])
    row = dfa.loc[identifier.replace("|", "/")]
    return row["savedname"] + ":" + str(row.iloc[8])


def _best_score(scored_list):
    """Return the entry with the highest score from a list of 'name:score' strings."""
    scores = [float(s.split(":")[1]) for s in scored_list]
    best_idx = scores.index(max(scores))
    return scored_list[best_idx]


def _removekey(d, key):
    """Return a copy of dict d without the given key."""
    return {k: v for k, v in d.items() if k != key}


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


def _process_tree_worker(args):
    """Worker: RF-distance based duplicate resolution for one marker tree."""
    filepath, table_path, species_tree_path, outdir, ls_refs, rf_outfile = args
    t = Tree(filepath)
    lst_nodes = [node.name for node in next(t.copy().traverse())]
    marker_name = os.path.basename(filepath).split(".")[0]

    # find duplicates
    seen = {}
    for x in lst_nodes:
        genome = x.split("|")[0]
        if genome not in seen:
            seen[genome] = x
        else:
            seen[genome] = seen[genome] + "," + x
    dups = {k: v.split(",") for k, v in seen.items()}

    # score all entries
    for key in dups:
        dups[key] = [_get_ascore(v, table_path) for v in dups[key]]

    speciestree = Tree(species_tree_path)
    ls_best_nodes = []
    bad_nodes = []

    for key, value in dups.items():
        # skip reference genomes
        if ls_refs is not None and key.split("|")[0] + ".faa" in ls_refs:
            continue
        if len(value) == 1:
            continue

        best_each = None
        best_protein = None
        rf_results = []
        candidates = []

        for each in value:
            ls_nodes_local = [each]
            dcopy = _removekey(dups, key)
            for k, v in dcopy.items():
                bs = _best_score(v)
                if bs is not None:
                    ls_nodes_local.append(bs)

            ls_nodes_local = [e for e in ls_nodes_local if e is not None]
            alt = [x.split(":")[0].replace("/", "|") for x in ls_nodes_local]

            speciestree_copy = speciestree.copy()
            speciestree_copy.prune([n.split("|")[0] for n in alt])

            t_prot = t.copy()
            t_prot.prune(alt)
            t_protcopy = t_prot.copy()

            for node in next(t_protcopy.traverse()):
                node.name = node.name.split("|")[0]

            rf, maxrf, _, _, _, _, _ = speciestree_copy.robinson_foulds(
                t_protcopy, unrooted_trees=True
            )
            maxrf = maxrf + 0.0001
            rf_dist = rf / maxrf

            rf_results.append((each.split(":")[0], rf_dist))

            candidates.append((rf_dist, each))

        if candidates:
            # Deterministic tie-breaker: lowest RF, then lexicographically smallest protein id.
            _, best_each = min(candidates, key=lambda item: (item[0], item[1].split(":")[0]))
            best_protein = best_each.split(":")[0]

        # write RF results
        with open(rf_outfile, "a") as f:
            for protein, rf_d in rf_results:
                status = "Kept" if protein == best_protein else "Removed"
                f.write(f"{protein} {marker_name} {rf_d:.6f} {status}\n")

        if best_each is not None:
            ls_best_nodes.append(best_each)
        bad_nodes.extend(
            [n.split(":")[0].replace("/", "|") for n in value if n != best_each]
        )

    bad_nodes = [s.split(":")[0].replace("/", "|") for s in bad_nodes]

    with open(os.path.join(outdir, "removed", marker_name), "w") as f:
        for item in bad_nodes:
            f.write(f"{item}\n")
        f.write(f"{len(bad_nodes)} {len(lst_nodes)}\n{'*' * 80}\n")

    t_final = t.copy()
    cleaned_nodes = [n for n in lst_nodes if n not in bad_nodes]
    t_final.prune(cleaned_nodes)
    t_final.write(
        format=1,
        outfile=os.path.join(
            outdir, "protTrees", "no_duplicates", "out",
            f"_no_dups_{marker_name}_.nw",
        ),
    )


def run_noperm(cfg: Config, ls_refs: list[str] | None):
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

    species_tree_path = os.path.join(cfg.outdir, "tree.nwk")
    ls_of_files = glob.glob(os.path.join(treeout_dir, "*"))
    args = [
        (f, table_path, species_tree_path, cfg.outdir, ls_refs, rf_outfile)
        for f in ls_of_files
    ]

    _map_with_fallback(_process_tree_worker, args, cfg.num_cpus)


def _score_func(ls_mod, ls_in):
    """Score neighbor list similarity between protein tree and species tree."""
    score = 0
    for i, init in enumerate(ls_in):
        if init in ls_mod:
            score += len(ls_in) - abs(ls_in.index(init) - ls_mod.index(init)) - i
    return score


def _remove_singles_worker(args):
    """Worker: remove singleton markers for one file."""
    filepath, species_tree_path, outdir = args
    tf = Tree(filepath)
    td = Tree(filepath)

    # RF comparison is genome-level; deduplicate leaves that collapse to the
    # same genome id (can happen for reference genomes).
    seen_genomes = set()
    keep_leaves = []
    for leaf in td.iter_leaves():
        genome = leaf.name.split("|")[0]
        if genome not in seen_genomes:
            seen_genomes.add(genome)
            keep_leaves.append(leaf.name)
    td.prune(keep_leaves)

    td_leaves = list(td.iter_leaves())
    for node in td_leaves:
        node.name = node.name.split("|")[0]

    lst_nodes = [node for node in next(tf.copy().traverse())]
    ti = Tree(species_tree_path)
    ti.prune([n.name for n in td_leaves])

    rf, maxrf, *_ = ti.robinson_foulds(td, unrooted_trees=True)
    maxrf = maxrf + 0.0001
    rdist = rf / maxrf
    num_nei = round(len(lst_nodes) * (1 - rdist))
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
    marker_name = os.path.basename(filepath).split(".")[0].split("_")[-2]

    for node in lst_nodes:
        ls_model = [n.split("|")[0] for n in dict_neighbors[node.name]]
        ls_init = dict_neighborsi[node.name.split("|")[0]]
        score = _score_func(ls_model, ls_init)
        if score > cutoff:
            flagged.append(node)
        else:
            removed_count += 1
            with open(os.path.join(outdir, "removed", marker_name), "a") as f:
                f.write(f"\nsin{node.name}\n")

    with open(os.path.join(outdir, "removed", marker_name), "a") as f:
        f.write(f"\n{removed_count}/{len(lst_nodes)} {rdist}")

    flagged_names = [x.name for x in flagged]
    t_prot = tf.copy()
    t_prot.prune(flagged_names)
    t_prot.write(
        format=1,
        outfile=os.path.join(
            outdir, "protTrees", "no_singles",
            os.path.basename(filepath).split(".")[0] + ".nw",
        ),
    )


def remove_singles(cfg: Config):
    """Remove singleton markers poorly placed by comparing neighbor topology."""
    files = glob.glob(os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out", "*"))
    species_tree_path = os.path.join(cfg.outdir, "tree.nwk")
    args = [(f, species_tree_path, cfg.outdir) for f in files]

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
    with open(out_path, "a") as out:
        for k in record_dict:
            SeqIO.write(record_dict[k], out, "fasta")


def write_cleaned_alignments(cfg: Config):
    """Write cleaned alignments to aligned_final/ based on marker selection results."""
    if cfg.singles:
        newick_dir = os.path.join(cfg.outdir, "protTrees", "no_singles", "*")
    else:
        newick_dir = os.path.join(cfg.outdir, "protTrees", "no_duplicates", "out", "*")

    aligned_dir = cfg.aligned_dir
    aligned_final_dir = os.path.join(cfg.outdir, "aligned_final")
    os.makedirs(aligned_final_dir, exist_ok=True)

    ls_of_files = glob.glob(newick_dir)
    args = [(f, aligned_dir, aligned_final_dir) for f in ls_of_files]

    _map_with_fallback(_write_cleaned_worker, args, cfg.num_cpus)
