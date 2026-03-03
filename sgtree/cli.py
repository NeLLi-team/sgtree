import os
import sys
import glob
import shutil
import datetime
import time
import argparse

from sgtree.config import Config
from sgtree import search, extract, align, duplicates, supermatrix, phylogeny
from sgtree import marker_selection, itol, render, sgtree_logging, cleanup, reference

os.environ["QT_QPA_PLATFORM"] = "offscreen"


def parse_args() -> Config:
    parser = argparse.ArgumentParser(description="SGTree - Species tree from marker gene phylogenies")

    parser.add_argument("genomedir", type=str,
                        help="directory containing .faa proteome files (or a concatenated fasta)")
    parser.add_argument("modeldir", type=str,
                        help="path to marker-set .hmm file (legacy: directory of per-marker .hmm files)")
    parser.add_argument("--ref_concat", type=str, default=None,
                        help="path to store reference directory concat files")
    parser.add_argument("--num_cpus", type=int, default=8,
                        help="number of CPUs to use")
    parser.add_argument("--marker_selection", type=str, default="no",
                        help="run marker selection (yes/no)")
    parser.add_argument("--ref", type=str, default=None,
                        help="reference genomes directory")
    parser.add_argument("--percent_models", type=int, default=10,
                        help="minimum percentage of models a genome must have")
    parser.add_argument("--save_dir", type=str, default=None,
                        help="output directory name")
    parser.add_argument("--singles", type=str, default="no",
                        help="remove singleton markers (yes/no)")
    parser.add_argument("--lflt", type=int, default=0,
                        help="remove sequences shorter than N%% of median length")
    parser.add_argument("--num_nei", type=int, default=15,
                        help="number of neighbors to check")
    parser.add_argument("--aln", type=str, default="hmmalign",
                        help="alignment method: mafft, mafft-linsi, or hmmalign")
    parser.add_argument("--tree_method", type=str, default="fasttree",
                        choices=["fasttree", "iqtree"],
                        help="tree builder for species and marker trees: fasttree or iqtree")
    parser.add_argument("--iqtree_fast", type=str, default="yes",
                        help="when --tree_method iqtree, use IQ-TREE -fast (yes/no)")
    parser.add_argument("--iqtree_model", type=str, default="LG+F+I+G4",
                        help="IQ-TREE model string (used when --tree_method iqtree)")
    parser.add_argument("--hmmsearch_cutoff", type=str, default="cut_ga",
                        choices=["cut_ga", "cut_tc", "cut_nc", "evalue"],
                        help="hmmsearch threshold mode")
    parser.add_argument("--hmmsearch_evalue", type=float, default=1e-5,
                        help="hmmsearch E-value threshold when --hmmsearch_cutoff evalue")
    parser.add_argument("--max_sdup", type=int, default=-1,
                        help="max copies allowed for any single marker per genome (-1 disables)")
    parser.add_argument("--max_dupl", type=float, default=-1.0,
                        help="max fraction of markers allowed in duplicate per genome (-1 disables)")
    parser.add_argument("--is_ref", type=str, default="no",
                        help="internal flag, not for user use")

    args = parser.parse_args()
    start_time = str(datetime.datetime.now())

    # build output directory path
    genomedir = args.genomedir.rstrip("/")
    modeldir = args.modeldir.rstrip("/")

    if args.save_dir:
        outdir = args.save_dir
    else:
        ref_part = args.ref.rstrip("/").split("/")[-1] if args.ref else "no_ref_directory"
        outdir = os.path.join(
            os.getcwd(), "runs", "python",
            f"SG_{genomedir.split('/')[-1]}_{ref_part}_{modeldir.split('/')[-1]}"
            f"_{start_time.replace(':', '-').replace(' ', '_').split('.')[0]}",
        )

    ref_concat = args.ref_concat if args.ref_concat else os.path.join(
        os.getcwd(), "runs", "reference_cache"
    )
    ref_concat = ref_concat.rstrip("/")
    iqtree_fast = str(args.iqtree_fast).strip().lower() in {"yes", "true", "1"}
    if args.max_dupl != -1.0 and not (0.0 <= args.max_dupl <= 1.0):
        raise ValueError("--max_dupl must be between 0 and 1, or -1 to disable")
    if args.hmmsearch_cutoff == "evalue" and args.hmmsearch_evalue <= 0:
        raise ValueError("--hmmsearch_evalue must be > 0 when cutoff mode is evalue")

    return Config(
        genomedir=genomedir,
        modeldir=modeldir,
        outdir=outdir,
        num_cpus=args.num_cpus,
        percent_models=args.percent_models,
        lflt_fraction=float(args.lflt) / 100,
        aln_method=args.aln,
        tree_method=args.tree_method,
        iqtree_fast=iqtree_fast,
        iqtree_model=args.iqtree_model,
        hmmsearch_cutoff=args.hmmsearch_cutoff,
        hmmsearch_evalue=args.hmmsearch_evalue,
        max_sdup=args.max_sdup,
        max_dupl=args.max_dupl,
        ref=args.ref.rstrip("/") if args.ref else None,
        ref_concat=ref_concat,
        marker_selection=args.marker_selection == "yes",
        singles=args.singles == "yes",
        is_ref=args.is_ref == "yes",
        start_time=start_time,
    )


def main():
    cfg = parse_args()
    os.makedirs(cfg.outdir, exist_ok=True)
    cfg.print_banner()

    print(f"-... Reference directory and arguments\n"
          f" Reference directory located at {cfg.ref_concat}\n ...starting sgtree...")

    # handle reference genomes
    ls_refs = reference.prepare_reference(cfg)

    # print arguments
    print(f"arguments:\n"
          f" proteomes, fasta {cfg.genomedir}\n"
          f" models, hmm {cfg.modeldir}\n"
          f" working directory {cfg.outdir}\n"
          f" number of CPUs {cfg.num_cpus}\n"
          f" tree method {cfg.tree_method}\n"
          f" iqtree fast {'yes' if cfg.iqtree_fast else 'no'}\n"
          f" iqtree model {cfg.iqtree_model}\n"
          f" hmmsearch cutoff {cfg.hmmsearch_cutoff}\n"
          f" hmmsearch evalue {cfg.hmmsearch_evalue}\n"
          f" minimum percentage of models {cfg.percent_models}\n"
          f" max single-marker copies {cfg.max_sdup}\n"
          f" max duplicated-marker fraction {cfg.max_dupl}\n"
          f" reference directory {cfg.ref}\n"
          f" --marker_selection {'yes' if cfg.marker_selection else 'no'}\n")
    if cfg.ref:
        print(f"--ref_concat {cfg.ref_dir_path()}\n")
    else:
        print("--ref_concat no reference directory\n")
    print("=" * 80)

    # clean previous runs
    if os.path.exists(os.path.join(cfg.outdir, "tree.nwk")) or \
       os.path.exists(os.path.join(cfg.outdir, "tree_final.nwk")):
        for f in glob.glob(os.path.join(cfg.outdir, "*")):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)

    # write iTOL color file
    itol.write_color_file(cfg)

    timings = {}

    try:
        # Step 1: Concatenate inputs
        cfg.model_count = search.concat_inputs(cfg)

        # Step 2: Run hmmsearch
        t0 = datetime.datetime.now()
        search_time = search.run_hmmsearch(cfg)
        timings["running hmmsearch"] = (t0, search_time)

        # Step 3: Parse results and build working df
        t0 = datetime.datetime.now()
        t_start = time.time()
        finaldf, dict_counts = search.parse_hmmsearch(cfg)
        df, df_fordups = search.build_working_df(cfg, finaldf)
        extract_time = time.time() - t_start

        # Step 4: Extract sequences
        extract.extract_hits(cfg, df)
        extract.write_extracted_sequences(cfg)
        extract_time = time.time() - t_start
        print(f"\nextraction of best hits done - runtime: {extract_time:.1f} seconds")
        print("=" * 80)
        timings["extracting best hits"] = (t0, extract_time)

        # Step 5: Alignment
        t0 = datetime.datetime.now()
        t_start = time.time()
        align.run_alignment(cfg)
        aln_time = time.time() - t_start
        print(f"\nalignment done - runtime: {aln_time:.1f} seconds")
        print("=" * 80 + "\n")
        timings["running alignment"] = (t0, aln_time)

        # Step 6: Eliminate duplicates
        t0 = datetime.datetime.now()
        t_start = time.time()
        print("- ... eliminating duplicates\n")
        duplicates.eliminate_duplicates(cfg, df_fordups)
        dup_time = time.time() - t_start
        print(f"aligned to aln_SpecTree runtime, duplicates eliminated - runtime: {dup_time:.1f} seconds")
        print("=" * 80 + "\n")

        # Step 7: Trim alignments
        t0 = datetime.datetime.now()
        t_start = time.time()
        print("- ...running trimal")
        supermatrix.run_trimal(cfg, cfg.aln_spectree_dir, cfg.trimmed_dir)
        trim_time = time.time() - t_start
        print(f"\ntrimming done - runtime: {trim_time:.1f} seconds")
        print("=" * 80 + "\n")
        timings["running trimal"] = (t0, trim_time)

        # Step 8: Build supermatrix
        t0 = datetime.datetime.now()
        t_start = time.time()
        print("- ...creating supermatrix")
        table_path = os.path.join(cfg.tables_dir, "table_df_concatenated_w_X")
        concat_path = os.path.join(cfg.concat_dir, "concatenated.faa")
        supermatrix.build_supermatrix(cfg.trimmed_dir, cfg.concat_dir, table_path, concat_path)
        concat_time = time.time() - t_start
        print(f"\nsupermatrix created - runtime: {concat_time:.1f} seconds")
        print("=" * 80 + "\n")
        timings["creating supermatrix"] = (t0, concat_time)

        # Step 9: Build species tree
        t0 = datetime.datetime.now()
        t_start = time.time()
        print(f"- ...running {cfg.tree_method}")
        tree_path = os.path.join(cfg.outdir, "tree.nwk")
        phylogeny.run_species_tree(cfg, concat_path, tree_path)
        tree_time = time.time() - t_start
        print(f"\n{cfg.tree_method} done - total runtime: {tree_time:.1f} seconds")
        print("=" * 80 + "\n")
        timings[f"running {cfg.tree_method}"] = (t0, tree_time)

        # iTOL heatmap (for basic run)
        if not cfg.marker_selection:
            itol.write_heatmap(cfg, tree_path, "marker_count.txt")

        # Write logfile
        print(sys.argv[:])
        sgtree_logging.write_logfile(cfg, timings)

        # Cleanup (basic run)
        if not cfg.marker_selection and not cfg.is_ref:
            cleanup.cleanup_basic(cfg.outdir)

    except Exception as e:
        print(f"ERROR: {e.__doc__}\n {e}")
        import traceback
        traceback.print_exc()
        raise

    # Marker selection phase
    if cfg.marker_selection:
        try:
            ms_timings = {}

            # Trim aligned files (with duplicates) for protein trees
            t0 = datetime.datetime.now()
            t_start = time.time()
            print("- ...running trimal (on cleaned proteomes)")
            trimmed_prot_dir = os.path.join(cfg.outdir, "trimmed_protTrees")
            supermatrix.run_trimal_simple(cfg, cfg.aligned_dir, trimmed_prot_dir)
            trim_time = time.time() - t_start
            print(f"\ntrimming done - runtime: {trim_time:.1f} seconds")
            print("=" * 80 + "\n")
            ms_timings["running trimal"] = (t0, trim_time)

            # Build per-marker protein trees
            t0 = datetime.datetime.now()
            t_start = time.time()
            print(f"- ...running {cfg.tree_method}, making protein trees for marker selection:")
            treeout_dir = os.path.join(cfg.outdir, "treeouts_protTrees")
            phylogeny.run_fasttree_per_marker(cfg, trimmed_prot_dir, treeout_dir)
            tree_time = time.time() - t_start
            print(f"\n{cfg.tree_method} done - total runtime: {tree_time:.1f} seconds")
            print("=" * 80 + "\n")

            # RF-distance marker selection
            t_start_ms = time.time()
            print("- ...starting marker selection (Noperm):\n")
            marker_selection.run_noperm(cfg, ls_refs)

            if cfg.singles:
                marker_selection.remove_singles(cfg)

            ms_time = time.time() - t_start_ms
            print(f"Marker selection runtime {ms_time:.1f}")
            print("=" * 80 + "\n")
            ms_timings["Marker selection"] = (t0, ms_time)

            sgtree_logging.append_logfile(cfg, ms_timings)

        except Exception as e:
            print(f"ERROR in marker selection: {e.__doc__}\n {e}")
            import traceback
            traceback.print_exc()
            raise

        try:
            # Write cleaned alignments
            marker_selection.write_cleaned_alignments(cfg)

            # Final trimal
            print("- ...running trimal for final alignment:")
            aligned_final_dir = os.path.join(cfg.outdir, "aligned_final")
            trimmed_final_dir = os.path.join(cfg.outdir, "trimmed_final")
            supermatrix.run_trimal_simple(cfg, aligned_final_dir, trimmed_final_dir)
            print("\ntrimming done")
            print("=" * 80 + "\n")

            # Final supermatrix
            print("- ...creating supermatrix")
            table_path = os.path.join(cfg.tables_dir, "table_df_concatenated_w_X_final")
            concat_final_dir = os.path.join(cfg.outdir, "concat_final")
            concat_final_path = os.path.join(concat_final_dir, "concatenated.faa")
            supermatrix.build_supermatrix(
                trimmed_final_dir, concat_final_dir, table_path, concat_final_path
            )
            print("\nsupermatrix created")
            print("=" * 80 + "\n")

            # Final tree
            print(f"- ...running {cfg.tree_method} for tree_final.nwk")
            tree_final_path = os.path.join(cfg.outdir, "tree_final.nwk")
            phylogeny.run_species_tree(cfg, concat_final_path, tree_final_path)

            # Render tree
            color_file = os.path.join(cfg.outdir, "color.txt")
            render.render_tree(cfg, color_file)

            # iTOL heatmap
            itol.write_heatmap(cfg, tree_final_path, "marker_counts.txt")

            # Cleanup
            cleanup.cleanup_marker_selection(cfg.outdir)

        except Exception as e:
            print(f"ERROR in final tree: {e.__doc__}\n {e}")
            import traceback
            traceback.print_exc()
            raise

    print(f"START: {cfg.start_time} END {datetime.datetime.now()}")


if __name__ == "__main__":
    main()
