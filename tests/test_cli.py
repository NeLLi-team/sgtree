import os
import sys
import tempfile
import unittest
from io import StringIO
from unittest.mock import DEFAULT, call, patch

from sgtree.cli import _clean_previous_run, main, parse_args


class CliTests(unittest.TestCase):
    def _config_for(self, outdir):
        argv = ["sgtree", "input_dir", "models.hmm", "--save_dir", outdir]
        with patch.object(sys, "argv", argv):
            return parse_args()

    def test_alignment_defaults_to_hmmalign(self):
        argv = ["sgtree", "input_dir", "models.hmm"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertEqual(cfg.aln_method, "hmmalign")

    def test_required_inputs_accept_flag_aliases(self):
        argv = ["sgtree", "--genomedir", "input_dir", "--modeldir", "models.hmm"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertEqual(cfg.genomedir, "input_dir")
        self.assertEqual(cfg.modeldir, "models.hmm")

    def test_singles_mode_accepts_hyphenated_preferred_alias(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--singles-mode", "gcp"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertEqual(cfg.singles_mode, "gcp")

    def test_singles_mode_accepts_underscore_legacy_alias(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--singles_mode", "gcp"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertEqual(cfg.singles_mode, "gcp")

    def test_singles_mode_accepts_report_only_loo_profile(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--singles-mode", "loo_profile"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertEqual(cfg.singles_mode, "loo_profile")

    def test_help_lists_preferred_and_legacy_singles_mode_spellings(self):
        output = StringIO()
        with patch.object(sys, "argv", ["sgtree", "--help"]), patch("sys.stdout", output):
            with self.assertRaises(SystemExit) as raised:
                parse_args()
        self.assertEqual(raised.exception.code, 0)
        help_text = output.getvalue()
        self.assertIn("--singles-mode", help_text)
        self.assertIn("--singles_mode", help_text)
        normalized_help = " ".join(help_text.split())
        self.assertIn("loo_profile", normalized_help)
        self.assertIn("reports evidence without pruning", normalized_help)

    def test_invalid_alignment_method_is_rejected(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--aln", "hmmaling"]
        with patch.object(sys, "argv", argv), patch("sys.stderr", StringIO()):
            with self.assertRaises(SystemExit) as raised:
                parse_args()
        self.assertNotEqual(raised.exception.code, 0)

    def test_previous_run_cleanup_removes_only_generated_outputs(self):
        with tempfile.TemporaryDirectory() as outdir:
            cfg = self._config_for(outdir)
            generated_dirs = (
                cfg.ani_dir,
                cfg.aligned_dir,
                cfg.tables_dir,
                os.path.join(outdir, "treeouts_protTrees"),
            )
            for directory in generated_dirs:
                os.makedirs(directory)
                open(os.path.join(directory, "stale"), "w").close()
            for name in ("tree.nwk", "tree_final.nwk", "table_elim_dups",
                         "logfile_2024_01_01_00_00_00.txt", "tree_round_1.nwk",
                         "tree.nwk.iqtree.log", "hits.hmmout.del.ls"):
                open(os.path.join(outdir, name), "w").close()

            os.makedirs(os.path.join(outdir, "my_analysis"))
            open(os.path.join(outdir, "my_analysis", "notes.md"), "w").close()
            open(os.path.join(outdir, "genomes_of_interest.txt"), "w").close()
            open(os.path.join(outdir, "logfile_notes.txt"), "w").close()

            _clean_previous_run(cfg)

            self.assertEqual(
                sorted(os.listdir(outdir)),
                ["genomes_of_interest.txt", "logfile_notes.txt", "my_analysis"],
            )
            self.assertTrue(os.path.exists(os.path.join(outdir, "my_analysis", "notes.md")))

    def test_previous_run_cleanup_is_skipped_without_a_previous_tree(self):
        with tempfile.TemporaryDirectory() as outdir:
            cfg = self._config_for(outdir)
            os.makedirs(cfg.ani_dir)
            open(cfg.ani_cluster_members_path, "w").close()

            _clean_previous_run(cfg)

            self.assertTrue(os.path.exists(cfg.ani_cluster_members_path))

    def test_previous_run_cleanup_precedes_ani_and_reference_preparation(self):
        class StopPipeline(Exception):
            pass

        def fake_prepare_ani(cfg):
            os.makedirs(cfg.ani_selected_query_dir, exist_ok=True)
            with open(cfg.ani_cluster_members_path, "w") as handle:
                handle.write("genome_id\n")
            cfg.genomedir = cfg.ani_selected_query_dir

        with tempfile.TemporaryDirectory() as outdir:
            open(os.path.join(outdir, "tree.nwk"), "w").close()
            os.makedirs(os.path.join(outdir, "ani"))
            open(os.path.join(outdir, "ani", "stale.tsv"), "w").close()

            argv = [
                "sgtree", "input_dir", "models.hmm",
                "--save_dir", outdir,
                "--ani_cluster", "yes",
            ]
            with (
                patch.object(sys, "argv", argv),
                patch("builtins.print"),
                patch(
                    "sgtree.cli.ani_clustering.prepare_ani_cluster_inputs",
                    side_effect=fake_prepare_ani,
                ),
                patch("sgtree.cli.reference.prepare_reference", side_effect=StopPipeline),
            ):
                with self.assertRaises(StopPipeline):
                    main()

            # cleanup ran (stale ani/ and the old tree are gone) but did not eat the
            # inputs the ANI step wrote for the rest of the pipeline.
            self.assertFalse(os.path.exists(os.path.join(outdir, "tree.nwk")))
            self.assertFalse(os.path.exists(os.path.join(outdir, "ani", "stale.tsv")))
            self.assertTrue(os.path.exists(os.path.join(outdir, "ani", "ani_clusters.tsv")))
            self.assertTrue(os.path.isdir(os.path.join(outdir, "ani", "query_representatives")))

    def test_snp_defaults_to_disabled(self):
        argv = ["sgtree", "input_dir", "models.hmm"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertFalse(cfg.snp)
        self.assertFalse(cfg.ani_cluster)

    def test_required_input_conflict_between_positional_and_flag_fails(self):
        argv = [
            "sgtree",
            "input_dir",
            "models.hmm",
            "--genomedir",
            "other_input_dir",
        ]
        with patch.object(sys, "argv", argv), patch("sys.stderr", StringIO()):
            with self.assertRaises(SystemExit):
                parse_args()

    def test_snp_requires_ani_cluster(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--snp", "yes"]
        with patch.object(sys, "argv", argv):
            with self.assertRaisesRegex(ValueError, "--snp requires --ani_cluster yes"):
                parse_args()

    def test_snp_enabled_with_ani_cluster(self):
        argv = [
            "sgtree",
            "input_dir",
            "models.hmm",
            "--ani_cluster",
            "yes",
            "--snp",
            "yes",
        ]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertTrue(cfg.ani_cluster)
        self.assertTrue(cfg.snp)

    def test_early_global_round_convergence_runs_singletons_once_then_stops(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = f"{tmpdir}/run"
            argv = [
                "sgtree",
                "input_dir",
                "models.hmm",
                "--save_dir",
                outdir,
                "--marker_selection",
                "yes",
                "--singles",
                "yes",
                "--selection_global_rounds",
                "4",
                "--keep_intermediates",
                "yes",
            ]
            kept = {("GenomeA", "MarkerX")}

            with (
                patch.object(sys, "argv", argv),
                patch("builtins.print"),
                patch("sgtree.cli.reference.prepare_reference", return_value=[]) as prepare_reference,
                patch.multiple(
                    "sgtree.cli.search",
                    concat_inputs=DEFAULT,
                    run_hmmsearch=DEFAULT,
                    parse_hmmsearch=DEFAULT,
                    build_working_df=DEFAULT,
                ) as search_mocks,
                patch.multiple(
                    "sgtree.cli.extract",
                    extract_hits=DEFAULT,
                    write_extracted_sequences=DEFAULT,
                ),
                patch("sgtree.cli.align.run_alignment"),
                patch("sgtree.cli.duplicates.eliminate_duplicates"),
                patch.multiple(
                    "sgtree.cli.supermatrix",
                    run_trimal=DEFAULT,
                    run_trimal_simple=DEFAULT,
                    build_supermatrix=DEFAULT,
                ),
                patch.multiple(
                    "sgtree.cli.phylogeny",
                    run_species_tree=DEFAULT,
                    run_fasttree_per_marker=DEFAULT,
                ) as phylogeny_mocks,
                patch.multiple(
                    "sgtree.cli.sgtree_logging",
                    write_logfile=DEFAULT,
                    append_logfile=DEFAULT,
                ),
                patch("sgtree.cli.shutil.copyfile"),
                patch.multiple(
                    "sgtree.marker_selection",
                    run_noperm=DEFAULT,
                    remove_singles=DEFAULT,
                    write_cleaned_sequences=DEFAULT,
                ) as marker_selection_mocks,
            ):
                search_mocks["concat_inputs"].return_value = 1
                search_mocks["run_hmmsearch"].return_value = 0.0
                search_mocks["parse_hmmsearch"].return_value = (object(), {})
                search_mocks["build_working_df"].return_value = (object(), object())
                marker_selection_mocks["run_noperm"].side_effect = [kept, kept]
                marker_selection_mocks["write_cleaned_sequences"].return_value = f"{tmpdir}/cleaned"
                main()

            cfg = prepare_reference.call_args.args[0]
            self.assertEqual(marker_selection_mocks["run_noperm"].call_count, 2)
            marker_selection_mocks["remove_singles"].assert_called_once_with(
                cfg,
                species_tree_path=f"{outdir}/tree_final.nwk",
            )
            self.assertEqual(
                marker_selection_mocks["write_cleaned_sequences"].call_args_list,
                [
                    call(cfg, use_singles=False),
                    call(cfg, use_singles=True),
                ],
            )
            self.assertEqual(phylogeny_mocks["run_species_tree"].call_count, 3)


if __name__ == "__main__":
    unittest.main()
