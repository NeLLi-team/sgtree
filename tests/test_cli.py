import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import DEFAULT, call, patch

from sgtree.cleanup import is_generated_entry
from sgtree.cli import _clean_previous_run, _mark_run_directory, main, parse_args
from sgtree.sgtree_logging import _logfile_path


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
        with (
            patch.object(sys, "argv", ["sgtree", "--help"]),
            patch("sys.stdout", output),
            self.assertRaises(SystemExit) as raised,
        ):
            parse_args()
        self.assertEqual(raised.exception.code, 0)
        help_text = output.getvalue()
        self.assertIn("--singles-mode", help_text)
        self.assertIn("--singles_mode", help_text)
        normalized_help = " ".join(help_text.split())
        self.assertIn("loo_profile", normalized_help)
        self.assertIn("reports evidence without pruning", normalized_help)
        self.assertIn(
            "mode-dependent singleton diagnostics or experimental pruning",
            normalized_help,
        )

    def test_invalid_alignment_method_is_rejected(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--aln", "hmmaling"]
        with (
            patch.object(sys, "argv", argv),
            patch("sys.stderr", StringIO()),
            self.assertRaises(SystemExit) as raised,
        ):
            parse_args()
        self.assertNotEqual(raised.exception.code, 0)

    def test_previous_run_cleanup_removes_only_generated_outputs(self):
        with tempfile.TemporaryDirectory() as outdir:
            output_dir = Path(outdir)
            cfg = self._config_for(outdir)
            generated_dirs = (
                cfg.ani_dir,
                cfg.aligned_dir,
                cfg.tables_dir,
                str(output_dir / "treeouts_protTrees"),
            )
            for directory in generated_dirs:
                Path(directory).mkdir()
                (Path(directory) / "stale").touch()
            for name in (
                "tree.nwk",
                "tree_final.nwk",
                "table_elim_dups",
                "logfile_2024_01_01_00_00_00.txt",
                "tree_round_1.nwk",
                "tree.nwk.iqtree.log",
                "hits.hmmout.del.ls",
            ):
                (output_dir / name).touch()

            analysis_dir = output_dir / "my_analysis"
            analysis_dir.mkdir()
            (analysis_dir / "notes.md").touch()
            (output_dir / "genomes_of_interest.txt").touch()
            (output_dir / "logfile_notes.txt").touch()

            _clean_previous_run(cfg)

            self.assertEqual(
                sorted(path.name for path in output_dir.iterdir()),
                ["genomes_of_interest.txt", "logfile_notes.txt", "my_analysis"],
            )
            self.assertTrue((analysis_dir / "notes.md").exists())

    def test_unmarked_reserved_output_is_preserved_and_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            cfg = self._config_for(outdir)
            Path(cfg.ani_dir).mkdir()
            Path(cfg.ani_cluster_members_path).touch()

            with self.assertRaisesRegex(FileExistsError, "ownership marker"):
                _clean_previous_run(cfg)

            self.assertTrue(Path(cfg.ani_cluster_members_path).exists())

    def test_marker_allows_cleanup_after_a_partial_run(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            cfg = self._config_for(outdir)
            _mark_run_directory(cfg)
            models_dir = Path(outdir) / "models"
            models_dir.mkdir()
            (models_dir / "partial").touch()

            _clean_previous_run(cfg)

            self.assertFalse(models_dir.exists())

    def test_ani_rerun_unlinks_nested_symlinks_without_touching_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_dir = Path(tmpdir)
            outdir = temp_dir / "run"
            outdir.mkdir()
            cfg = self._config_for(str(outdir))
            _mark_run_directory(cfg)
            representatives = Path(cfg.ani_dir) / "query_representatives"
            representatives.mkdir(parents=True)
            target = temp_dir / "genome.fna"
            target.write_text(">contig\nACGT\n")
            (representatives / "genome.fna").symlink_to(target)

            _clean_previous_run(cfg)

            self.assertFalse(Path(cfg.ani_dir).exists())
            self.assertEqual(target.read_text(), ">contig\nACGT\n")

    def test_top_level_generated_symlink_is_preserved_and_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_dir = Path(tmpdir)
            outdir = temp_dir / "run"
            outdir.mkdir()
            cfg = self._config_for(str(outdir))
            _mark_run_directory(cfg)
            target = temp_dir / "external"
            target.mkdir()
            linked = outdir / "models"
            linked.symlink_to(target, target_is_directory=True)

            with self.assertRaisesRegex(FileExistsError, "linked generated path"):
                _clean_previous_run(cfg)

            self.assertTrue(linked.is_symlink())
            self.assertTrue(target.is_dir())

    def test_previous_run_cleanup_leaves_user_suffix_notes(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            cfg = self._config_for(outdir)
            output_dir = Path(outdir)
            (output_dir / "tree.nwk").touch()
            notes = output_dir / "hits.hmmout.notes"
            notes.write_text("keep\n")

            _clean_previous_run(cfg)

            self.assertEqual(notes.read_text(), "keep\n")

    def test_actual_logfile_name_is_generated_but_notes_are_not(self) -> None:
        with tempfile.TemporaryDirectory() as outdir:
            cfg = self._config_for(outdir)
            cfg.start_time = "2026-09-14 12:30:32.123456"
            basename = _logfile_path(cfg).name

        self.assertEqual(basename, "logfile_2026_09_14_12:30:32.txt")
        self.assertTrue(is_generated_entry(basename))
        self.assertTrue(is_generated_entry("logfile_2026_09_14_12_30_32.txt"))
        self.assertFalse(is_generated_entry("logfile_2026_09_14_12:30_32.txt"))
        self.assertFalse(is_generated_entry(f"{basename}.notes"))

    def test_previous_run_cleanup_precedes_ani_and_reference_preparation(self):
        class StopPipelineError(Exception):
            pass

        def fake_prepare_ani(cfg):
            Path(cfg.ani_selected_query_dir).mkdir(parents=True, exist_ok=True)
            with Path(cfg.ani_cluster_members_path).open("w") as handle:
                handle.write("genome_id\n")
            cfg.genomedir = cfg.ani_selected_query_dir

        with tempfile.TemporaryDirectory() as outdir:
            output_dir = Path(outdir)
            (output_dir / "tree.nwk").touch()
            ani_dir = output_dir / "ani"
            ani_dir.mkdir()
            (ani_dir / "stale.tsv").touch()

            argv = [
                "sgtree",
                "input_dir",
                "models.hmm",
                "--save_dir",
                outdir,
                "--ani_cluster",
                "yes",
            ]
            with (
                patch.object(sys, "argv", argv),
                patch("builtins.print"),
                patch(
                    "sgtree.cli.ani_clustering.prepare_ani_cluster_inputs",
                    side_effect=fake_prepare_ani,
                ),
                patch(
                    "sgtree.cli.reference.prepare_reference",
                    side_effect=StopPipelineError,
                ),
                self.assertRaises(StopPipelineError),
            ):
                main()

            # cleanup ran (stale ani/ and the old tree are gone) but did not eat the
            # inputs the ANI step wrote for the rest of the pipeline.
            self.assertFalse((output_dir / "tree.nwk").exists())
            self.assertFalse((ani_dir / "stale.tsv").exists())
            self.assertTrue((ani_dir / "ani_clusters.tsv").exists())
            self.assertTrue((ani_dir / "query_representatives").is_dir())

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
        with (
            patch.object(sys, "argv", argv),
            patch("sys.stderr", StringIO()),
            self.assertRaises(SystemExit),
        ):
            parse_args()

    def test_snp_requires_ani_cluster(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--snp", "yes"]
        with (
            patch.object(sys, "argv", argv),
            self.assertRaisesRegex(ValueError, "--snp requires --ani_cluster yes"),
        ):
            parse_args()

    def test_num_cpus_must_be_positive(self) -> None:
        argv = ["sgtree", "input_dir", "models.hmm", "--num_cpus", "0"]
        with (
            patch.object(sys, "argv", argv),
            self.assertRaisesRegex(ValueError, "num_cpus must be >= 1"),
        ):
            parse_args()

    def test_percent_models_must_be_a_percentage(self) -> None:
        argv = ["sgtree", "input_dir", "models.hmm", "--percent_models", "101"]
        with (
            patch.object(sys, "argv", argv),
            self.assertRaisesRegex(ValueError, "percent_models must be between"),
        ):
            parse_args()

    def test_lflt_must_be_a_percentage(self) -> None:
        argv = ["sgtree", "input_dir", "models.hmm", "--lflt", "-1"]
        with (
            patch.object(sys, "argv", argv),
            self.assertRaisesRegex(ValueError, "lflt must be between"),
        ):
            parse_args()

    def test_max_sdup_rejects_values_below_disabled_sentinel(self) -> None:
        argv = ["sgtree", "input_dir", "models.hmm", "--max_sdup", "-2"]
        with (
            patch.object(sys, "argv", argv),
            self.assertRaisesRegex(ValueError, "max_sdup must be >= -1"),
        ):
            parse_args()

    def test_root_input_paths_are_preserved(self) -> None:
        argv = [
            "sgtree",
            "/",
            "/",
            "--ref",
            "/",
            "--ref_concat",
            "/cache",
        ]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()

        self.assertEqual(cfg.genomedir, "/")
        self.assertEqual(cfg.modeldir, "/")
        self.assertEqual(cfg.ref, "/")
        self.assertEqual(cfg.ref_dir_path(), "/cache/root_root")

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
                patch(
                    "sgtree.cli.reference.prepare_reference", return_value=[]
                ) as prepare_reference,
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
                marker_selection_mocks[
                    "write_cleaned_sequences"
                ].return_value = f"{tmpdir}/cleaned"
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
