import copy
import csv
import tempfile
import unittest
from pathlib import Path

from sgtree.benchmarks import loo_sequence_benchmark as benchmark
from sgtree.benchmarks.loo_tree_fixtures import load_cmtv_rf_weights


class LeaveOneOutSequenceBenchmarkTests(unittest.TestCase):
    def test_manifest_fixes_the_small_held_out_design(self):
        panels = benchmark.build_panels()
        manifest = benchmark.build_manifest(panels)

        self.assertEqual(len(panels), 4)
        self.assertEqual(len(manifest["cases"]), 12)
        self.assertEqual(manifest["development_seeds"], [101, 401, 503])
        self.assertEqual(manifest["held_out_seeds"], [607, 809])
        self.assertEqual(
            manifest["primary_estimand"],
            "truth_topology_far_source_single_copy_contamination",
        )
        self.assertEqual(
            manifest["near_source_role"],
            "identifiability_control_expected_abstention",
        )
        self.assertFalse(manifest["production_pruning_enabled"])
        self.assertEqual(
            manifest["scorer_simplicity_order"],
            ["loo", "cmtv_weighted", "recipient_consensus"],
        )
        self.assertEqual(
            {
                (row["regime"], row["seed"], row["event_class"])
                for row in manifest["cases"]
            },
            {
                (regime, seed, event_class)
                for regime in benchmark.REGIMES
                for seed in benchmark.SEEDS
                for event_class in benchmark.EVENT_CLASSES
            },
        )
        self.assertEqual(
            [row["distance_stratum"] for row in manifest["panels"]].count("near"),
            2,
        )
        self.assertEqual(
            [row["distance_stratum"] for row in manifest["panels"]].count("far"),
            2,
        )

        event_ids = []
        event_keys = []
        for panel, row in zip(panels, manifest["panels"], strict=True):
            self.assertEqual(len(panel["genomes"]), 12)
            self.assertEqual(len(panel["markers"]), 8)
            self.assertNotIn(panel["source_id"], panel["genomes"])
            self.assertNotIn(
                panel["source_id"],
                benchmark._reference_proteomes(panel),
            )
            self.assertEqual(len(row["sequence_design_sha256"]), 64)
            self.assertTrue(row["source_excluded_from_panel"])
            for event in (panel["main_event"], panel["sentinel_event"]):
                self.assertEqual(len(event["candidate_genes"]), 3)
                event_ids.append(event["event_id"])
                event_keys.append((event["marker"], event["observed_record_id"]))
            main = panel["main_event"]
            self.assertEqual(
                main["source_generation"],
                "joint_13_taxon_sister_simulation",
            )
            self.assertEqual(
                main["source_anchor_rule"],
                "truth_topology_extreme_then_blake2_tiebreak",
            )
            self.assertGreater(main["source_tree_edge_distance"], 0)
            self.assertGreater(main["source_anchor_marker_difference_count"], 0)
            self.assertTrue(
                all(
                    difference > 0
                    for difference in main["source_anchor_gene_difference_counts"]
                )
            )

        self.assertEqual(len(event_ids), len(set(event_ids)))
        self.assertEqual(len(event_keys), len(set(event_keys)))
        self.assertEqual(manifest, benchmark.build_manifest(benchmark.build_panels()))

    def test_manifest_hash_covers_inputs_used_for_tree_inference(self):
        panel = benchmark.build_panel("alpha_like", 607)
        original = benchmark.build_manifest([panel])
        changed = copy.deepcopy(panel)
        marker = changed["markers"][0]
        genome = changed["genomes"][0]
        changed["marker_sequences"][marker][genome] = benchmark._mutate(
            changed["marker_sequences"][marker][genome],
            seed_parts=("manifest-test",),
            count=1,
        )
        modified = benchmark.build_manifest([changed])

        self.assertNotEqual(
            original["panels"][0]["sequence_design_sha256"],
            modified["panels"][0]["sequence_design_sha256"],
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            outdir = Path(temp_dir)
            benchmark._write_frozen_manifest(outdir, original)
            with self.assertRaisesRegex(ValueError, "manifest differs"):
                benchmark._write_frozen_manifest(outdir, modified)

    def test_thread_ceiling_fails_before_work_starts(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            outdir = Path(temp_dir)
            with self.assertRaisesRegex(ValueError, "between 1 and 4"):
                benchmark.run_benchmark(outdir, threads=5)
            self.assertEqual(list(outdir.iterdir()), [])

    def test_cmtv_rf_loader_parses_four_columns_into_positive_weights(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "marker_selection_rf_values.txt"
            path.write_text(
                "ProteinID MarkerGene RFdistance Status\n"
                "g1/c1/p1 M1 0.250000 Kept\n"
                "g2/c2/p2 M1 0.250000 Removed\n"
                "g3/c3/p3 M2 1.000000 Kept\n"
            )

            weights = load_cmtv_rf_weights(path)

        self.assertEqual(weights["M1"], 0.75)
        self.assertGreater(weights["M2"], 0.0)

    def test_cmtv_rf_loader_rejects_malformed_or_inconsistent_rows(self):
        fixtures = {
            "malformed": (
                "ProteinID MarkerGene RFdistance Status\n"
                "g1/c1/p1 M1 not-a-number Kept\n"
            ),
            "inconsistent": (
                "ProteinID MarkerGene RFdistance Status\n"
                "g1/c1/p1 M1 0.1 Kept\n"
                "g2/c2/p2 M1 0.2 Removed\n"
            ),
        }
        for label, content in fixtures.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp_dir:
                path = Path(temp_dir) / "marker_selection_rf_values.txt"
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_cmtv_rf_weights(path)

    def test_full_bounded_sequence_benchmark_passes_and_reuses_cache(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            outdir = Path(temp_dir)
            report = benchmark.run_benchmark(outdir, threads=1)

            self.assertEqual(benchmark.check_benchmark(report), [])
            self.assertEqual(report["case_count"], 12)
            self.assertEqual(report["panel_count"], 4)
            self.assertEqual(report["marker_tree_instance_count"], 96)
            self.assertEqual(report["distinct_marker_tree_build_count"], 40)
            self.assertEqual(report["marker_tree_builds_this_run"], 40)
            self.assertLessEqual(
                report["distinct_marker_tree_build_count"],
                benchmark.MAX_MARKER_TREE_BUILDS,
            )
            self.assertLessEqual(report["threads"], benchmark.MAX_THREADS)
            self.assertLessEqual(
                report["runtime_seconds"],
                benchmark.MAX_RUNTIME_SECONDS,
            )
            self.assertFalse(report["production_pruning_enabled"])
            self.assertTrue(report["manifest_written_before_tree_inference"])
            self.assertEqual(
                report["cmtv_baseline"],
                "CMTV-current with positive RF-quality voter weights",
            )
            self.assertEqual(
                report["fair_pipeline_stages"],
                [
                    "raw_call",
                    "common_attachment_and_contig_audit",
                    "candidate_budget",
                    "rf_guard",
                ],
            )
            self.assertEqual(
                report["compared_scorers"],
                ["loo", "cmtv_weighted", "recipient_consensus"],
            )
            self.assertEqual(
                report["metrics"]["loo_gene_rich_true_positive_count"],
                2,
            )
            self.assertEqual(report["metrics"]["loo_truth_false_positive_count"], 0)
            self.assertEqual(
                report["metrics"]["loo_gene_rich_f1"],
                report["metrics"]["cmtv_weighted_gene_rich_f1"],
            )
            self.assertEqual(report["wp1_decision"]["status"], "select_simpler_tied_scorer")
            self.assertEqual(report["wp1_decision"]["selected_scorer"], "loo")
            self.assertTrue(report["wp1_decision"]["empirical_confirmation_ready"])
            self.assertEqual(
                report["wp1_decision"]["tie_break_order"],
                ["loo", "cmtv_weighted", "recipient_consensus"],
            )
            self.assertGreaterEqual(
                report["metrics"]["loo_truth_marker_precision"]["rate"],
                0.90,
            )
            self.assertEqual(report["metrics"]["loo_gene_rich_recall"]["rate"], 0.5)
            self.assertEqual(report["metrics"]["loo_far_source_recall"]["rate"], 1.0)
            self.assertEqual(
                report["metrics"]["loo_near_source_abstention"]["rate"],
                1.0,
            )
            self.assertEqual(
                report["metrics"]["loo_marker_rf_guard_improvement"]["rate"],
                1.0,
            )
            self.assertEqual(
                report["metrics"]["loo_far_source_patristic_improvement"]["rate"],
                1.0,
            )
            self.assertIn(
                "recipient_consensus_truth_marker_precision",
                report["metrics"],
            )
            self.assertIn(
                "cmtv_weighted_unsupported_action_count",
                report["metrics"],
            )
            comparison_path = Path(report["per_event_comparison_path"])
            with comparison_path.open() as handle:
                serialized_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                len(serialized_rows),
                len(report["per_event_comparison_rows"]),
            )
            self.assertTrue(
                all(
                    row["reason"]
                    in {
                        "no_raw_call",
                        "attachment_abstention",
                        "evidence_abstention",
                        "budget_rejected",
                        "rf_guard_rejected",
                        "removed",
                    }
                    for row in report["per_event_comparison_rows"]
                )
            )
            self.assertTrue(
                all(
                    row["raw_call"]
                    and row["contig_gate_pass"]
                    and row["budget_pass"]
                    and row["rf_guard_decision"] == "pruned"
                    for row in report["per_event_comparison_rows"]
                    if row["final_action"] == "removed"
                )
            )

            gene_rich_removals = 0
            for case in report["cases"]:
                self.assertTrue(case["loo_final_taxa_match"])
                self.assertTrue(case["cmtv_weighted_final_taxa_match"])
                self.assertTrue(case["recipient_consensus_final_taxa_match"])
                self.assertEqual(len(case["cmtv_marker_rf_weights"]), 8)
                self.assertTrue(
                    all(
                        weight > 0.0
                        for weight in case["cmtv_marker_rf_weights"].values()
                    )
                )
                self.assertEqual(case["loo_unaffected_record_loss_count"], 0)
                self.assertEqual(case["loo_unaffected_record_losses"], [])
                self.assertEqual(case["loo_unknown_removal_count"], 0)
                self.assertGreaterEqual(case["loo_min_remaining_marker_count"], 1)
                self.assertTrue(case["truth_reference_preserves_sentinel"])
                if case["event_class"] == "clean":
                    self.assertEqual(case["evaluation_role"], "clean_negative_control")
                    self.assertFalse(case["primary_positive"])
                    self.assertEqual(case["loo_removed_count"], 0)
                elif case["event_class"] == "solo_marker_contaminant":
                    self.assertEqual(
                        case["evaluation_role"],
                        "solo_marker_abstention_control",
                    )
                    self.assertFalse(case["primary_positive"])
                    self.assertEqual(case["loo_removed_count"], 0)
                else:
                    events = {event["event_kind"]: event for event in case["events"]}
                    main = events["source_replacement"]
                    expected_primary = case["distance_stratum"] == "far"
                    self.assertEqual(case["primary_positive"], expected_primary)
                    self.assertEqual(
                        case["evaluation_role"],
                        (
                            "primary_far_source_positive"
                            if expected_primary
                            else "near_source_identifiability_control"
                        ),
                    )
                    gene_rich_removals += int(main["loo_removed"])
                    self.assertEqual(
                        main["loo_removed"],
                        case["distance_stratum"] == "far",
                    )
                    if main["loo_removed"]:
                        self.assertEqual(main["loo_class"], "discordant_marker")
                        self.assertGreaterEqual(
                            main["sequence_informative_vote_count"],
                            3,
                        )
                        self.assertTrue(main["contig_gate_pass"])
                    self.assertFalse(events["native_contig_sentinel"]["loo_removed"])
            self.assertEqual(gene_rich_removals, 2)

            cached = benchmark.run_benchmark(outdir, threads=1)
            self.assertEqual(cached["marker_tree_builds_this_run"], 0)
            self.assertEqual(cached["species_tree_builds_this_run"], 0)
            self.assertEqual(cached["metrics"], report["metrics"])
            self.assertEqual(cached["cases"], report["cases"])

            marker_tree = next((outdir / "marker_cache").glob("*/trees/*.nwk"))
            marker_tree.write_text("(wrong_a:0.1,wrong_b:0.1);\n")
            repaired = benchmark.run_benchmark(outdir, threads=1)
            self.assertEqual(repaired["marker_tree_builds_this_run"], 1)
            self.assertEqual(repaired["metrics"], report["metrics"])


if __name__ == "__main__":
    unittest.main()
