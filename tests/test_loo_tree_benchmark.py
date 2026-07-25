import unittest
import warnings

warnings.simplefilter("ignore", SyntaxWarning)

from sgtree.benchmarks import loo_tree_fixtures as benchmark


class LeaveOneOutTreeBenchmarkTests(unittest.TestCase):
    def test_fixture_matrix_and_truth_ids_are_fixed_before_scoring(self):
        first = benchmark.build_fixtures()
        second = benchmark.build_fixtures()

        def signatures(fixtures):
            return [
                (
                    fixture["fixture_id"],
                    tuple(
                        (
                            event["event_id"],
                            event["marker_name"],
                            event["observed_record_id"],
                        )
                        for event in fixture["events"]
                    ),
                    tuple(
                        (marker, tree.write(format=1))
                        for marker, tree in sorted(fixture["trees"].items())
                    ),
                    tuple(
                        (
                            record_id,
                            tuple(
                                tuple(sorted(vote.items()))
                                for vote in votes
                            ),
                        )
                        for record_id, votes in sorted(
                            fixture["votes_by_record"].items()
                        )
                    ),
                )
                for fixture in fixtures
            ]

        self.assertEqual(len(first), 24)
        self.assertEqual(
            {
                (fixture["shape"], fixture["event_class"], fixture["seed"])
                for fixture in first
            },
            {
                (shape, event_class, seed)
                for shape in benchmark.SHAPES
                for event_class in benchmark.EVENT_CLASSES
                for seed in benchmark.SEEDS
            },
        )
        self.assertEqual(signatures(first), signatures(second))
        self.assertTrue(all(fixture["fixture_tier"] == "mechanism" for fixture in first))
        self.assertEqual({len(fixture["genomes"]) for fixture in first}, {10, 12})
        self.assertTrue(all(len(fixture["trees"]) == 8 for fixture in first))
        for shape in benchmark.SHAPES:
            for seed in benchmark.SEEDS:
                paired = [
                    fixture
                    for fixture in first
                    if fixture["shape"] == shape and fixture["seed"] == seed
                ]
                self.assertEqual(len({fixture["genomes"] for fixture in paired}), 1)
                self.assertEqual(len({fixture["markers"] for fixture in paired}), 1)
        event_ids = [
            event["event_id"]
            for fixture in first
            for event in fixture["events"]
        ]
        event_keys = [
            (event["marker_name"], event["observed_record_id"])
            for fixture in first
            for event in fixture["events"]
        ]
        self.assertEqual(len(event_keys), 30)
        self.assertEqual(len(event_keys), len(set(event_keys)))
        self.assertEqual(len(event_ids), len(set(event_ids)))
        for fixture in first:
            leaves_by_marker = {
                marker: {str(leaf.name) for leaf in tree.iter_leaves()}
                for marker, tree in fixture["trees"].items()
            }
            for event in fixture["events"]:
                self.assertIn(
                    event["observed_record_id"],
                    leaves_by_marker[event["marker_name"]],
                )
                if event["is_contaminant"]:
                    self.assertNotIn(
                        event["native_record_id"],
                        leaves_by_marker[event["marker_name"]],
                    )
                if fixture["event_class"] in {
                    "gene_rich_contaminant",
                    "solo_marker_contaminant",
                }:
                    self.assertTrue(event["cmtv_k5_neighbors_preserved"])

    def test_scale_fixture_matrix_and_far_source_truth_are_fixed(self):
        first = benchmark.build_scale_fixtures()
        second = benchmark.build_scale_fixtures()

        def signatures(fixtures):
            return [
                (
                    fixture["fixture_id"],
                    fixture["genomes"],
                    fixture["markers"],
                    tuple(
                        (marker, tree.write(format=1))
                        for marker, tree in sorted(fixture["trees"].items())
                    ),
                    tuple(
                        (event["event_id"], event["observed_record_id"])
                        for event in fixture["events"]
                    ),
                    tuple(
                        (
                            record_id,
                            tuple(vote["assigned_clade"] for vote in votes),
                        )
                        for record_id, votes in sorted(
                            fixture["votes_by_record"].items()
                        )
                    ),
                )
                for fixture in fixtures
            ]

        self.assertEqual(len(first), 8)
        self.assertEqual(signatures(first), signatures(second))
        self.assertEqual(
            {
                (
                    fixture["fixture_tier"],
                    len(fixture["genomes"]),
                    fixture["shape"],
                    fixture["event_class"],
                    fixture["seed"],
                )
                for fixture in first
            },
            {
                ("scale", taxa_count, shape, event_class, benchmark.SCALE_SEED)
                for taxa_count in benchmark.SCALE_TAXA_COUNTS
                for shape in benchmark.SHAPES
                for event_class in benchmark.SCALE_EVENT_CLASSES
            },
        )
        self.assertTrue(all(len(fixture["trees"]) == 8 for fixture in first))
        for taxa_count in benchmark.SCALE_TAXA_COUNTS:
            for shape in benchmark.SHAPES:
                paired = [
                    fixture
                    for fixture in first
                    if len(fixture["genomes"]) == taxa_count
                    and fixture["shape"] == shape
                ]
                self.assertEqual(len({fixture["genomes"] for fixture in paired}), 1)
                self.assertEqual(len({fixture["markers"] for fixture in paired}), 1)

        truth_events = [
            (fixture, event)
            for fixture in first
            for event in fixture["events"]
        ]
        self.assertEqual(len(truth_events), 4)
        self.assertEqual(
            len({event["event_id"] for _fixture, event in truth_events}),
            4,
        )
        for fixture in first:
            if fixture["event_class"] == "scale_clean":
                self.assertEqual(fixture["events"], [])
                continue
            self.assertEqual(len(fixture["events"]), 1)
            event = fixture["events"][0]
            donor = fixture["genomes"][-1]
            recipient_node = fixture["reference_tree"] & fixture["genomes"][2]
            donor_distance = recipient_node.get_distance(
                fixture["reference_tree"] & donor
            )
            self.assertAlmostEqual(
                donor_distance,
                max(
                    recipient_node.get_distance(other)
                    for other in fixture["reference_tree"].iter_leaves()
                    if other is not recipient_node
                ),
            )
            self.assertEqual(event["marker_name"], fixture["markers"][0])
            self.assertTrue(
                event["observed_record_id"].startswith(
                    f"{fixture['genomes'][2]}|"
                )
            )
            self.assertEqual(event["donor_genome"], donor)
            self.assertTrue(event["is_contaminant"])
            self.assertTrue(event["expected_loo"])
            self.assertTrue(event["expected_screen_candidate"])
            self.assertEqual(
                [
                    vote["assigned_clade"]
                    for vote in fixture["votes_by_record"][
                        event["observed_record_id"]
                    ]
                ],
                [donor] * 3,
            )

    def test_full_bounded_screen_passes_all_fixed_expectations(self):
        report = benchmark.run_benchmark()

        self.assertEqual(benchmark.check_benchmark(report), [])
        self.assertEqual(report["fixture_count"], 24)
        self.assertEqual(report["marker_tree_instance_count"], 192)
        self.assertEqual(report["distinct_marker_tree_count"], 78)
        self.assertEqual(report["scale"]["fixture_count"], 8)
        self.assertEqual(report["scale"]["marker_tree_instance_count"], 64)
        self.assertEqual(report["scale"]["distinct_marker_tree_count"], 36)
        self.assertEqual(report["scale"]["taxa_counts"], [50, 100])
        self.assertEqual(report["combined_fixture_count"], 32)
        self.assertEqual(report["combined_marker_tree_instance_count"], 256)
        self.assertEqual(report["combined_distinct_marker_tree_count"], 114)
        self.assertEqual(report["combined_panel_count"], 10)
        self.assertEqual(report["combined_truth_event_count"], 34)
        self.assertEqual(report["benchmark_kind"], "synthetic_tree_mechanism_screen")
        self.assertFalse(report["production_pruning_enabled"])
        self.assertEqual(
            report["metrics"]["base_panel_expectations"]["unit"],
            "paired_shape_seed_panel",
        )
        self.assertEqual(
            report["metrics"]["base_panel_expectations"]["total"],
            6,
        )
        self.assertLessEqual(report["runtime_seconds"], benchmark.MAX_RUNTIME_SECONDS)
        self.assertTrue(
            all(
                metric["low"] is not None and metric["high"] is not None
                for metric in report["metrics"].values()
            )
        )
        self.assertTrue(
            all(
                "low" not in metric and "high" not in metric
                for metric in report["scale"]["metrics"].values()
            )
        )
        self.assertEqual(
            report["scale"]["metrics"]["panel_expectations"]["total"],
            4,
        )
        for metric in (
            "panel_expectations",
            "clean_loo_safety",
            "clean_screen_safety",
            "clean_voter_dispersion",
            "far_source_loo_detection",
            "far_source_screen_detection",
            "far_source_rf_improvement",
        ):
            self.assertEqual(report["scale"]["metrics"][metric]["rate"], 1.0)
        self.assertEqual(
            report["scale"]["metrics"]["cmtv_far_source_detection"]["successes"],
            4,
        )
        self.assertIn(
            "not comparative superiority",
            report["scale"]["comparison_scope"],
        )
        self.assertTrue(
            any(
                "one fixed seed" in limitation
                for limitation in report["scale"]["limitations"]
            )
        )
        self.assertTrue(
            all("cmtv_call_count" in case for case in report["cases"])
        )
        self.assertGreater(
            report["metrics"]["loo_gene_rich_detection"]["successes"],
            report["metrics"]["cmtv_gene_rich_detection"]["successes"],
        )
        self.assertGreaterEqual(
            report["metrics"]["screen_marker_precision"]["rate"],
            0.90,
        )
        self.assertGreaterEqual(
            report["metrics"]["screen_gene_rich_recall"]["rate"],
            0.75,
        )
        clean_cases = [
            case for case in report["cases"] if case["event_class"] == "clean"
        ]
        self.assertTrue(
            all(case["nonzero_voter_dispersion_count"] > 0 for case in clean_cases)
        )
        gene_rich_scores = {
            shape: {
                case["events"][0]["loo_score"]
                for case in report["cases"]
                if case["event_class"] == "gene_rich_contaminant"
                and case["shape"] == shape
            }
            for shape in benchmark.SHAPES
        }
        self.assertTrue(all(len(scores) > 1 for scores in gene_rich_scores.values()))
        for case in report["cases"]:
            for event in case["events"]:
                if event["screen_candidate"]:
                    self.assertTrue(event["rf_guard_pass"])
                    self.assertLess(event["rf_after"], event["rf_before"])
        for case in report["cases"]:
            if case["event_class"] == "safety_controls":
                self.assertEqual(case["truth_contaminant_count"], 2)
                self.assertEqual(case["loo_true_positive_count"], 1)
                self.assertEqual(case["loo_false_positive_count"], 1)
                self.assertEqual(case["screen_candidate_count"], 0)
        for case in report["scale"]["cases"]:
            self.assertEqual(case["fixture_tier"], "scale")
            if case["event_class"] == "scale_clean":
                self.assertEqual(case["loo_call_count"], 0)
                self.assertEqual(case["screen_candidate_count"], 0)
                self.assertGreater(case["nonzero_voter_dispersion_count"], 0)
                continue
            self.assertEqual(case["loo_true_positive_count"], 1)
            self.assertEqual(case["loo_false_positive_count"], 0)
            self.assertEqual(case["screen_true_positive_count"], 1)
            self.assertEqual(case["screen_false_positive_count"], 0)
            self.assertEqual(case["cmtv_true_positive_count"], 1)
            event = case["events"][0]
            self.assertEqual(event["loo_marker_rank"], 1)
            self.assertAlmostEqual(event["loo_target_support"], 0.95)
            self.assertGreater(event["loo_conflict_margin"], 0.0)
            self.assertGreater(event["loo_marker_margin"], 0.0)
            self.assertTrue(event["rf_guard_pass"])
            self.assertLess(
                event["rf_after"],
                event["rf_before"],
            )


if __name__ == "__main__":
    unittest.main()
