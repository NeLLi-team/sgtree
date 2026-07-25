import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

warnings.simplefilter("ignore", SyntaxWarning)

import pandas as pd
from sgtree import marker_selection

Tree = marker_selection.Tree


class MarkerSelectionTests(unittest.TestCase):
    def _synthetic_gcp_panel(self, marker_count, genome_count=5):
        proposals = []
        for genome_index in range(genome_count):
            for marker_index in range(marker_count):
                is_outlier = genome_index == 0 and marker_index == 0
                proposal = {
                    "marker_name": f"Marker{marker_index:02d}",
                    "genome": f"Genome{genome_index:02d}",
                    "leaf_name": f"Genome{genome_index:02d}|contig1|gene{marker_index:02d}",
                    "contig_id": "contig1",
                    "taxa_count": genome_count,
                }
                for feature in marker_selection.GCP_KEY_FEATURES:
                    proposal[feature] = 10.0 if is_outlier else 0.1
                proposal["anchor_knn_agreement"] = 0.0 if is_outlier else 0.8
                proposals.append(proposal)
        return proposals

    def test_choose_best_candidate_prefers_bitscore_when_rf_is_uninformative(self):
        candidates = [
            {
                "protein_id": "GenomeA/protein_002",
                "rf_distance": 0.0,
                "informative_splits": 0,
                "bitscore": 120.0,
            },
            {
                "protein_id": "GenomeA/protein_001",
                "rf_distance": 0.0,
                "informative_splits": 0,
                "bitscore": 180.0,
            },
        ]

        best = marker_selection.choose_best_candidate(candidates)

        self.assertEqual(best["protein_id"], "GenomeA/protein_001")

    def test_choose_best_candidate_prefers_previous_choice_on_exact_rf_tie(self):
        candidates = [
            {
                "protein_id": "GenomeA/protein_002",
                "rf_distance": 0.0,
                "informative_splits": 0,
                "bitscore": 200.0,
            },
            {
                "protein_id": "GenomeA/protein_001",
                "rf_distance": 0.0,
                "informative_splits": 0,
                "bitscore": 120.0,
            },
        ]

        best = marker_selection.choose_best_candidate(
            candidates,
            preferred_protein_id="GenomeA/protein_001",
        )

        self.assertEqual(best["protein_id"], "GenomeA/protein_001")

    def test_choose_best_candidate_prefers_higher_contig_support_before_bitscore(self):
        candidates = [
            {
                "protein_id": "GenomeA/contig1/protein_001",
                "rf_distance": 0.0,
                "informative_splits": 4,
                "contig_marker_support": 1,
                "bitscore": 200.0,
            },
            {
                "protein_id": "GenomeA/contig2/protein_001",
                "rf_distance": 0.0,
                "informative_splits": 4,
                "contig_marker_support": 3,
                "bitscore": 100.0,
            },
        ]

        best = marker_selection.choose_best_candidate(candidates)

        self.assertEqual(best["protein_id"], "GenomeA/contig2/protein_001")

    def test_resolve_marker_tree_uses_seeded_assignment_on_exact_rf_tie(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "marker.nwk"
            species_tree = tmp / "species.nwk"
            table_path = tmp / "table.csv"

            marker_tree.write_text("((A|p1,B|b1),A|p2);\n")
            species_tree.write_text("(A,B);\n")
            pd.DataFrame(
                [
                    {"savedname": "A/p1", "score_bits": 100.0},
                    {"savedname": "A/p2", "score_bits": 200.0},
                    {"savedname": "B/b1", "score_bits": 150.0},
                ]
            ).to_csv(table_path, index=False)

            _, records = marker_selection.resolve_marker_tree(
                marker_tree_path=str(marker_tree),
                species_tree_path=str(species_tree),
                table_path=str(table_path),
                marker_name="MarkerX",
                ls_refs=None,
                selection_mode="coordinate",
                max_rounds=5,
                lock_references=False,
                initial_kept={("MarkerX", "A"): "A/p1"},
            )

            kept = {
                (row["genome"], row["protein_id"])
                for row in records
                if row["status"] == "Kept"
            }

            self.assertEqual(kept, {("A", "A/p1")})

    def test_coordinate_mode_recovers_jointly_better_assignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "marker.nwk"
            species_tree = tmp / "species.nwk"
            table_path = tmp / "table.csv"

            marker_tree.write_text("(((A|p1,B|b1),(C|p2,D|d1)),(A|p2,C|p1));\n")
            species_tree.write_text("((A,B),(C,D));\n")
            pd.DataFrame(
                [
                    {"savedname": "A/p1", "score_bits": 100.0},
                    {"savedname": "A/p2", "score_bits": 200.0},
                    {"savedname": "B/b1", "score_bits": 150.0},
                    {"savedname": "C/p1", "score_bits": 200.0},
                    {"savedname": "C/p2", "score_bits": 100.0},
                    {"savedname": "D/d1", "score_bits": 150.0},
                ]
            ).to_csv(table_path, index=False)

            _, legacy_records = marker_selection.resolve_marker_tree(
                marker_tree_path=str(marker_tree),
                species_tree_path=str(species_tree),
                table_path=str(table_path),
                marker_name="MarkerX",
                ls_refs=None,
                selection_mode="legacy",
                max_rounds=5,
                lock_references=False,
            )
            _, coordinate_records = marker_selection.resolve_marker_tree(
                marker_tree_path=str(marker_tree),
                species_tree_path=str(species_tree),
                table_path=str(table_path),
                marker_name="MarkerX",
                ls_refs=None,
                selection_mode="coordinate",
                max_rounds=5,
                lock_references=False,
            )

            legacy_kept = {
                (row["genome"], row["protein_id"])
                for row in legacy_records
                if row["status"] == "Kept"
            }
            coordinate_kept = {
                (row["genome"], row["protein_id"])
                for row in coordinate_records
                if row["status"] == "Kept"
            }

            self.assertEqual(legacy_kept, {("A", "A/p1"), ("C", "C/p2")})
            self.assertEqual(coordinate_kept, {("A", "A/p1"), ("C", "C/p1")})

    def test_process_tree_worker_preserves_internal_branch_support(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "MarkerX.nwk"
            species_tree = tmp / "species.nwk"
            table_path = tmp / "table.csv"
            outdir = tmp / "out"
            (outdir / "removed").mkdir(parents=True)
            (outdir / "protTrees" / "no_duplicates" / "out").mkdir(parents=True)
            marker_tree.write_text(
                "((A|ctg|a:1,B|ctg|b:1):1,"
                "(C|ctg|c:1,D|ctg|d:1)1.0:1,"
                "(E|ctg|e:1,F|ctg|f:1)0.73:1);\n"
            )
            species_tree.write_text("((A,B),(C,D),(E,F));\n")
            pd.DataFrame(
                [
                    {"savedname": f"{genome}/ctg/{gene}", "score_bits": 100.0}
                    for genome, gene in (
                        ("A", "a"),
                        ("B", "b"),
                        ("C", "c"),
                        ("D", "d"),
                        ("E", "e"),
                        ("F", "f"),
                    )
                ]
            ).to_csv(table_path, index=False)

            marker_selection._process_tree_worker(
                (
                    str(marker_tree),
                    str(table_path),
                    str(species_tree),
                    str(outdir),
                    None,
                    "coordinate",
                    5,
                    False,
                    None,
                )
            )

            cached_path = outdir / "protTrees" / "no_duplicates" / "out" / "_no_dups_MarkerX_.nw"
            cached_labels = Tree(str(cached_path), format=1)
            internal_labels = {
                node.name
                for node in cached_labels.traverse()
                if not node.is_leaf()
            }
            self.assertTrue({"", "1.0", "0.73"}.issubset(internal_labels))

            cached_tree = Tree(str(cached_path))
            supports = {
                round(float(node.support), 2)
                for node in cached_tree.traverse()
                if not node.is_leaf()
            }
            self.assertTrue({0.73, 1.0}.issubset(supports))

    def test_choose_tree_by_rf_prefers_original_when_pruning_worsens_rf(self):
        species = Tree("((A,B),(C,D));")
        original = Tree("((A,B),(C,D));")
        candidate = Tree("((A,C),(B,D));")

        chosen = marker_selection.choose_tree_by_rf(
            species_tree=species,
            original_tree=original,
            candidate_tree=candidate,
        )

        self.assertEqual(chosen.write(format=9), original.write(format=9))

    def test_effective_singleton_mode_returns_requested_runtime_mode(self):
        for mode in ("delta_rf", "composite", "contig_consensus", "recipient_consensus", "neighbor_clade", "neighbor_ml"):
            self.assertEqual(
                marker_selection.effective_singleton_mode(
                    mode,
                    0.40,
                    duplicate_resolution_present=True,
                ),
                mode,
            )

    def test_every_singleton_mode_uses_global_rf_gate(self):
        modes = (
            "delta_rf",
            "topoknn",
            "hybrid",
            "composite",
            "contig_consensus",
            "recipient_consensus",
            "neighbor_clade",
            "neighbor_ml",
            "gcp",
        )
        for mode in modes:
            with self.subTest(mode=mode):
                self.assertTrue(marker_selection.singleton_mode_uses_global_rf_gate(mode))

    def test_choose_singleton_prune_prefers_highest_delta_rf(self):
        species = Tree("((A,B),(C,(D,E)));")
        working = Tree("((A,B),(C,(D,E)));")
        candidate_a = Tree("(B,(C,(D,E)));")
        candidate_c = Tree("((A,B),(D,E));")

        with patch.object(
            marker_selection,
            "_score_singleton_candidates",
            return_value=[
                {
                    "leaf_name": "A",
                    "delta_rf": 0.30,
                    "topoknn_score": 0.40,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "candidate_tree": candidate_a,
                },
                {
                    "leaf_name": "C",
                    "delta_rf": 0.10,
                    "topoknn_score": 1.10,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "candidate_tree": candidate_c,
                },
            ],
        ):
            chosen = marker_selection.choose_singleton_prune(
                species_tree=species,
                working_tree=working,
                mode="delta_rf",
                k=3,
            )

        self.assertIsNotNone(chosen)
        self.assertEqual(chosen["leaf_name"], "A")

    def test_choose_singleton_prune_composite_abstains_on_low_signal(self):
        species = Tree("((A,B),(C,(D,E)));")
        working = Tree("((A,B),(C,(D,E)));")
        candidate_a = Tree("(B,(C,(D,E)));")
        candidate_c = Tree("((A,B),(D,E));")

        with patch.object(
            marker_selection,
            "_score_singleton_candidates",
            return_value=[
                {
                    "leaf_name": "A",
                    "delta_rf": 0.02,
                    "topoknn_score": 0.10,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "recipient_consensus_score": 0.0,
                    "candidate_tree": candidate_a,
                },
                {
                    "leaf_name": "C",
                    "delta_rf": 0.01,
                    "topoknn_score": 0.08,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "recipient_consensus_score": 0.0,
                    "candidate_tree": candidate_c,
                },
            ],
        ):
            chosen = marker_selection.choose_singleton_prune(
                species_tree=species,
                working_tree=working,
                mode="composite",
                k=3,
            )

        self.assertIsNone(chosen)

    def test_choose_singleton_prune_reuses_candidates_without_mutating_input(self):
        species = Tree("((A,B),(C,(D,E)));")
        working = Tree("((A,B),(C,(D,E)));")
        candidate_a = Tree("(B,(C,(D,E)));")
        candidate_c = Tree("((A,B),(D,E));")
        scored_candidates = [
            {
                "leaf_name": "A",
                "delta_rf": 0.12,
                "topoknn_score": 1.8,
                "branch_outlier": 0.2,
                "bitscore_outlier": 1.4,
                "recipient_consensus_score": 3.8,
                "candidate_tree": candidate_a,
            },
            {
                "leaf_name": "C",
                "delta_rf": 0.1,
                "topoknn_score": 0.2,
                "branch_outlier": 0.0,
                "bitscore_outlier": 0.0,
                "recipient_consensus_score": 0.3,
                "candidate_tree": candidate_c,
            },
        ]
        original_keys = [set(candidate) for candidate in scored_candidates]

        with patch.object(marker_selection, "_score_singleton_candidates") as scorer:
            chosen = marker_selection.choose_singleton_prune(
                species_tree=species,
                working_tree=working,
                mode="recipient_consensus",
                k=3,
                scored_candidates=scored_candidates,
            )

        scorer.assert_not_called()
        self.assertIsNotNone(chosen)
        self.assertEqual(chosen["leaf_name"], "A")
        self.assertGreater(chosen["score"], 0.0)
        self.assertEqual([set(candidate) for candidate in scored_candidates], original_keys)
        self.assertNotIn("recipient_rank_score", scored_candidates[0])

    def test_gcp_fallback_uses_recipient_consensus_ranking_for_small_panels(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "leaf_name": "Genome1|contig1|gene1",
                "contig_id": "contig1",
                "delta_rf": 0.12,
                "topoknn_score": 1.8,
                "recipient_consensus_score": 3.8,
                "taxa_count": 2,
            },
            {
                "marker_name": "MarkerA",
                "genome": "Genome2",
                "leaf_name": "Genome2|contig1|gene1",
                "contig_id": "contig1",
                "delta_rf": 0.01,
                "topoknn_score": 0.1,
                "recipient_consensus_score": 0.1,
                "taxa_count": 2,
            },
        ]

        with self.assertLogs("sgtree", level="WARNING") as logs:
            classified = marker_selection._classify_gcp_proposals(
                proposals,
                contig_marker_context={},
            )

        by_leaf = {proposal["leaf_name"]: proposal for proposal in classified}
        self.assertIn("Using recipient_consensus classification.", "\n".join(logs.output))
        self.assertEqual(
            by_leaf["Genome1|contig1|gene1"]["singleton_class"],
            "contamination_candidate",
        )
        self.assertEqual(by_leaf["Genome2|contig1|gene1"]["singleton_class"], "clean")
        self.assertGreater(by_leaf["Genome1|contig1|gene1"]["score"], by_leaf["Genome2|contig1|gene1"]["score"])

    def test_gcp_uses_recipient_consensus_fallback_for_three_and_four_marker_panels(self):
        for marker_count in (3, 4):
            with self.subTest(marker_count=marker_count):
                with self.assertLogs("sgtree", level="WARNING") as logs:
                    classified = marker_selection._classify_gcp_proposals(
                        self._synthetic_gcp_panel(marker_count),
                        contig_marker_context={},
                    )

                hits = [
                    proposal
                    for proposal in classified
                    if proposal.get("singleton_class") == "contamination_candidate"
                ]
                self.assertIn("Using recipient_consensus classification.", "\n".join(logs.output))
                self.assertEqual([proposal["leaf_name"] for proposal in hits], ["Genome00|contig1|gene00"])
                self.assertNotIn("gcp_score", hits[0])
                self.assertIn("recipient_rank_score", hits[0])

    def test_gcp_scores_clear_outlier_for_five_and_ten_marker_panels(self):
        for marker_count in (5, 10):
            with self.subTest(marker_count=marker_count):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", SyntaxWarning)
                    with self.assertNoLogs("sgtree", level="WARNING"):
                        classified = marker_selection._classify_gcp_proposals(
                            self._synthetic_gcp_panel(marker_count),
                            contig_marker_context={},
                        )

                hits = [
                    proposal
                    for proposal in classified
                    if proposal.get("singleton_class") == "contamination_candidate"
                ]
                self.assertEqual([proposal["leaf_name"] for proposal in hits], ["Genome00|contig1|gene00"])
                self.assertIn("gcp_score", hits[0])
                self.assertGreaterEqual(hits[0]["gcp_score"], marker_selection.GCP_COMBINED_THRESHOLD)

    def test_choose_singleton_prune_neighbor_clade_accepts_local_outlier_without_positive_rf(self):
        species = Tree("((A,B),(C,(D,E)));")
        working = Tree("((A,B),(C,(D,E)));")
        candidate_a = Tree("(B,(C,(D,E)));")
        candidate_c = Tree("((A,B),(D,E));")

        with patch.object(
            marker_selection,
            "_score_singleton_candidates",
            return_value=[
                {
                    "leaf_name": "A",
                    "delta_rf": 0.0,
                    "topoknn_score": 0.2,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "recipient_consensus_score": 1.2,
                    "species_anchor_score": 0.9,
                    "neighbor_anchor_purity": 1.0,
                    "join_purity": 0.4,
                    "purity_drop": 0.6,
                    "anchor_knn_agreement": 0.0,
                    "attachment_gap": 1.5,
                    "present_neighbor_count": 4,
                    "neighbor_clade_score": 4.8,
                    "candidate_tree": candidate_a,
                },
                {
                    "leaf_name": "C",
                    "delta_rf": -0.02,
                    "topoknn_score": 0.1,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "recipient_consensus_score": 0.3,
                    "species_anchor_score": 0.7,
                    "neighbor_anchor_purity": 0.8,
                    "join_purity": 0.7,
                    "purity_drop": 0.1,
                    "anchor_knn_agreement": 0.75,
                    "attachment_gap": 0.1,
                    "present_neighbor_count": 4,
                    "neighbor_clade_score": 1.4,
                    "candidate_tree": candidate_c,
                },
            ],
        ):
            chosen = marker_selection.choose_singleton_prune(
                species_tree=species,
                working_tree=working,
                mode="neighbor_clade",
                k=4,
            )

        self.assertIsNotNone(chosen)
        self.assertEqual(chosen["leaf_name"], "A")

    def test_choose_singleton_prune_hybrid_requires_consistent_support(self):
        species = Tree("((A,B),(C,(D,E)));")
        working = Tree("((A,B),(C,(D,E)));")
        candidate_a = Tree("(B,(C,(D,E)));")
        candidate_c = Tree("((A,B),(D,E));")

        with patch.object(
            marker_selection,
            "_score_singleton_candidates",
            return_value=[
                {
                    "leaf_name": "A",
                    "delta_rf": 0.25,
                    "topoknn_score": 0.15,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "candidate_tree": candidate_a,
                },
                {
                    "leaf_name": "C",
                    "delta_rf": 0.02,
                    "topoknn_score": 1.40,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "candidate_tree": candidate_c,
                },
            ],
        ):
            chosen = marker_selection.choose_singleton_prune(
                species_tree=species,
                working_tree=working,
                mode="hybrid",
                k=3,
            )

        self.assertIsNone(chosen)

    def test_prune_singletons_removes_only_selected_leaf(self):
        species = Tree("((A,B),(C,(D,E)));")
        working = Tree("((A,B),(C,(D,E)));")

        with patch.object(
            marker_selection,
            "choose_singleton_prune",
            return_value={
                "leaf_name": "A",
                "candidate_tree": Tree("(B,(C,(D,E)));"),
            },
        ):
            pruned = marker_selection.prune_singletons(
                species_tree=species,
                working_tree=working,
                mode="delta_rf",
                k=3,
            )

        self.assertEqual(sorted(leaf.name for leaf in pruned.iter_leaves()), ["B", "C", "D", "E"])

    def test_choose_singleton_prune_outlier_picks_obvious_misplaced_leaf(self):
        species = Tree("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);")
        working = Tree("((((A:1,B:1):1,D:3):1,C:1):1,E:1);")

        chosen = marker_selection.choose_singleton_prune(
            species_tree=species,
            working_tree=working,
            mode="outlier",
            k=3,
        )

        self.assertIsNotNone(chosen)
        self.assertEqual(chosen["leaf_name"], "D")
        self.assertEqual(chosen["genome"], "D")
        self.assertGreater(chosen["score"], 0.0)

    def test_select_singleton_proposals_keeps_last_marker_for_genome(self):
        proposals = [
            {"marker_name": "MarkerA", "genome": "Genome1", "score": 5.0},
            {"marker_name": "MarkerB", "genome": "Genome1", "score": 4.0},
            {"marker_name": "MarkerC", "genome": "Genome2", "score": 3.0},
        ]

        accepted = marker_selection.select_singleton_proposals(
            proposals,
            genome_marker_counts={"Genome1": 2, "Genome2": 1},
            min_markers_per_genome=1,
        )

        self.assertEqual(
            {proposal["marker_name"] for proposal in accepted},
            {"MarkerA"},
        )

    def test_select_singleton_proposals_defaults_to_one_prune_per_genome(self):
        proposals = [
            {"marker_name": "MarkerA", "genome": "Genome1", "score": 5.0},
            {"marker_name": "MarkerB", "genome": "Genome1", "score": 4.0},
            {"marker_name": "MarkerC", "genome": "Genome2", "score": 3.0},
            {"marker_name": "MarkerD", "genome": "Genome2", "score": 2.0},
        ]

        accepted = marker_selection.select_singleton_proposals(
            proposals,
            genome_marker_counts={"Genome1": 10, "Genome2": 10},
            min_markers_per_genome=1,
        )

        self.assertEqual(
            {proposal["marker_name"] for proposal in accepted},
            {"MarkerA", "MarkerC"},
        )

    def test_singleton_proposals_from_results_neighbor_clade_keeps_all_candidates(self):
        results = [
            {
                "marker_name": "MarkerA",
                "chosen": {"leaf_name": "Genome1|contig1|gene1"},
                "candidates": [
                    {"leaf_name": "Genome1|contig1|gene1", "genome": "Genome1"},
                    {"leaf_name": "Genome2|contig1|gene1", "genome": "Genome2"},
                ],
            }
        ]

        proposals, proposal_keys = marker_selection.singleton_proposals_from_results(
            results,
            mode="neighbor_clade",
        )

        self.assertEqual(
            {(proposal["marker_name"], proposal["leaf_name"]) for proposal in proposals},
            {
                ("MarkerA", "Genome1|contig1|gene1"),
                ("MarkerA", "Genome2|contig1|gene1"),
            },
        )
        self.assertEqual(
            proposal_keys,
            {
                ("MarkerA", "Genome1|contig1|gene1"),
                ("MarkerA", "Genome2|contig1|gene1"),
            },
        )

    def test_singleton_proposals_from_results_legacy_keeps_only_chosen_candidate(self):
        results = [
            {
                "marker_name": "MarkerA",
                "chosen": {"leaf_name": "Genome1|contig1|gene1", "genome": "Genome1"},
                "candidates": [
                    {"leaf_name": "Genome1|contig1|gene1", "genome": "Genome1"},
                    {"leaf_name": "Genome2|contig1|gene1", "genome": "Genome2"},
                ],
            }
        ]

        proposals, proposal_keys = marker_selection.singleton_proposals_from_results(
            results,
            mode="recipient_consensus",
        )

        self.assertEqual(
            {(proposal["marker_name"], proposal["leaf_name"]) for proposal in proposals},
            {("MarkerA", "Genome1|contig1|gene1")},
        )
        self.assertEqual(
            proposal_keys,
            {("MarkerA", "Genome1|contig1|gene1")},
        )

    def test_singleton_proposals_from_results_legacy_with_refs_keeps_all_query_candidates(self):
        results = [
            {
                "marker_name": "MarkerA",
                "chosen": {"leaf_name": "Ref1|contig1|gene1", "genome": "Ref1"},
                "candidates": [
                    {
                        "leaf_name": "Ref1|contig1|gene1",
                        "genome": "Ref1",
                        "delta_rf": 0.9,
                        "topoknn_score": 2.0,
                        "recipient_consensus_score": 2.5,
                    },
                    {
                        "leaf_name": "Query1|contig1|gene1",
                        "genome": "Query1",
                        "delta_rf": 0.8,
                        "topoknn_score": 1.9,
                        "recipient_consensus_score": 2.1,
                    },
                    {
                        "leaf_name": "Query2|contig1|gene1",
                        "genome": "Query2",
                        "delta_rf": 0.4,
                        "topoknn_score": 1.6,
                        "recipient_consensus_score": 1.7,
                    },
                ],
            }
        ]

        proposals, proposal_keys = marker_selection.singleton_proposals_from_results(
            results,
            mode="recipient_consensus",
            reference_genomes={"Ref1"},
        )

        self.assertEqual(
            {(proposal["marker_name"], proposal["leaf_name"]) for proposal in proposals},
            {
                ("MarkerA", "Query1|contig1|gene1"),
                ("MarkerA", "Query2|contig1|gene1"),
            },
        )
        self.assertEqual(
            proposal_keys,
            {
                ("MarkerA", "Query1|contig1|gene1"),
                ("MarkerA", "Query2|contig1|gene1"),
            },
        )
        scores = {proposal["leaf_name"]: float(proposal["score"]) for proposal in proposals}
        self.assertGreater(scores["Query1|contig1|gene1"], scores["Query2|contig1|gene1"])

    def test_filter_reference_singleton_proposals_removes_reference_candidates(self):
        proposals = [
            {"marker_name": "MarkerA", "leaf_name": "Ref1|contig1|gene1", "genome": "Ref1"},
            {"marker_name": "MarkerA", "leaf_name": "Query1|contig1|gene1", "genome": "Query1"},
        ]
        proposal_keys = {
            ("MarkerA", "Ref1|contig1|gene1"),
            ("MarkerA", "Query1|contig1|gene1"),
        }

        filtered, filtered_keys = marker_selection._filter_reference_singleton_proposals(
            proposals,
            proposal_keys,
            reference_genomes={"Ref1"},
        )

        self.assertEqual(
            {(proposal["marker_name"], proposal["leaf_name"]) for proposal in filtered},
            {("MarkerA", "Query1|contig1|gene1")},
        )
        self.assertEqual(filtered_keys, {("MarkerA", "Query1|contig1|gene1")})

    def test_select_neighbor_ml_proposals_uses_genome_first_policy_for_small_panels(self):
        scored = [
            {
                "marker_name": f"Marker{i:02d}",
                "leaf_name": f"Genome{i:02d}|contig1|gene1",
                "genome": f"Genome{i:02d}",
                "taxa_count": 50,
                "high_confidence": True,
                "shape_penalized_score_v4": float(100 - i),
                "genome_first_score_v8": float(100 - i),
            }
            for i in range(40)
        ]

        with patch.object(marker_selection, "score_neighbor_ml_proposals", return_value=scored):
            proposals, accepted = marker_selection.select_neighbor_ml_proposals([{"marker_name": "unused"}])

        self.assertTrue(all(proposal["policy_variant"] == "genome_first_score_v8" for proposal in proposals))
        self.assertEqual(len(accepted), 40)
        self.assertEqual(
            {proposal["genome"] for proposal in accepted},
            {f"Genome{i:02d}" for i in range(40)},
        )

    def test_select_neighbor_ml_proposals_uses_shape_policy_for_large_panels(self):
        scored = [
            {
                "marker_name": f"Marker{i:03d}",
                "leaf_name": f"Genome{i:03d}|contig1|gene1",
                "genome": f"Genome{i:03d}",
                "taxa_count": 100,
                "high_confidence": True,
                "shape_penalized_score_v4": float(200 - i),
                "genome_first_score_v8": float(i),
            }
            for i in range(70)
        ]

        with patch.object(marker_selection, "score_neighbor_ml_proposals", return_value=scored):
            proposals, accepted = marker_selection.select_neighbor_ml_proposals([{"marker_name": "unused"}])

        self.assertTrue(all(proposal["policy_variant"] == "shape_penalized_score_v4" for proposal in proposals))
        self.assertEqual(len(accepted), 70)
        self.assertEqual(
            {proposal["genome"] for proposal in accepted},
            {f"Genome{i:03d}" for i in range(70)},
        )

    def test_prune_tree_to_query_genomes_removes_references(self):
        tree = Tree("((Ref1|r1,Query1|q1),(Ref2|r2,Query2|q2));")

        pruned = marker_selection._prune_tree_to_query_genomes(tree, {"Ref1", "Ref2"})

        self.assertEqual(
            sorted(leaf.name for leaf in pruned.iter_leaves()),
            ["Query1|q1", "Query2|q2"],
        )

    def test_propose_singleton_prune_worker_keeps_references_in_scoring_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "aln_MarkerA_out.nw"
            species_tree = tmp / "species.nwk"
            ref_dir = tmp / "refs"
            ref_dir.mkdir()
            (ref_dir / "Ref1.faa").write_text(">Ref1\nM\n")
            marker_tree.write_text("((Ref1|r1,Query1|q1),(Query2|q2,Query3|q3));\n")
            species_tree.write_text("((Ref1,Query1),(Query2,Query3));\n")

            seen = {}

            def _fake_score_singleton_candidates(species_tree, working_tree, **_kwargs):
                seen["species_leaves"] = sorted(leaf.name for leaf in species_tree.iter_leaves())
                seen["working_leaves"] = sorted(leaf.name for leaf in working_tree.iter_leaves())
                return []

            with patch.object(
                marker_selection,
                "_score_singleton_candidates",
                side_effect=_fake_score_singleton_candidates,
            ):
                result = marker_selection._propose_singleton_prune_worker(
                    (
                        str(marker_tree),
                        str(species_tree),
                        0,
                        0.25,
                        "neighbor_ml",
                        str(tmp / "missing_table.csv"),
                        set(),
                        str(tmp / "missing_alignment.faa"),
                        str(ref_dir),
                    )
                )

        self.assertEqual(
            seen["species_leaves"],
            ["Query1", "Query2", "Query3", "Ref1"],
        )
        self.assertEqual(
            seen["working_leaves"],
            ["Query1|q1", "Query2|q2", "Query3|q3", "Ref1|r1"],
        )
        self.assertEqual(result["candidates"], [])

    def test_propose_singleton_prune_worker_scores_candidates_once(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "_no_dups_MarkerA_.nw"
            species_tree = tmp / "species.nwk"
            marker_tree.write_text(
                "(((A|ctg|a,B|ctg|b),C|ctg|c),(D|ctg|d,E|ctg|e));\n"
            )
            species_tree.write_text("(((A,B),C),(D,E));\n")

            with patch.object(
                marker_selection,
                "_score_singleton_candidates",
                wraps=marker_selection._score_singleton_candidates,
            ) as score_candidates:
                marker_selection._propose_singleton_prune_worker(
                    (
                        str(marker_tree),
                        str(species_tree),
                        0,
                        0.0,
                        "delta_rf",
                        str(tmp / "missing_table.csv"),
                        set(),
                        str(tmp / "missing_alignment.faa"),
                        None,
                    )
                )

            self.assertEqual(score_candidates.call_count, 1)

    def test_propose_loo_profile_worker_bypasses_legacy_candidate_scoring(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "_no_dups_MarkerA_.nw"
            species_tree = tmp / "species.nwk"
            leaf_names = [f"G{index}|ctg{index}|gene{index}" for index in range(7)]
            marker_tree.write_text(
                "(((G0|ctg0|gene0,G1|ctg1|gene1),(G2|ctg2|gene2,G3|ctg3|gene3)),"
                "((G4|ctg4|gene4,G5|ctg5|gene5),G6|ctg6|gene6));\n"
            )
            species_tree.write_text("(((G0,G1),(G2,G3)),((G4,G5),G6));\n")

            with (
                patch.object(marker_selection, "_score_singleton_candidates") as scorer,
                patch.object(marker_selection, "choose_singleton_prune") as chooser,
            ):
                result = marker_selection._propose_singleton_prune_worker(
                    (
                        str(marker_tree),
                        str(species_tree),
                        0,
                        0.0,
                        "loo_profile",
                        str(tmp / "missing_table.csv"),
                        set(),
                        str(tmp / "missing_alignment.faa"),
                        None,
                    )
                )

            scorer.assert_not_called()
            chooser.assert_not_called()
            self.assertIsNone(result["chosen"])
            self.assertEqual(result["mode"], "loo_profile")
            self.assertEqual(
                result["candidates"],
                [
                    {
                        "leaf_name": leaf_name,
                        "taxa_count": 7,
                        "genome": f"G{index}",
                        "contig_id": f"ctg{index}",
                        "gene_id": f"gene{index}",
                    }
                    for index, leaf_name in enumerate(sorted(leaf_names))
                ],
            )

    def test_remove_singles_loo_profile_reports_without_pruning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            outdir = tmp / "run"
            cached_dir = outdir / "protTrees" / "no_duplicates" / "out"
            no_singles_dir = outdir / "protTrees" / "no_singles"
            cached_dir.mkdir(parents=True)
            no_singles_dir.mkdir(parents=True)
            cached_tree = cached_dir / "_no_dups_MarkerA_.nw"
            cached_tree.write_text(
                "(((G0|ctg0|gene0:1,G1|ctg1|gene1:1)0.95:1,"
                "(G2|ctg2|gene2:1,G3|ctg3|gene3:1)0.95:1)0.95:1,"
                "((G4|ctg4|gene4:1,G5|ctg5|gene5:1)0.95:1,"
                "(G6|ctg6|gene6:1,G7|ctg7|gene7:1)0.95:1)0.95:1);\n"
            )
            original_bytes = cached_tree.read_bytes()
            species_tree = outdir / "tree.nwk"
            species_tree.write_text("(((G0,G1),(G2,G3)),((G4,G5),(G6,G7)));\n")
            candidates = [
                {
                    "leaf_name": f"G{index}|ctg{index}|gene{index}",
                    "taxa_count": 8,
                    "genome": f"G{index}",
                    "contig_id": f"ctg{index}",
                    "gene_id": f"gene{index}",
                }
                for index in range(8)
            ]
            worker_result = {
                "filepath": str(cached_tree),
                "marker_name": "MarkerA",
                "mode": "loo_profile",
                "rdist": 0.0,
                "num_nei": 7,
                "chosen": None,
                "candidates": candidates,
            }
            cfg = SimpleNamespace(
                outdir=str(outdir),
                num_nei=0,
                singles_min_rfdist=0.25,
                singles_mode="loo_profile",
                num_cpus=1,
                ref=None,
            )
            tree_class = marker_selection.Tree
            score_impl = marker_selection.score_loo_profiles
            writer_impl = marker_selection._write_singleton_candidate_table

            with (
                patch.object(marker_selection, "map_processed", return_value=[worker_result]),
                patch.object(marker_selection, "Tree", wraps=tree_class) as tree_loader,
                patch.object(marker_selection, "score_loo_profiles", wraps=score_impl) as scorer,
                patch.object(
                    marker_selection,
                    "_write_singleton_candidate_table",
                    wraps=writer_impl,
                ) as table_writer,
                patch.object(
                    marker_selection,
                    "singleton_proposals_from_results",
                    side_effect=AssertionError("legacy proposal path called"),
                ),
                patch.object(
                    marker_selection,
                    "classify_singleton_proposals",
                    side_effect=AssertionError("legacy classification path called"),
                ),
                patch.object(
                    marker_selection,
                    "select_singleton_proposals",
                    side_effect=AssertionError("legacy selection path called"),
                ),
                patch.object(
                    marker_selection,
                    "build_singleton_output_tree",
                    side_effect=AssertionError("RF/pruning path called"),
                ),
            ):
                marker_selection.remove_singles(cfg, species_tree_path=str(species_tree))

            tree_loader.assert_called_once_with(str(cached_tree), format=1)
            scorer.assert_called_once()
            self.assertEqual(scorer.call_args.kwargs["excluded_target_genomes"], set())
            self.assertEqual(table_writer.call_args.kwargs["accepted_keys"], set())
            self.assertEqual(table_writer.call_args.kwargs["proposal_keys"], set())

            copied_tree = no_singles_dir / cached_tree.name
            self.assertEqual(copied_tree.read_bytes(), original_bytes)
            self.assertEqual(
                sorted(leaf.name for leaf in tree_class(str(copied_tree), format=1).iter_leaves()),
                sorted(leaf.name for leaf in tree_class(str(cached_tree), format=1).iter_leaves()),
            )

            table = pd.read_csv(
                outdir / "singleton_candidates.tsv",
                sep="\t",
                dtype=str,
                keep_default_na=False,
            )
            self.assertEqual(len(table), 8)
            self.assertEqual(set(table["decision"]), {"kept_report_only"})
            self.assertEqual(set(table["loo_decision"]), {"kept_report_only"})
            self.assertEqual(set(table["loo_class"]), {"ambiguous"})
            self.assertEqual(set(table["loo_abstention_reason"]), {"insufficient_voters"})
            self.assertEqual(set(table["loo_voter_count"]), {"0"})
            for column in (
                "delta_rf",
                "loo_target_discordance",
                "loo_voter_center",
                "loo_voter_mad",
                "loo_voter_upper",
                "loo_conflict_margin",
                "loo_score",
                "loo_marker_rank",
                "loo_marker_margin",
                "loo_conflict_beyond_dispersion",
            ):
                self.assertEqual(set(table[column]), {""})

    def test_classify_singleton_proposals_marks_hgt_when_contig_has_other_clean_markers(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "contig_id": "contig1",
                "leaf_name": "Genome1|contig1|gene1",
            }
        ]

        classified = marker_selection.classify_singleton_proposals(
            proposals,
            contig_marker_context={("Genome1", "contig1"): {"MarkerA", "MarkerB"}},
        )

        self.assertEqual(classified[0]["singleton_class"], "hgt_candidate")

    def test_classify_singleton_proposals_marks_contamination_when_all_markers_on_contig_are_suspect(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "contig_id": "contig1",
                "leaf_name": "Genome1|contig1|gene1",
            },
            {
                "marker_name": "MarkerB",
                "genome": "Genome1",
                "contig_id": "contig1",
                "leaf_name": "Genome1|contig1|gene2",
            },
        ]

        classified = marker_selection.classify_singleton_proposals(
            proposals,
            contig_marker_context={("Genome1", "contig1"): {"MarkerA", "MarkerB"}},
        )

        self.assertEqual(
            {proposal["singleton_class"] for proposal in classified},
            {"contamination_candidate"},
        )

    def test_classify_singleton_proposals_contig_consensus_rescues_native_contig_contaminant(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "contig_id": "contig1",
                "leaf_name": "Genome1|contig1|gene1",
                "recipient_neighbor_overlap": 0.0,
                "recipient_consensus_score": 4.0,
            }
        ]

        classified = marker_selection.classify_singleton_proposals(
            proposals,
            contig_marker_context={("Genome1", "contig1"): {"MarkerA", "MarkerB", "MarkerC"}},
            marker_neighbor_context={
                ("Genome1", "contig1", "MarkerB"): {"recipient_neighbor_overlap": 0.9},
                ("Genome1", "contig1", "MarkerC"): {"recipient_neighbor_overlap": 0.8},
            },
            mode="contig_consensus",
        )

        self.assertEqual(classified[0]["singleton_class"], "contamination_candidate")
        self.assertGreaterEqual(classified[0]["contig_consensus_score"], 0.8)

    def test_classify_singleton_proposals_contig_consensus_marks_fully_discordant_contig_ambiguous(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "contig_id": "contig1",
                "leaf_name": "Genome1|contig1|gene1",
                "recipient_neighbor_overlap": 0.1,
                "recipient_consensus_score": 1.0,
            }
        ]

        classified = marker_selection.classify_singleton_proposals(
            proposals,
            contig_marker_context={("Genome1", "contig1"): {"MarkerA", "MarkerB", "MarkerC"}},
            marker_neighbor_context={
                ("Genome1", "contig1", "MarkerB"): {"recipient_neighbor_overlap": 0.2},
                ("Genome1", "contig1", "MarkerC"): {"recipient_neighbor_overlap": 0.3},
            },
            mode="contig_consensus",
        )

        self.assertEqual(classified[0]["singleton_class"], "ambiguous")

    def test_build_singleton_output_tree_neighbor_clade_applies_rf_guard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "marker.nwk"
            species_tree = tmp / "species.nwk"
            marker_tree.write_text("((A|a1,B|b1),(C|c1,D|d1));\n")
            species_tree.write_text("((A,B),(C,D));\n")

            chosen_tree, decision = marker_selection.build_singleton_output_tree(
                marker_tree_path=str(marker_tree),
                species_tree_path=str(species_tree),
                accepted_leaf_names=["A|a1", "C|c1"],
            )

            self.assertEqual(
                sorted(leaf.name for leaf in chosen_tree.iter_leaves()),
                ["A|a1", "B|b1", "C|c1", "D|d1"],
            )
            self.assertEqual(decision, "kept_rf_guard")

    def test_build_singleton_output_tree_legacy_prunes_multiple_accepted_leaves(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            marker_tree = tmp / "marker.nwk"
            species_tree = tmp / "species.nwk"
            marker_tree.write_text("(((A|a1,B|b1),E|e1),(C|c1,D|d1));\n")
            species_tree.write_text("((((B,E),D),A),C);\n")

            chosen_tree, decision = marker_selection.build_singleton_output_tree(
                marker_tree_path=str(marker_tree),
                species_tree_path=str(species_tree),
                accepted_leaf_names=["A|a1", "C|c1"],
            )

            self.assertEqual(
                sorted(leaf.name for leaf in chosen_tree.iter_leaves()),
                ["B|b1", "D|d1", "E|e1"],
            )
            self.assertEqual(decision, "pruned")

    def test_classify_singleton_proposals_neighbor_clade_marks_high_confidence_candidate(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "contig_id": "unknown_contig",
                "leaf_name": "Genome1|unknown_contig|gene1",
                "species_anchor_score": 0.9,
                "present_neighbor_count": 4,
                "target_neighbor_count": 5,
                "neighbor_anchor_purity": 1.0,
                "join_purity": 0.4,
                "purity_drop": 0.6,
                "anchor_knn_agreement": 0.0,
                "recipient_consensus_score": 1.4,
                "neighbor_clade_score": 4.2,
            }
        ]

        classified = marker_selection.classify_singleton_proposals(
            proposals,
            contig_marker_context={},
            mode="neighbor_clade",
        )

        self.assertEqual(classified[0]["singleton_class"], "contamination_candidate")

    def test_classify_singleton_proposals_neighbor_clade_marks_weak_anchor_ambiguous(self):
        proposals = [
            {
                "marker_name": "MarkerA",
                "genome": "Genome1",
                "contig_id": "unknown_contig",
                "leaf_name": "Genome1|unknown_contig|gene1",
                "species_anchor_score": 0.35,
                "present_neighbor_count": 2,
                "target_neighbor_count": 5,
                "neighbor_anchor_purity": 1.0,
                "join_purity": 0.5,
                "purity_drop": 0.5,
                "anchor_knn_agreement": 0.0,
                "recipient_consensus_score": 1.8,
                "neighbor_clade_score": 3.5,
            }
        ]

        classified = marker_selection.classify_singleton_proposals(
            proposals,
            contig_marker_context={},
            mode="neighbor_clade",
        )

        self.assertEqual(classified[0]["singleton_class"], "ambiguous")


class CachedLoaderTests(unittest.TestCase):
    """Cover the process-local caches used by RF-distance workers."""

    def setUp(self):
        marker_selection._SCORE_TABLE_CACHE.clear()
        marker_selection._SPECIES_TREE_CACHE.clear()

    def tearDown(self):
        marker_selection._SCORE_TABLE_CACHE.clear()
        marker_selection._SPECIES_TREE_CACHE.clear()

    def test_cached_score_table_parses_once_and_returns_same_object(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = str(Path(tmpdir) / "table.csv")
            pd.DataFrame(
                [
                    {"savedname": "A/p1", "score_bits": 100.0},
                    {"savedname": "B/p1", "score_bits": 150.0},
                ]
            ).to_csv(path, index=False)

            with patch.object(
                marker_selection,
                "_load_score_table",
                wraps=marker_selection._load_score_table,
            ) as wrapped_load:
                first_df, first_col = marker_selection._cached_score_table(path)
                second_df, second_col = marker_selection._cached_score_table(path)

            self.assertEqual(wrapped_load.call_count, 1)
            self.assertIs(first_df, second_df)
            self.assertEqual(first_col, second_col)

    def test_cached_species_tree_parses_once_and_returns_fresh_copies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = str(Path(tmpdir) / "species.nwk")
            Path(path).write_text("((A,B),(C,D));\n")

            self.assertNotIn(path, marker_selection._SPECIES_TREE_CACHE)
            first = marker_selection._cached_species_tree(path)
            self.assertIn(path, marker_selection._SPECIES_TREE_CACHE)
            cached_ref = marker_selection._SPECIES_TREE_CACHE[path]

            second = marker_selection._cached_species_tree(path)
            # Same underlying cached object, but distinct copies returned.
            self.assertIs(marker_selection._SPECIES_TREE_CACHE[path], cached_ref)
            self.assertIsNot(first, second)
            self.assertIsNot(first, cached_ref)
            self.assertIsNot(second, cached_ref)

            # Mutating one copy does not affect the other or the cache.
            first.prune(["A", "B"])
            self.assertEqual(
                sorted(leaf.name for leaf in second.iter_leaves()),
                ["A", "B", "C", "D"],
            )
            self.assertEqual(
                sorted(leaf.name for leaf in cached_ref.iter_leaves()),
                ["A", "B", "C", "D"],
            )


class BuildMarkerNeighborContextTests(unittest.TestCase):
    """Cover serial- vs threaded-execution equivalence for neighbor context."""

    def setUp(self):
        marker_selection._SPECIES_TREE_CACHE.clear()

    def tearDown(self):
        marker_selection._SPECIES_TREE_CACHE.clear()

    def _build_fixture(self, tmpdir: Path) -> dict:
        species_tree = tmpdir / "species.nwk"
        species_tree.write_text("(((A,B),C),(D,E));\n")

        # Two marker trees: MarkerA (congruent with species), MarkerB (shuffled).
        marker_a = tmpdir / "aln_MarkerA_out.nw"
        marker_b = tmpdir / "aln_MarkerB_out.nw"
        marker_a.write_text("(((A|contig1|g1,B|contig1|g1),C|contig1|g1),(D|contig1|g1,E|contig1|g1));\n")
        marker_b.write_text("(((A|contig1|g2,D|contig1|g2),C|contig1|g2),(B|contig1|g2,E|contig1|g2));\n")

        # Minimal table with the contig-hit rows for both markers. The
        # marker id is the last '/'-separated segment of ``namemodel`` per
        # _load_contig_marker_hit_map.
        table = tmpdir / "table_elim_dups"
        rows = []
        for marker, suffix in [("MarkerA", "g1"), ("MarkerB", "g2")]:
            for genome in ["A", "B", "C", "D", "E"]:
                rows.append(
                    {
                        "savedname": f"{genome}/contig1/{suffix}",
                        "score_bits": 100.0,
                        "namemodel": f"pfam/{marker}",
                    }
                )
        pd.DataFrame(rows).to_csv(table, index=False)

        return {
            "files": [str(marker_a), str(marker_b)],
            "species_tree_path": str(species_tree),
            "table_path": str(table),
        }

    def test_parallel_output_matches_serial(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = self._build_fixture(Path(tmpdir))

            serial = marker_selection._build_marker_neighbor_context(
                fixture["files"],
                species_tree_path=fixture["species_tree_path"],
                table_path=fixture["table_path"],
                k=3,
                num_cpus=1,
            )
            parallel = marker_selection._build_marker_neighbor_context(
                fixture["files"],
                species_tree_path=fixture["species_tree_path"],
                table_path=fixture["table_path"],
                k=3,
                num_cpus=4,
            )

            self.assertEqual(set(serial.keys()), set(parallel.keys()))
            self.assertGreater(len(serial), 0)
            for key in serial:
                self.assertAlmostEqual(
                    serial[key]["recipient_neighbor_overlap"],
                    parallel[key]["recipient_neighbor_overlap"],
                    places=12,
                )

    def test_empty_hit_map_returns_empty_dict(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fixture = self._build_fixture(Path(tmpdir))
            # Drop the namemodel column so _load_contig_marker_hit_map
            # returns an empty dict, short-circuiting the builder.
            table = Path(fixture["table_path"])
            pd.DataFrame([{"savedname": "A/contig1/g1", "score_bits": 100.0}]).to_csv(table, index=False)
            result = marker_selection._build_marker_neighbor_context(
                fixture["files"],
                species_tree_path=fixture["species_tree_path"],
                table_path=str(table),
                k=3,
                num_cpus=4,
            )
            self.assertEqual(result, {})


def _synthetic_ml_proposals():
    proposals = []
    for g_idx in range(6):
        for m_idx in range(3):
            is_outlier_genome = g_idx == 0
            is_outlier_marker = g_idx == 1 and m_idx == 0
            base = 2.0 if (is_outlier_genome or is_outlier_marker) else 0.3
            proposals.append(
                {
                    "genome": f"G{g_idx:02d}",
                    "marker_name": f"Marker{m_idx:02d}",
                    "leaf_name": f"G{g_idx:02d}|c1|gene{m_idx}",
                    "delta_rf": base * 0.1,
                    "neighbor_overlap": 0.5,
                    "topoknn_score": base * 0.5,
                    "branch_outlier": 0.0,
                    "bitscore_outlier": 0.0,
                    "recipient_consensus_score": base,
                    "species_anchor_score": 0.8,
                    "species_anchor_purity": 0.9,
                    "species_anchor_compactness_score": 0.5,
                    "species_long_branch_z": 0.0,
                    "species_long_branch_support": 0.0,
                    "present_neighbor_count": 5,
                    "present_neighbor_fraction": 0.7,
                    "neighbor_anchor_purity": 1.0,
                    "join_purity": 0.3 if is_outlier_genome else 0.8,
                    "purity_drop": base * 0.4,
                    "anchor_knn_agreement": 0.1 if is_outlier_genome else 0.7,
                    "attachment_gap": base,
                    "taxa_count": 6,
                }
            )
    return proposals


_ML_GOLDEN_ROWS = {
    ("Marker00", "G00|c1|gene0"): {
        "iforest_anomaly": 0.6907509445465396,
        "mahalanobis": 0.0,
        "ensemble_score_v1": 0.7152777777777778,
        "shape_penalized_score_v4": 0.3088352623456791,
        "genome_first_score_v8": 1.0833333333333335,
    },
    ("Marker00", "G01|c1|gene0"): {
        "iforest_anomaly": 0.7651128938196522,
        "mahalanobis": 0.0,
        "ensemble_score_v1": 0.8402777777777778,
        "shape_penalized_score_v4": 0.405994212962963,
        "genome_first_score_v8": 1.2083333333333333,
    },
    ("Marker00", "G02|c1|gene0"): {
        "iforest_anomaly": 0.42462132004490966,
        "mahalanobis": 0.0,
        "ensemble_score_v1": 0.4652777777777778,
        "shape_penalized_score_v4": 0.17977314814814815,
        "genome_first_score_v8": 0.5208333333333334,
    },
    ("Marker01", "G00|c1|gene1"): {
        "iforest_anomaly": 0.6168607356912175,
        "mahalanobis": 0.0,
        "ensemble_score_v1": 0.7152777777777778,
        "shape_penalized_score_v4": 0.2983780864197531,
        "genome_first_score_v8": 0.0,
    },
    ("Marker01", "G02|c1|gene1"): {
        "iforest_anomaly": 0.42462132004490966,
        "mahalanobis": 0.0,
        "ensemble_score_v1": 0.4652777777777778,
        "shape_penalized_score_v4": 0.17977314814814815,
        "genome_first_score_v8": 0.0,
    },
    ("Marker02", "G00|c1|gene2"): {
        "iforest_anomaly": 0.6168607356912175,
        "mahalanobis": 0.0,
        "ensemble_score_v1": 0.7152777777777778,
        "shape_penalized_score_v4": 0.2983780864197531,
        "genome_first_score_v8": 0.0,
    },
}


class ScoreNeighborMlProposalsRegression(unittest.TestCase):
    """Pin score_neighbor_ml_proposals outputs on a fixed synthetic panel.

    Guards against behavior drift when score_neighbor_ml_proposals is
    decomposed into helpers (Task 4.3 of the code-review remediation plan).
    Uses a 6-genome x 3-marker fixture with G00 as a global outlier and
    G01 as a single-marker outlier on Marker00. IsolationForest's
    random_state=42 and MinCovDet's support_fraction=0.8 make the output
    deterministic. MinCovDet is expected to fail on this rank-deficient
    fixture, exercising the mahalanobis=0.0 fallback.
    """

    def _score_synthetic_ml_proposals(self):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The covariance matrix associated to your dataset is not full rank",
                category=UserWarning,
            )
            with self.assertLogs("sgtree", level="WARNING") as logs:
                scored = marker_selection.score_neighbor_ml_proposals(
                    _synthetic_ml_proposals()
                )
        self.assertTrue(any("MinCovDet fit failed" in msg for msg in logs.output))
        return scored

    def test_pinned_columns_match_golden_snapshot(self):
        scored = self._score_synthetic_ml_proposals()
        self.assertEqual(len(scored), 18)

        df = (
            pd.DataFrame(scored)
            .sort_values(["marker_name", "leaf_name"])
            .reset_index(drop=True)
        )
        self.assertTrue(df["high_confidence"].all())

        for (marker, leaf), expected in _ML_GOLDEN_ROWS.items():
            mask = (df["marker_name"] == marker) & (df["leaf_name"] == leaf)
            self.assertEqual(int(mask.sum()), 1, f"row for {marker}/{leaf}")
            row = df.loc[mask].iloc[0]
            for col, value in expected.items():
                self.assertAlmostEqual(
                    float(row[col]),
                    value,
                    places=10,
                    msg=f"{col} mismatch for {marker}/{leaf}",
                )

    def test_genome_first_score_identifies_global_outlier(self):
        scored = self._score_synthetic_ml_proposals()
        df = pd.DataFrame(scored)
        # genome_first_score_v8 is assigned per genome on the top row; look
        # at the maximum value across all rows of a genome.
        by_genome = df.groupby("genome")["genome_first_score_v8"].max()
        self.assertEqual(by_genome.idxmax(), "G01")  # strongest single-marker signal
        self.assertGreater(by_genome["G00"], by_genome["G02"])

    def test_returns_zeros_when_no_high_confidence_rows(self):
        proposals = [
            {
                "genome": "G01",
                "marker_name": "MarkerA",
                "leaf_name": "G01|c1|gene1",
                "delta_rf": 0.1,
                "neighbor_overlap": 0.5,
                "topoknn_score": 0.1,
                "branch_outlier": 0.0,
                "bitscore_outlier": 0.0,
                "recipient_consensus_score": 0.3,
                "species_anchor_score": 0.8,
                "species_anchor_purity": 0.9,
                "species_anchor_compactness_score": 0.5,
                "species_long_branch_z": 0.0,
                "species_long_branch_support": 0.0,
                # Drop below NEIGHBOR_ML_HIGH_CONF_MIN_PRESENT=3 to force empty train.
                "present_neighbor_count": 0,
                "present_neighbor_fraction": 0.0,
                "neighbor_anchor_purity": 1.0,
                "join_purity": 0.8,
                "purity_drop": 0.2,
                "anchor_knn_agreement": 0.7,
                "attachment_gap": 0.3,
                "taxa_count": 3,
            }
        ]
        scored = marker_selection.score_neighbor_ml_proposals(proposals)
        self.assertEqual(len(scored), 1)
        for col in (
            "iforest_anomaly",
            "mahalanobis",
            "ensemble_score_v1",
            "shape_penalized_score_v4",
            "genome_first_score_v8",
        ):
            self.assertEqual(scored[0][col], 0.0, f"{col} should be zero on empty-train path")


if __name__ == "__main__":
    unittest.main()
