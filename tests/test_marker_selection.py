import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from ete3 import Tree
from sgtree import marker_selection


class MarkerSelectionTests(unittest.TestCase):
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

    def test_effective_singleton_mode_pins_runtime_cleanup_to_delta_rf(self):
        self.assertEqual(
            marker_selection.effective_singleton_mode(
                "neighbor",
                0.40,
                duplicate_resolution_present=False,
            ),
            "delta_rf",
        )
        self.assertEqual(
            marker_selection.effective_singleton_mode(
                "neighbor",
                0.20,
                duplicate_resolution_present=False,
            ),
            "delta_rf",
        )
        self.assertEqual(
            marker_selection.effective_singleton_mode(
                "delta_rf",
                0.40,
                duplicate_resolution_present=False,
            ),
            "delta_rf",
        )
        self.assertEqual(
            marker_selection.effective_singleton_mode(
                "outlier",
                0.40,
                duplicate_resolution_present=True,
            ),
            "delta_rf",
        )

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


if __name__ == "__main__":
    unittest.main()
