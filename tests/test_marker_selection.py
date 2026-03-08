import tempfile
import unittest
from pathlib import Path

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


if __name__ == "__main__":
    unittest.main()
