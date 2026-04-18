import unittest

import pandas as pd

from sgtree.search import _count_proteins_per_marker, _drop_rows_for_genomes


class CountProteinsPerMarkerTests(unittest.TestCase):
    def test_dedups_multi_domain_hits_on_same_protein_and_model(self):
        # GenomeA|p1 against M1 is hit twice (two domains) — must count as 1.
        # GenomeA also has p2/M1 (new protein, same model) and p3/M2.
        # GenomeB has p1/M1 and p2/M2.
        df = pd.DataFrame(
            {
                0: [
                    "GenomeA|p1",
                    "GenomeA|p1",
                    "GenomeA|p2",
                    "GenomeA|p3",
                    "GenomeB|p1",
                    "GenomeB|p2",
                ],
                3: ["M1", "M1", "M1", "M2", "M1", "M2"],
            }
        )
        result = _count_proteins_per_marker(df)
        self.assertEqual(result["GenomeA"], {"M1": 2, "M2": 1})
        self.assertEqual(result["GenomeB"], {"M1": 1, "M2": 1})

    def test_empty_input_yields_empty_dict(self):
        df = pd.DataFrame({0: pd.Series(dtype=str), 3: pd.Series(dtype=str)})
        self.assertEqual(_count_proteins_per_marker(df), {})

    def test_model_with_only_one_genome(self):
        df = pd.DataFrame({0: ["A|p1"], 3: ["M1"]})
        self.assertEqual(_count_proteins_per_marker(df), {"A": {"M1": 1}})

    def test_hit_on_ten_complete_plus_one_incomplete_genome(self):
        # 10 complete genomes (2 models each) + 1 incomplete (1 model only).
        # dict_counts len() by model tells incomplete detection apart.
        rows_0 = []
        rows_3 = []
        for g in range(10):
            rows_0.extend([f"G{g}|p1", f"G{g}|p2"])
            rows_3.extend(["M1", "M2"])
        # Incomplete genome hit only M1.
        rows_0.append("G10|p1")
        rows_3.append("M1")
        df = pd.DataFrame({0: rows_0, 3: rows_3})

        result = _count_proteins_per_marker(df)
        for g in range(10):
            self.assertEqual(set(result[f"G{g}"].keys()), {"M1", "M2"})
        self.assertEqual(set(result["G10"].keys()), {"M1"})


class DropRowsForGenomesTests(unittest.TestCase):
    def test_drops_only_matching_genomes_and_reports_count(self):
        df = pd.DataFrame(
            {
                0: ["A|p1", "B|p1", "A|p2", "C|p1"],
                3: ["M1", "M1", "M2", "M1"],
            }
        )
        filtered, n = _drop_rows_for_genomes(df, {"A", "C"})
        self.assertEqual(n, 3)
        self.assertEqual(list(filtered[0]), ["B|p1"])

    def test_no_removed_genomes_passes_through(self):
        df = pd.DataFrame({0: ["A|p1"], 3: ["M1"]})
        filtered, n = _drop_rows_for_genomes(df, set())
        self.assertEqual(n, 0)
        self.assertEqual(len(filtered), 1)

    def test_empty_dataframe(self):
        df = pd.DataFrame({0: pd.Series(dtype=str), 3: pd.Series(dtype=str)})
        filtered, n = _drop_rows_for_genomes(df, {"A"})
        self.assertEqual(n, 0)
        self.assertEqual(len(filtered), 0)


if __name__ == "__main__":
    unittest.main()
