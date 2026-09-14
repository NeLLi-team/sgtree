import tempfile
import unittest
from pathlib import Path

import pandas as pd

from sgtree.config import Config
from sgtree.search import (
    _count_proteins_per_marker,
    _drop_rows_for_genomes,
    build_working_df,
    parse_hmmsearch,
)


def _make_cfg(tmp: Path, **overrides) -> Config:
    kwargs = {
        "genomedir": str(tmp / "input"),
        "modeldir": str(tmp / "models.hmm"),
        "outdir": str(tmp / "run"),
        "num_cpus": 1,
        "percent_models": 0,
        "input_format": "faa",
        "lflt_fraction": 0.0,
        "aln_method": "hmmalign",
        "tree_method": "fasttree",
        "iqtree_fast": True,
        "iqtree_model": "LG+F+I+G4",
        "hmmsearch_cutoff": "cut_ga",
        "hmmsearch_evalue": 1e-5,
        "selection_mode": "coordinate",
        "selection_max_rounds": 5,
        "selection_global_rounds": 1,
        "lock_references": False,
        "max_sdup": -1,
        "max_dupl": -1.0,
        "ref": None,
        "ref_concat": str(tmp / "ref_cache"),
        "marker_selection": False,
        "singles": False,
        "singles_mode": "delta_rf",
        "num_nei": 0,
        "singles_min_rfdist": 0.25,
        "keep_intermediates": True,
        "is_ref": False,
        "start_time": "now",
    }
    kwargs.update(overrides)
    cfg = Config(**kwargs)
    Path(cfg.tables_dir).mkdir(parents=True, exist_ok=True)
    return cfg


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

    def test_underscores_do_not_merge_distinct_protein_model_pairs(self):
        df = pd.DataFrame(
            {
                0: ["GenomeA|c1|a_b", "GenomeA|c1|a"],
                3: ["C", "b_C"],
            }
        )

        result = _count_proteins_per_marker(df)

        self.assertEqual(result, {"GenomeA": {"C": 1, "b_C": 1}})

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


class BestScoringMarkerTests(unittest.TestCase):
    """A protein that hits several markers must land on its best-scoring one."""

    def test_best_scoring_marker_wins_over_hmm_file_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = _make_cfg(Path(tmpdir))
            # MarkerFirst comes first in HMM file order but scores 40 bits;
            # MarkerBest scores 300 and must win.
            finaldf = pd.DataFrame(
                [
                    {0: "GenomeA|c1|g1", 3: "MarkerFirst", 7: "40.0"},
                    {0: "GenomeB|c1|g9", 3: "MarkerFirst", 7: "90.0"},
                    {0: "GenomeA|c1|g1", 3: "MarkerBest", 7: "300.0"},
                ]
            )

            df, _ = build_working_df(cfg, finaldf)

            self.assertEqual(len(df), 2)
            kept = df.set_index("savedname")
            self.assertEqual(kept.loc["GenomeA/c1/g1", 3], "MarkerBest")
            self.assertEqual(kept.loc["GenomeA/c1/g1", "score_bits"], 300.0)
            self.assertEqual(kept.loc["GenomeB/c1/g9", 3], "MarkerFirst")

    def test_equal_scores_resolve_on_marker_name_whatever_the_row_order(self):
        rows = [
            {0: "GenomeA|c1|g1", 3: "MarkerZ", 7: "100.0"},
            {0: "GenomeA|c1|g1", 3: "MarkerA", 7: "100.0"},
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = _make_cfg(Path(tmpdir))
            for order in (rows, list(reversed(rows))):
                df, _ = build_working_df(cfg, pd.DataFrame(order))
                self.assertEqual(len(df), 1)
                self.assertEqual(df.iloc[0][3], "MarkerA")

    def test_multi_domain_rows_of_one_marker_keep_the_first_row(self):
        # Domains of the same (protein, marker) share the full-sequence score in
        # column 7, so the stable sort must leave the first domain row on top.
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = _make_cfg(Path(tmpdir))
            finaldf = pd.DataFrame(
                [
                    {0: "GenomeA|c1|g1", 3: "MarkerX", 7: "250.0", 13: "80.0"},
                    {0: "GenomeA|c1|g1", 3: "MarkerX", 7: "250.0", 13: "120.0"},
                ]
            )

            df, _ = build_working_df(cfg, finaldf)

            self.assertEqual(list(df[13]), ["80.0"])


HITS_TABLE = (
    "GenomeA|c1|g1 - 100 M1 - 100 1e-50 250.0\n"
    "GenomeA|c1|g2 - 300 M1 - 300 1e-50 300.0\n"
    "GenomeA|c1|g3 - 320 M1 - 320 1e-50 310.0\n"
)


class LengthFilterTests(unittest.TestCase):
    """The optional length filter removes exact target IDs."""

    def _cfg_with_hits(self, tmpdir: str):
        cfg = _make_cfg(Path(tmpdir), lflt_fraction=0.5)
        Path(cfg.hitsoutdir).write_text(HITS_TABLE)
        return cfg

    def test_filters_short_hits(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = self._cfg_with_hits(tmpdir)

            finaldf, _ = parse_hmmsearch(cfg)

            # g1 is 100 aa against a 300 aa median, below 0.5 * median.
            self.assertEqual(list(finaldf[0]), ["GenomeA|c1|g2", "GenomeA|c1|g3"])
            self.assertEqual(len(Path(cfg.hitsoutdir).read_text().splitlines()), 2)

    def test_keeps_id_with_removed_id_as_hyphenated_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = self._cfg_with_hits(tmpdir)
            Path(cfg.hitsoutdir).write_text(
                HITS_TABLE + "GenomeA|c1|g1-extra - 300 M1 - 300 1e-50 290.0\n"
            )

            finaldf, _ = parse_hmmsearch(cfg)

            self.assertEqual(
                list(finaldf[0]),
                ["GenomeA|c1|g2", "GenomeA|c1|g3", "GenomeA|c1|g1-extra"],
            )


class ParseHmmsearchAssignmentTests(unittest.TestCase):
    def test_multimarker_protein_does_not_make_incomplete_genome_pass(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = _make_cfg(Path(tmpdir), percent_models=100, model_count=2)
            Path(cfg.hitsoutdir).write_text(
                "GenomeSpurious|c1|p1 - 100 M1 - 100 1e-10 50.0\n"
                "GenomeSpurious|c1|p1 - 100 M2 - 100 1e-20 100.0\n"
                "GenomeComplete|c1|p1 - 100 M1 - 100 1e-20 100.0\n"
                "GenomeComplete|c1|p2 - 100 M2 - 100 1e-20 100.0\n"
            )

            finaldf, counts = parse_hmmsearch(cfg)

            self.assertEqual(
                set(finaldf[0]), {"GenomeComplete|c1|p1", "GenomeComplete|c1|p2"}
            )
            self.assertEqual(counts["GenomeSpurious"], {"M2": 1})
            count_matrix = pd.read_csv(
                Path(cfg.outdir) / "marker_count_matrix.csv", index_col=0
            )
            self.assertEqual(list(count_matrix.columns), ["GenomeComplete"])
            self.assertEqual(
                count_matrix["GenomeComplete"].to_dict(), {"M1": 1, "M2": 1}
            )
            self.assertEqual(
                (Path(cfg.outdir) / "log_genomes_removed.txt").read_text(),
                "GenomeSpurious\tminmarker:0.5000\n",
            )


if __name__ == "__main__":
    unittest.main()
