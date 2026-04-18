import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from sgtree.supermatrix import (
    _fill_nan_gaps,
    _trim_alignment_fallback,
    _trimal_simple_worker,
    build_supermatrix,
)


class SupermatrixTests(unittest.TestCase):
    def test_build_supermatrix_rejects_duplicate_genome_ids_in_alignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            trimmed_dir = tmp / "trimmed"
            trimmed_dir.mkdir()

            (trimmed_dir / "M1.faa").write_text(
                ">GenomeA|native\nAAAA\n>GenomeA|contam\nCCCC\n"
            )
            (trimmed_dir / "M2.faa").write_text(
                ">GenomeA|native\nGGGG\n>GenomeB|native\nTTTT\n"
            )

            with self.assertRaises(ValueError):
                build_supermatrix(
                    str(trimmed_dir),
                    str(tmp / "out"),
                    str(tmp / "table.csv"),
                    str(tmp / "concat.faa"),
                )

    def test_trim_alignment_fallback_drops_all_gap_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.faa"
            output_path = tmp / "output.faa"
            input_path.write_text(
                ">A\nA-\n>B\n--\n>C\nG-\n"
            )

            _trim_alignment_fallback(str(input_path), str(output_path), gap_threshold=0.1)

            self.assertEqual(
                output_path.read_text(),
                ">A\nA\n>B\n-\n>C\nG\n",
            )

    def test_fill_nan_gaps_preserves_alignment_widths(self):
        # Col M1: row widths 100, NaN, 100, NaN — forward pass only.
        # Col M2: NaN, NaN, 150, 150 — leading NaNs exercise the backward pass,
        # which uses the BOTTOM-most non-NaN row as ref_len (quirk preserved).
        df = pd.DataFrame(
            {
                "SeqID": ["g1", "g2", "g3", "g4"],
                "M1": ["A" * 100, np.nan, "B" * 100, np.nan],
                "M2": [np.nan, np.nan, "C" * 150, "D" * 150],
            }
        )

        _fill_nan_gaps(df)

        self.assertEqual(df.loc[0, "M1"], "A" * 100)
        self.assertEqual(df.loc[1, "M1"], "X" * 100)
        self.assertEqual(df.loc[2, "M1"], "B" * 100)
        self.assertEqual(df.loc[3, "M1"], "X" * 100)

        self.assertEqual(df.loc[0, "M2"], "X" * 150)
        self.assertEqual(df.loc[1, "M2"], "X" * 150)
        self.assertEqual(df.loc[2, "M2"], "C" * 150)
        self.assertEqual(df.loc[3, "M2"], "D" * 150)

        for col in ("M1", "M2"):
            widths = df[col].map(lambda s: len(s.replace("\n", "")))
            self.assertEqual(widths.nunique(), 1, f"column {col} has mixed widths")

    def test_fill_nan_gaps_uses_last_non_nan_width_for_leading_gaps(self):
        # Backward pass quirk: leading NaNs are filled with the BOTTOM-most
        # non-NaN row's width, not the immediately-following one.
        df = pd.DataFrame(
            {
                "SeqID": ["g1", "g2", "g3"],
                "M1": [np.nan, "A" * 50, "B" * 200],
            }
        )

        _fill_nan_gaps(df)

        # row 0 uses g3's width (200), not g2's (50).
        self.assertEqual(df.loc[0, "M1"], "X" * 200)

    def test_trimal_simple_worker_falls_back_when_trimal_is_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.faa"
            output_path = tmp / "output.faa"
            input_path.write_text(
                ">A\nA-\n>B\n--\n>C\nG-\n"
            )

            with patch("sgtree.supermatrix.subprocess.run", side_effect=FileNotFoundError()):
                _trimal_simple_worker((str(input_path), str(output_path)))

            self.assertEqual(
                output_path.read_text(),
                ">A\nA\n>B\n-\n>C\nG\n",
            )


if __name__ == "__main__":
    unittest.main()
