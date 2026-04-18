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

    def test_build_supermatrix_artifact_equivalence(self):
        # 3 genomes x 3 markers with missing cells. Pins: table CSV columns,
        # concat FASTA record ids and per-record length, X-gap fill widths per
        # marker, column order (alphabetical after SeqID).
        from Bio import SeqIO

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            trimmed_dir = tmp / "trimmed"
            trimmed_dir.mkdir()
            (trimmed_dir / "M1.faa").write_text(">g1|p1\nAAAA\n>g2|p2\nCCCC\n")
            (trimmed_dir / "M2.faa").write_text(">g1|p1\nGGGGG\n>g3|p1\nTTTTT\n")
            (trimmed_dir / "M3.faa").write_text(">g2|p1\nNN\n>g3|p1\nMM\n")

            table_path = str(tmp / "table.csv")
            concat_path = str(tmp / "concat.faa")
            build_supermatrix(
                str(trimmed_dir), str(tmp / "out"), table_path, concat_path
            )

            table = pd.read_csv(table_path).set_index("SeqID")
            self.assertEqual(sorted(table.index), ["g1", "g2", "g3"])
            marker_cols = [c for c in table.columns if c.endswith(".faa")]
            self.assertEqual(marker_cols, ["M1.faa", "M2.faa", "M3.faa"])

            def cell(genome: str, marker: str) -> str:
                return str(table.loc[genome, marker]).replace("\n", "")

            # Real cells preserved (ignore legacy trailing newline inside cells)
            self.assertEqual(cell("g1", "M1.faa"), "AAAA")
            self.assertEqual(cell("g2", "M1.faa"), "CCCC")
            self.assertEqual(cell("g1", "M2.faa"), "GGGGG")
            self.assertEqual(cell("g3", "M2.faa"), "TTTTT")

            # NaN cells filled with X-gaps of the column's alignment width
            self.assertEqual(cell("g3", "M1.faa"), "X" * 4)
            self.assertEqual(cell("g2", "M2.faa"), "X" * 5)
            self.assertEqual(cell("g1", "M3.faa"), "X" * 2)

            with open(concat_path) as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            self.assertEqual(sorted(r.id for r in records), ["g1", "g2", "g3"])
            # Each genome's supermatrix row has width 4 + 5 + 2 = 11.
            for rec in records:
                self.assertEqual(len(rec.seq), 11)

            seqs = {r.id: str(r.seq) for r in records}
            self.assertEqual(seqs["g1"], "AAAA" + "GGGGG" + "X" * 2)
            self.assertEqual(seqs["g2"], "CCCC" + "X" * 5 + "NN")
            self.assertEqual(seqs["g3"], "X" * 4 + "TTTTT" + "MM")

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
