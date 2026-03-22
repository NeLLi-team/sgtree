import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sgtree.supermatrix import _trim_alignment_fallback, _trimal_simple_worker, build_supermatrix


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
