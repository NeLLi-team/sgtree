import tempfile
import unittest
from pathlib import Path

from sgtree.supermatrix import build_supermatrix


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


if __name__ == "__main__":
    unittest.main()
