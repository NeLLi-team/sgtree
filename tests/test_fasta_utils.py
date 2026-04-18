import tempfile
import unittest
from pathlib import Path

from sgtree._fasta_utils import fasta_contig_bases_stats


class FastaContigBasesStatsTests(unittest.TestCase):
    def test_counts_contigs_and_bases(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "genome.fna"
            path.write_text(
                ">contig1\nACGT\n>contig2\nAAAACCCC\n>contig3\nGG\n"
            )
            contigs, bases = fasta_contig_bases_stats(path)
            self.assertEqual(contigs, 3)
            self.assertEqual(bases, 4 + 8 + 2)

    def test_accepts_str_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "g.faa"
            path.write_text(">p1\nMKT\n")
            self.assertEqual(fasta_contig_bases_stats(str(path)), (1, 3))

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "empty.fna"
            path.write_text("")
            self.assertEqual(fasta_contig_bases_stats(path), (0, 0))


if __name__ == "__main__":
    unittest.main()
