import tempfile
import unittest
from pathlib import Path

from sgtree._fasta_utils import fasta_contig_bases_stats
from sgtree.fasta_normalize import normalize_and_concat_proteomes
from sgtree.id_schema import parse_sequence_id, sanitize_token


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


class SanitizeTokenTests(unittest.TestCase):
    def test_token_table(self):
        cases = [
            # token, expected
            ("", "FALLBACK"),
            (" ", "FALLBACK"),
            ("\t\n ", "FALLBACK"),
            (None, "FALLBACK"),
            ("!!!", "FALLBACK"),
            ("gene_1", "gene_1"),
            ("GCF_000005845.2", "GCF_000005845.2"),
            ("NC_000913.3 Escherichia coli", "NC_000913.3"),
            ("a|b/c", "a_b_c"),
        ]
        for token, expected in cases:
            with self.subTest(token=token):
                self.assertEqual(sanitize_token(token, "FALLBACK"), expected)


class HeaderWithEmptyIdFieldTests(unittest.TestCase):
    def test_headers_with_empty_fields_normalize_instead_of_crashing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            (tmp / "genome.faa").write_text(
                ">GENOME|\nMKT\n"
                ">a||b\nMKT\n"
                ">GENOME|c1|\nMKT\n"
                ">\nMKT\n"
            )
            out_fasta = tmp / "proteomes"

            stats = normalize_and_concat_proteomes(str(tmp), str(out_fasta))

            self.assertEqual(stats["records"], 4)
            ids = [
                line[1:].strip()
                for line in out_fasta.read_text().splitlines()
                if line.startswith(">")
            ]
            self.assertEqual(
                ids,
                [
                    "genome|unknown_contig|protein_000001",
                    "genome|unknown_contig|b",
                    "genome|c1|protein_000003",
                    "genome|unknown_contig|protein_000004",
                ],
            )
            for identifier in ids:
                genome, contig, gene = parse_sequence_id(identifier)
                self.assertTrue(all([genome, contig, gene]), identifier)


if __name__ == "__main__":
    unittest.main()
