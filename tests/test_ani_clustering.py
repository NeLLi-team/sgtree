import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from sgtree.ani import (
    GenomeRecord,
    _merge_intervals,
    choose_cluster_representative,
)
from sgtree.ani_clustering import _discover_inputs


class AniClusteringTests(unittest.TestCase):
    def test_input_discovery_ignores_hidden_sidecars(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample.fna").write_text(">contig\nACGT\n", encoding="utf-8")
            (root / ".metadata").write_text("sample metadata\n", encoding="utf-8")

            inputs = _discover_inputs(str(root), "query")

        self.assertEqual([row["genome_id"] for row in inputs], ["sample"])

    def test_input_discovery_keeps_content_detected_fasta_files(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for filename in ("sample_a.fna", "sample_b.fsa"):
                (root / filename).write_text(">contig\nACGTACGT\n", encoding="utf-8")

            inputs = _discover_inputs(str(root), "query")

        self.assertEqual([row["genome_id"] for row in inputs], ["sample_a", "sample_b"])

    def test_input_discovery_rejects_mixed_fasta_formats(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_a.fna").write_text(">contig\nACGT\n", encoding="utf-8")
            (root / "sample_b.faa").write_text(">protein\nMPEPTIDE\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Mixed nucleotide and protein"):
                _discover_inputs(str(root), "query")

    def test_input_discovery_rejects_non_fasta_files(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample_a.fna").write_text(">contig\nACGT\n", encoding="utf-8")
            (root / "sample_b.txt").write_text("ACGT\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "FASTA header"):
                _discover_inputs(str(root), "query")

    def test_input_discovery_rejects_empty_fasta(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "sample.fna").write_text(">contig\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "no sequence"):
                _discover_inputs(str(root), "query")

    def test_input_discovery_preserves_symlink_sample_names(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            assembly = root / "original.fna"
            assembly.write_text(">contig\nACGTACGT\n", encoding="utf-8")
            queries = root / "queries"
            queries.mkdir()
            (queries / "sample_a.fna").symlink_to(assembly)
            (queries / "sample_b.FNA").symlink_to(assembly)

            inputs = _discover_inputs(str(queries), "query")

        self.assertEqual(len(inputs), 2)
        self.assertEqual(inputs[0]["genome_id"], "sample_a")
        self.assertEqual(inputs[1]["genome_id"], "sample_b")

    def test_choose_cluster_representative_prefers_reference_and_lower_contig_count(
        self,
    ):
        representative = choose_cluster_representative(
            [
                GenomeRecord(
                    "QueryA", "query", "fna", "QueryA.fna", "QueryA.fna", None, 1000, 5
                ),
                GenomeRecord(
                    "RefA", "ref", "fna", "RefA.fna", "RefA.fna", None, 900, 8
                ),
                GenomeRecord(
                    "RefB", "ref", "fna", "RefB.fna", "RefB.fna", None, 850, 2
                ),
            ]
        )
        self.assertEqual(representative.genome_id, "RefB")

    def test_merge_intervals_normalizes_overlapping_ranges(self):
        self.assertEqual(
            _merge_intervals([(5, 10), (1, 3), (2, 6), (12, 20), (18, 21)]),
            [(1, 10), (12, 21)],
        )


if __name__ == "__main__":
    unittest.main()
