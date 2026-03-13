import unittest

from sgtree.ani import (
    GenomeRecord,
    _merge_intervals,
    choose_cluster_representative,
)


class AniClusteringTests(unittest.TestCase):
    def test_choose_cluster_representative_prefers_reference_and_lower_contig_count(self):
        representative = choose_cluster_representative(
            [
                GenomeRecord("QueryA", "query", "fna", "QueryA.fna", "QueryA.fna", None, 1000, 5),
                GenomeRecord("RefA", "ref", "fna", "RefA.fna", "RefA.fna", None, 900, 8),
                GenomeRecord("RefB", "ref", "fna", "RefB.fna", "RefB.fna", None, 850, 2),
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
