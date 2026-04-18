import unittest

from sgtree import marker_selection
from sgtree.config import DEFAULT_SINGLETON_THRESHOLDS, SingletonThresholds


class SingletonThresholdsTests(unittest.TestCase):
    def test_defaults_match_legacy_module_constants(self):
        t = DEFAULT_SINGLETON_THRESHOLDS
        self.assertEqual(t.composite_score_threshold, marker_selection.COMPOSITE_SCORE_THRESHOLD)
        self.assertEqual(t.composite_score_margin, marker_selection.COMPOSITE_SCORE_MARGIN)
        self.assertEqual(t.recipient_consensus_z_threshold, marker_selection.RECIPIENT_CONSENSUS_Z_THRESHOLD)
        self.assertEqual(t.recipient_consensus_min_score, marker_selection.RECIPIENT_CONSENSUS_MIN_SCORE)
        self.assertEqual(t.recipient_consensus_rank_margin, marker_selection.RECIPIENT_CONSENSUS_RANK_MARGIN)
        self.assertEqual(t.contig_consensus_high_overlap, marker_selection.CONTIG_CONSENSUS_HIGH_OVERLAP)
        self.assertEqual(t.contig_consensus_low_overlap, marker_selection.CONTIG_CONSENSUS_LOW_OVERLAP)
        self.assertEqual(t.contig_consensus_min_support, marker_selection.CONTIG_CONSENSUS_MIN_SUPPORT)
        self.assertEqual(t.neighbor_clade_min_present, marker_selection.NEIGHBOR_CLADE_MIN_PRESENT)
        self.assertEqual(t.neighbor_clade_min_present_fraction, marker_selection.NEIGHBOR_CLADE_MIN_PRESENT_FRACTION)
        self.assertEqual(t.neighbor_clade_min_species_anchor, marker_selection.NEIGHBOR_CLADE_MIN_SPECIES_ANCHOR)
        self.assertEqual(t.neighbor_clade_min_neighbor_purity, marker_selection.NEIGHBOR_CLADE_MIN_NEIGHBOR_PURITY)
        self.assertEqual(t.neighbor_clade_min_purity_drop, marker_selection.NEIGHBOR_CLADE_MIN_PURITY_DROP)
        self.assertEqual(t.neighbor_clade_max_knn_agreement, marker_selection.NEIGHBOR_CLADE_MAX_KNN_AGREEMENT)
        self.assertEqual(t.neighbor_clade_min_score, marker_selection.NEIGHBOR_CLADE_MIN_SCORE)
        self.assertEqual(t.neighbor_clade_min_recipient_support, marker_selection.NEIGHBOR_CLADE_MIN_RECIPIENT_SUPPORT)
        self.assertEqual(t.neighbor_ml_high_conf_min_present, marker_selection.NEIGHBOR_ML_HIGH_CONF_MIN_PRESENT)
        self.assertEqual(t.neighbor_ml_high_conf_min_fraction, marker_selection.NEIGHBOR_ML_HIGH_CONF_MIN_FRACTION)
        self.assertEqual(t.neighbor_ml_high_conf_min_anchor, marker_selection.NEIGHBOR_ML_HIGH_CONF_MIN_ANCHOR)
        self.assertEqual(t.unknown_contig_delta_floor, marker_selection.UNKNOWN_CONTIG_DELTA_FLOOR)
        self.assertEqual(t.unknown_contig_topoknn_floor, marker_selection.UNKNOWN_CONTIG_TOPOKNN_FLOOR)
        self.assertEqual(t.unknown_contig_recipient_floor, marker_selection.UNKNOWN_CONTIG_RECIPIENT_FLOOR)
        self.assertEqual(t.gcp_combined_threshold, marker_selection.GCP_COMBINED_THRESHOLD)
        self.assertEqual(t.gcp_outlier_count_threshold, marker_selection.GCP_OUTLIER_COUNT_THRESHOLD)
        self.assertEqual(t.gcp_iforest_threshold, marker_selection.GCP_IFOREST_THRESHOLD)
        self.assertEqual(t.gcp_z_threshold, marker_selection.GCP_Z_THRESHOLD)
        self.assertEqual(t.gcp_min_genomes, marker_selection.GCP_MIN_GENOMES)
        self.assertEqual(t.gcp_min_markers, marker_selection.GCP_MIN_MARKERS)

    def test_instance_is_frozen(self):
        t = DEFAULT_SINGLETON_THRESHOLDS
        with self.assertRaises(Exception):
            # dataclass(frozen=True) raises FrozenInstanceError on mutation
            t.composite_score_threshold = 99.0  # type: ignore[misc]

    def test_overrides_via_constructor(self):
        t = SingletonThresholds(gcp_min_genomes=10, neighbor_clade_min_score=3.5)
        self.assertEqual(t.gcp_min_genomes, 10)
        self.assertEqual(t.neighbor_clade_min_score, 3.5)
        # Other values remain at defaults.
        self.assertEqual(t.composite_score_threshold, 2.0)


if __name__ == "__main__":
    unittest.main()
