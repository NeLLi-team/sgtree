import unittest

from sgtree.marker_selection.contig_evidence import contig_gene_vote_gate


def _vote(
    gene_id: str,
    clade: str,
    *,
    informative: bool = True,
) -> dict:
    return {
        "gene_id": gene_id,
        "assigned_clade": clade,
        "informative": informative,
    }


class ContigGeneVoteGateTests(unittest.TestCase):
    def test_three_matching_genes_pass(self):
        result = contig_gene_vote_gate(
            [_vote("g1", "donor"), _vote("g2", "donor"), _vote("g3", "donor")],
            "donor",
        )

        self.assertTrue(result["contig_gate_pass"])
        self.assertIsNone(result["contig_abstention_reason"])
        self.assertEqual(result["informative_gene_count"], 3)
        self.assertEqual(result["agreement_count"], 3)
        self.assertEqual(result["agreement_fraction"], 1.0)

    def test_exact_two_thirds_passes_without_strong_conflict(self):
        result = contig_gene_vote_gate(
            [_vote("g1", "donor"), _vote("g2", "donor"), _vote("g3", "other")],
            "donor",
        )

        self.assertTrue(result["contig_gate_pass"])
        self.assertEqual(result["agreement_count"], 2)
        self.assertEqual(result["agreement_fraction"], 2 / 3)

    def test_coherent_conflicting_clade_vetoes_two_thirds_agreement(self):
        result = contig_gene_vote_gate(
            [
                _vote("g1", "donor"),
                _vote("g2", "donor"),
                _vote("g3", "donor"),
                _vote("g4", "donor"),
                _vote("g5", "other"),
                _vote("g6", "other"),
            ],
            "donor",
        )

        self.assertFalse(result["contig_gate_pass"])
        self.assertEqual(result["contig_abstention_reason"], "strong_conflicting_clade")
        self.assertEqual(result["strongest_conflict_clade"], "other")
        self.assertEqual(result["strong_conflict_count"], 2)

    def test_solo_marker_contig_abstains(self):
        result = contig_gene_vote_gate([], "donor")

        self.assertFalse(result["contig_gate_pass"])
        self.assertEqual(
            result["contig_abstention_reason"],
            "insufficient_informative_genes",
        )

    def test_native_contig_votes_do_not_confirm_donor(self):
        result = contig_gene_vote_gate(
            [
                _vote("g1", "recipient"),
                _vote("g2", "recipient"),
                _vote("g3", "recipient"),
            ],
            "donor",
        )

        self.assertFalse(result["contig_gate_pass"])
        self.assertEqual(result["top_clade"], "recipient")
        self.assertEqual(result["contig_abstention_reason"], "strong_conflicting_clade")

    def test_weak_split_vote_abstains_below_two_thirds(self):
        result = contig_gene_vote_gate(
            [
                _vote("g1", "donor"),
                _vote("g2", "other"),
                _vote("g3", "third"),
            ],
            "donor",
        )

        self.assertFalse(result["contig_gate_pass"])
        self.assertEqual(
            result["contig_abstention_reason"],
            "insufficient_donor_agreement",
        )

    def test_noninformative_rows_do_not_count(self):
        result = contig_gene_vote_gate(
            [
                _vote("g1", "donor"),
                _vote("g2", "donor"),
                _vote("g3", "donor", informative=False),
            ],
            "donor",
        )

        self.assertEqual(result["informative_gene_count"], 2)
        self.assertEqual(
            result["contig_abstention_reason"],
            "insufficient_informative_genes",
        )

    def test_null_assignment_fails_closed_without_becoming_literal_none(self):
        result = contig_gene_vote_gate(
            [
                _vote("g1", "donor"),
                _vote("g2", "donor"),
                {"gene_id": "g3", "assigned_clade": None, "informative": True},
            ],
            "donor",
        )

        self.assertEqual(result["informative_gene_count"], 2)
        self.assertEqual(result["invalid_gene_vote_count"], 1)
        self.assertEqual(result["contig_abstention_reason"], "invalid_gene_votes")
        self.assertNotEqual(result["top_clade"], "None")

    def test_duplicate_gene_ids_fail_closed(self):
        result = contig_gene_vote_gate(
            [
                _vote("g1", "donor"),
                _vote("g1", "donor"),
                _vote("g2", "donor"),
                _vote("g3", "donor"),
            ],
            "donor",
        )

        self.assertEqual(result["informative_gene_count"], 3)
        self.assertEqual(result["duplicate_gene_vote_count"], 1)
        self.assertFalse(result["contig_gate_pass"])
        self.assertEqual(result["contig_abstention_reason"], "duplicate_gene_votes")

    def test_vote_order_does_not_change_result(self):
        votes = [
            _vote("g1", "donor"),
            _vote("g2", "other"),
            _vote("g3", "third"),
        ]

        self.assertEqual(
            contig_gene_vote_gate(votes, "donor"),
            contig_gene_vote_gate(reversed(votes), "donor"),
        )

    def test_missing_proposed_clade_fails_closed(self):
        result = contig_gene_vote_gate(
            [_vote("g1", "donor"), _vote("g2", "donor"), _vote("g3", "donor")],
            None,
        )

        self.assertFalse(result["contig_gate_pass"])
        self.assertEqual(result["contig_abstention_reason"], "missing_proposed_clade")


if __name__ == "__main__":
    unittest.main()
