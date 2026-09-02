import random
import unittest

from sgtree.benchmarks.sequence_evidence import (
    PHMMER_MAX_EVALUE,
    PHMMER_MAX_QUERY_COUNT,
    PHMMER_MIN_QUERY_COVERAGE,
    PHMMER_MIN_QUERY_COUNT,
    PHMMER_MIN_SCORE_MARGIN,
    assign_contig_gene_split_votes,
)
from sgtree.marker_selection.contig_evidence import contig_gene_vote_gate


class PhmmerSplitVoteTests(unittest.TestCase):
    @staticmethod
    def _protein(seed: int) -> str:
        rng = random.Random(seed)
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        return "".join(rng.choice(alphabet) for _ in range(150))

    @staticmethod
    def _mutate(sequence: str, every: int = 10) -> str:
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        mutated = list(sequence)
        for index in range(0, len(mutated), every):
            mutated[index] = alphabet[(alphabet.index(mutated[index]) + 7) % len(alphabet)]
        return "".join(mutated)

    def _queries(self) -> dict[str, str]:
        return {
            f"R|ctg|g{index}": self._protein(index)
            for index in range(1, PHMMER_MIN_QUERY_COUNT + 1)
        }

    @staticmethod
    def _records(genome: str, sequences: list[str], suffix: str = "ref") -> dict[str, str]:
        return {
            f"{genome}|{suffix}|g{index}": sequence
            for index, sequence in enumerate(sequences, start=1)
        }

    def _score(self, queries: dict[str, str], panel: dict, **overrides) -> dict:
        arguments = {
            "recipient_genome": "R",
            "candidate_contig_id": "ctg",
            "marker_record_ids": set(),
            "attachment_taxa": {"A"},
            "background_taxa": {"B"},
        }
        arguments.update(overrides)
        return assign_contig_gene_split_votes(queries, panel, **arguments)

    def test_attachment_side_wins_by_fixed_margin_and_collapses_paralogs(self):
        queries = self._queries()
        sequences = list(queries.values())
        panel = {
            "A": {
                **self._records("A", sequences),
                **self._records("A", [self._mutate(value, 5) for value in sequences], "paralog"),
            },
            "B": self._records("B", [self._mutate(value) for value in sequences]),
        }

        report = self._score(queries, panel)

        self.assertEqual(PHMMER_MAX_EVALUE, 1e-5)
        self.assertEqual(PHMMER_MIN_QUERY_COVERAGE, 0.5)
        self.assertEqual(PHMMER_MIN_SCORE_MARGIN, 5.0)
        self.assertEqual(PHMMER_MIN_QUERY_COUNT, 3)
        self.assertEqual(PHMMER_MAX_QUERY_COUNT, 12)
        self.assertEqual(report["informative_vote_count"], 3)
        self.assertEqual({vote["assigned_clade"] for vote in report["votes"]}, {"A"})
        self.assertTrue(
            all(vote["score_margin"] >= PHMMER_MIN_SCORE_MARGIN for vote in report["votes"])
        )
        self.assertTrue(
            all(vote["attachment_qualifying_genome_count"] == 1 for vote in report["votes"])
        )
        self.assertTrue(contig_gene_vote_gate(report["votes"], "A")["contig_gate_pass"])

    def test_background_side_becomes_deterministic_conflict_vote(self):
        queries = self._queries()
        sequences = list(queries.values())
        panel = {
            "A": self._records("A", [self._mutate(value) for value in sequences]),
            "B": self._records("B", sequences),
        }

        report = self._score(queries, panel)

        self.assertEqual(report["background_clade"], "complement:B")
        self.assertEqual(report["background_vote_count"], 3)
        self.assertEqual(
            {vote["assigned_clade"] for vote in report["votes"]},
            {"complement:B"},
        )
        self.assertFalse(contig_gene_vote_gate(report["votes"], "A")["contig_gate_pass"])

    def test_near_ties_are_uninformative(self):
        queries = self._queries()
        sequences = list(queries.values())
        panel = {
            "A": self._records("A", sequences),
            "B": self._records("B", sequences),
        }

        report = self._score(queries, panel)

        self.assertEqual(report["informative_vote_count"], 0)
        self.assertEqual(
            {vote["match_status"] for vote in report["votes"]},
            {"score_margin_below_threshold"},
        )
        self.assertTrue(
            all(abs(vote["score_margin"]) < PHMMER_MIN_SCORE_MARGIN for vote in report["votes"])
        )

    def test_invalid_or_missing_split_side_fails_closed(self):
        queries = self._queries()
        panel = {"A": self._records("A", list(queries.values()))}
        cases = (
            ({"attachment_taxa": set()}, "invalid_attachment_taxa"),
            ({"background_taxa": set()}, "invalid_background_taxa"),
            ({"attachment_taxa": {"A"}, "background_taxa": {"A"}}, "overlapping_split_taxa"),
            ({"attachment_taxa": {"R"}, "background_taxa": {"B"}}, "recipient_in_split_taxa"),
            ({}, "missing_split_reference_side"),
        )
        for overrides, expected in cases:
            with self.subTest(expected=expected):
                report = self._score(queries, panel, **overrides)
                self.assertEqual(report["input_status"], expected)
                self.assertEqual(report["informative_vote_count"], 0)
                self.assertTrue(all(not vote["informative"] for vote in report["votes"]))

    def test_requires_at_least_three_non_marker_queries_and_caps_work(self):
        queries = self._queries()
        sequences = list(queries.values())
        panel = {
            "A": self._records("A", sequences),
            "B": self._records("B", [self._mutate(value) for value in sequences]),
        }

        two_queries = dict(list(queries.items())[:2])
        wrong_count = self._score(two_queries, panel)
        self.assertEqual(
            wrong_count["input_status"],
            "requires_at_least_three_non_marker_queries",
        )

        many_queries = {
            f"R|ctg|g{index}": self._protein(index)
            for index in range(1, PHMMER_MAX_QUERY_COUNT + 2)
        }
        many_sequences = list(many_queries.values())
        many_panel = {
            "A": self._records("A", many_sequences),
            "B": self._records(
                "B",
                [self._mutate(value) for value in many_sequences],
            ),
        }
        bounded = self._score(many_queries, many_panel)
        self.assertEqual(bounded["input_status"], "ok")
        self.assertEqual(
            bounded["eligible_non_marker_query_count"],
            PHMMER_MAX_QUERY_COUNT + 1,
        )
        self.assertEqual(bounded["query_count"], PHMMER_MAX_QUERY_COUNT)
        self.assertEqual(bounded["query_selection_truncated_count"], 1)
        self.assertEqual(bounded["informative_vote_count"], PHMMER_MAX_QUERY_COUNT)

    def test_candidate_genes_must_share_recipient_contig_and_clean_id(self):
        queries = self._queries()
        sequences = list(queries.values())
        panel = {
            "A": self._records("A", sequences),
            "B": self._records("B", [self._mutate(value) for value in sequences]),
        }
        invalid_cases = (
            (
                {**queries, "R|other|g4": self._protein(4)},
                "candidate_contig_mismatch",
            ),
            (
                {**queries, "X|ctg|g4": self._protein(4)},
                "candidate_recipient_mismatch",
            ),
            (
                {
                    **queries,
                    " R|ctg|g1 ": queries["R|ctg|g1"],
                },
                "duplicate_candidate_record_id",
            ),
        )
        for candidate_genes, expected in invalid_cases:
            with self.subTest(expected=expected):
                report = self._score(candidate_genes, panel)
                self.assertEqual(report["input_status"], expected)
                self.assertEqual(report["informative_vote_count"], 0)
                self.assertTrue(
                    all(not vote["informative"] for vote in report["votes"])
                )

    def test_recipient_hits_are_excluded_when_source_is_absent(self):
        queries = self._queries()
        panel = {
            "R": self._records("R", list(queries.values())),
            "A": self._records("A", [value[:60] for value in queries.values()]),
            "B": self._records("B", [self._protein(200 + index) for index in range(3)]),
        }

        report = self._score(queries, panel)

        self.assertEqual(report["excluded_recipient_reference_count"], 3)
        self.assertEqual(report["informative_vote_count"], 0)
        self.assertEqual(
            {vote["match_status"] for vote in report["votes"]},
            {"missing_split_side_hit"},
        )
        self.assertTrue(
            all(vote["coverage_filtered_hit_count"] >= 1 for vote in report["votes"])
        )

    def test_marker_and_candidate_self_records_are_excluded(self):
        queries = self._queries()
        sequences = list(queries.values())
        marker_id = "R|ctg|marker"
        panel_a = {
            **queries,
            marker_id: sequences[0],
            **self._records("A", [self._mutate(value) for value in sequences]),
        }
        panel = {
            "A": panel_a,
            "B": self._records("B", sequences),
        }

        report = self._score(
            {**queries, marker_id: sequences[0]},
            panel,
            marker_record_ids={marker_id},
        )

        self.assertEqual(report["query_count"], 3)
        self.assertEqual(report["candidate_marker_record_count"], 1)
        self.assertEqual(report["excluded_candidate_reference_count"], 3)
        self.assertEqual(report["excluded_marker_reference_count"], 1)
        self.assertEqual(report["background_vote_count"], 3)
        self.assertEqual(
            {vote["assigned_clade"] for vote in report["votes"]},
            {"complement:B"},
        )


if __name__ == "__main__":
    unittest.main()
