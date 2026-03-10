import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sgtree import benchmark


class BenchmarkTests(unittest.TestCase):
    def test_make_contaminant_record_rehomes_donor_under_recipient(self):
        donor = SeqRecord(Seq("MPEPTIDE"), id="DonorA|prot123", description="DonorA|prot123")

        record = benchmark.make_contaminant_record(
            recipient_genome="RecipientB",
            donor_record=donor,
            marker="MarkerX",
            donor_genome="DonorA",
            event_index=1,
        )

        self.assertTrue(record.id.startswith("RecipientB|contam__MarkerX__DonorA__e001"))
        self.assertEqual(str(record.seq), "MPEPTIDE")

    def test_apply_replacement_event_removes_native_and_adds_contaminant(self):
        recipient_records = {
            "RecipientB|native_marker": SeqRecord(
                Seq("AAAA"),
                id="RecipientB|native_marker",
                description="RecipientB|native_marker",
            ),
            "RecipientB|background": SeqRecord(
                Seq("TTTT"),
                id="RecipientB|background",
                description="RecipientB|background",
            ),
        }
        contaminant = SeqRecord(
            Seq("CCCC"),
            id="RecipientB|contam__MarkerX__DonorA__e001",
            description="RecipientB|contam__MarkerX__DonorA__e001",
        )

        updated = benchmark.apply_replacement_event(
            recipient_records,
            native_record_id="RecipientB|native_marker",
            contaminant_record=contaminant,
        )

        self.assertNotIn("RecipientB|native_marker", updated)
        self.assertIn("RecipientB|background", updated)
        self.assertIn("RecipientB|contam__MarkerX__DonorA__e001", updated)

    def test_drop_native_marker_removes_record_without_replacement(self):
        recipient_records = {
            "RecipientB|native_marker": SeqRecord(
                Seq("AAAA"),
                id="RecipientB|native_marker",
                description="RecipientB|native_marker",
            ),
            "RecipientB|background": SeqRecord(
                Seq("TTTT"),
                id="RecipientB|background",
                description="RecipientB|background",
            ),
        }

        updated = benchmark.drop_native_marker(
            recipient_records,
            native_record_id="RecipientB|native_marker",
        )

        self.assertNotIn("RecipientB|native_marker", updated)
        self.assertIn("RecipientB|background", updated)

    def test_generate_cross_benchmark_accepts_different_model_sets_when_markers_overlap(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_a_dir = root / "bench_a"
            benchmark_b_dir = root / "bench_b"
            outdir = root / "cross"
            benchmark_a_dir.mkdir()
            benchmark_b_dir.mkdir()

            (benchmark_a_dir / "benchmark_manifest.json").write_text(
                json.dumps(
                    {
                        "models_path": "resources/models/UNI56.hmm",
                        "marker_ranking": [
                            {"marker": "COG0090"},
                            {"marker": "COG0092"},
                        ],
                    }
                )
            )
            (benchmark_b_dir / "benchmark_manifest.json").write_text(
                json.dumps(
                    {
                        "models_path": "resources/models/RProt16.hmm",
                        "marker_ranking": [
                            {"marker": "COG0092"},
                            {"marker": "COG0094"},
                        ],
                    }
                )
            )

            truth_records_a = {
                "FlavoA": {
                    "FlavoA|prot1": SeqRecord(
                        Seq("MPEPTIDE"),
                        id="FlavoA|prot1",
                        description="FlavoA|prot1",
                    )
                }
            }
            truth_records_b = {
                "GammaB": {
                    "GammaB|prot1": SeqRecord(
                        Seq("MPEPTIDE"),
                        id="GammaB|prot1",
                        description="GammaB|prot1",
                    )
                }
            }
            captured = {}

            def fake_materialize(**kwargs):
                captured.update(kwargs)

            with (
                patch.object(
                    benchmark,
                    "_read_normalized_proteomes",
                    side_effect=[truth_records_a, truth_records_b],
                ),
                patch.object(
                    benchmark,
                    "_materialize_benchmark_from_truth",
                    side_effect=fake_materialize,
                ),
            ):
                benchmark.generate_cross_benchmark_dataset(
                    benchmark_a_dir=benchmark_a_dir,
                    benchmark_b_dir=benchmark_b_dir,
                    outdir=outdir,
                    n_markers=2,
                    seed=42,
                    num_cpus=8,
                )

        self.assertEqual(captured["truth_markers"], ["COG0092"])
        self.assertEqual(captured["models_path"], Path("resources/models/UNI56.hmm"))
        self.assertEqual(captured["group_labels"], {"FlavoA": "flavo", "GammaB": "gamma"})
        self.assertTrue(captured["cross_group_only"])

    def test_read_normalized_proteomes_accepts_directory_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            proteomes_dir = Path(tmpdir) / "truth_inputs"
            proteomes_dir.mkdir()
            (proteomes_dir / "GenomeA.faa").write_text(
                ">GenomeA|prot1\nMPEPTIDE\n>GenomeA|prot2\nMSEQ\n"
            )
            (proteomes_dir / "GenomeB.faa").write_text(
                ">GenomeB|prot1\nMOTHER\n"
            )

            records = benchmark._read_normalized_proteomes(proteomes_dir)

        self.assertEqual(set(records), {"GenomeA", "GenomeB"})
        self.assertEqual(set(records["GenomeA"]), {"GenomeA|prot1", "GenomeA|prot2"})
        self.assertEqual(set(records["GenomeB"]), {"GenomeB|prot1"})

    def test_write_genome_summary_tsv_counts_duplicate_and_replacement_events(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "genome_summary.tsv"
            summary = benchmark._write_genome_summary_tsv(
                out_path,
                [
                    {
                        "event_index": 1,
                        "event_type": "duplicate",
                        "recipient_genome": "GenomeA",
                        "recipient_group": "flavo",
                        "source_relation": "within_group",
                    },
                    {
                        "event_index": 2,
                        "event_type": "replacement",
                        "recipient_genome": "GenomeA",
                        "recipient_group": "flavo",
                        "source_relation": "cross_group",
                    },
                    {
                        "event_index": 3,
                        "event_type": "replacement",
                        "recipient_genome": "GenomeB",
                        "recipient_group": "gamma",
                        "source_relation": "cross_group",
                    },
                ],
            )

        rows = summary.set_index("recipient_genome")
        self.assertEqual(int(rows.loc["GenomeA", "contaminant_markers_added"]), 2)
        self.assertEqual(int(rows.loc["GenomeA", "duplicate_events"]), 1)
        self.assertEqual(int(rows.loc["GenomeA", "replacement_events"]), 1)
        self.assertEqual(int(rows.loc["GenomeA", "within_group_events"]), 1)
        self.assertEqual(int(rows.loc["GenomeA", "cross_group_events"]), 1)
        self.assertEqual(int(rows.loc["GenomeB", "replacement_events"]), 1)

    def test_evaluate_benchmark_run_reports_missing_taxa(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_dir = root / "benchmark"
            scenario_dir = benchmark_dir / "scenarios" / "replacement_only"
            run_dir = root / "run"
            scenario_dir.mkdir(parents=True)
            run_dir.mkdir()
            (run_dir / "aligned_final").mkdir()

            reference_tree = scenario_dir / "reference_tree.nwk"
            reference_tree.write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree_final.nwk").write_text("(GenomeA,GenomeB);\n")
            (run_dir / "marker_selection_rf_values.txt").write_text(
                "ProteinID MarkerGene RFdistance Status\n"
            )
            (run_dir / "aligned_final" / "MarkerX.faa").write_text(
                ">GenomeB|native_marker\nMPEPTIDE\n"
            )
            (scenario_dir / "events.tsv").write_text(
                "\t".join(
                    [
                        "event_index",
                        "scenario",
                        "event_type",
                        "recipient_genome",
                        "recipient_group",
                        "marker",
                        "native_record_id",
                        "donor_genome",
                        "donor_group",
                        "source_relation",
                        "donor_record_id",
                        "contaminant_record_id",
                        "expected_replacement_outcome",
                        "native_degrade_fraction",
                    ]
                )
                + "\n"
                + "\t".join(
                    [
                        "1",
                        "replacement_only",
                        "replacement",
                        "GenomeA",
                        "flavo",
                        "MarkerX",
                        "GenomeA|native_marker",
                        "GenomeB",
                        "gamma",
                        "cross_group",
                        "GenomeB|native_marker",
                        "GenomeA|contam__MarkerX__GenomeB__e001",
                        "DropMarkerOrRemoveContaminant",
                        "0.12",
                    ]
                )
                + "\n"
            )
            (benchmark_dir / "benchmark_manifest.json").write_text(
                json.dumps(
                    {
                        "scenarios": [
                            {
                                "name": "replacement_only",
                                "reference_tree_path": str(reference_tree),
                            }
                        ]
                    }
                )
            )

            result = benchmark.evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name="replacement_only",
                run_dir=run_dir,
                runtime_seconds=1.0,
            )

        self.assertFalse(result["final_taxa_match_reference"])
        self.assertEqual(result["final_missing_taxa_count"], 1)
        self.assertEqual(result["final_missing_taxa"], "GenomeC")
        self.assertEqual(result["collateral_genome_loss_count"], 1)
        self.assertEqual(result["collateral_genomes_lost"], "GenomeC")

    def test_evaluate_benchmark_run_uses_manifest_reference_taxa_over_pruned_reference_tree(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_dir = root / "benchmark"
            scenario_dir = benchmark_dir / "scenarios" / "combined"
            run_dir = root / "run"
            scenario_dir.mkdir(parents=True)
            run_dir.mkdir()
            (run_dir / "aligned_final").mkdir()

            reference_tree = scenario_dir / "reference_tree.nwk"
            reference_tree.write_text("(GenomeA,GenomeB);\n")
            (run_dir / "tree.nwk").write_text("(GenomeA,GenomeB,GenomeC);\n")
            (run_dir / "tree_final.nwk").write_text("(GenomeA,GenomeB,GenomeC);\n")
            (run_dir / "marker_selection_rf_values.txt").write_text(
                "ProteinID MarkerGene RFdistance Status\n"
            )
            (scenario_dir / "events.tsv").write_text(
                "\t".join(
                    [
                        "event_index",
                        "scenario",
                        "event_type",
                        "recipient_genome",
                        "recipient_group",
                        "marker",
                        "native_record_id",
                        "donor_genome",
                        "donor_group",
                        "source_relation",
                        "donor_record_id",
                        "contaminant_record_id",
                        "expected_replacement_outcome",
                        "native_degrade_fraction",
                    ]
                )
                + "\n"
            )
            (benchmark_dir / "benchmark_manifest.json").write_text(
                json.dumps(
                    {
                        "selected_genomes": ["GenomeA", "GenomeB", "GenomeC"],
                        "scenarios": [
                            {
                                "name": "combined",
                                "reference_tree_path": str(reference_tree),
                                "reference_taxa": ["GenomeA", "GenomeB", "GenomeC"],
                            }
                        ],
                    }
                )
            )

            result = benchmark.evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name="combined",
                run_dir=run_dir,
                runtime_seconds=1.0,
            )

        self.assertTrue(result["final_taxa_match_reference"])
        self.assertEqual(result["final_extra_taxa_count"], 0)
        self.assertEqual(result["final_missing_taxa_count"], 0)


if __name__ == "__main__":
    unittest.main()
