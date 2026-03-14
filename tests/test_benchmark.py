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

        self.assertTrue(record.id.startswith("RecipientB|contig__contam__MarkerX__DonorA__e001|contam__MarkerX__DonorA__e001"))
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
            id="RecipientB|contig__contam__MarkerX__DonorA__e001|contam__MarkerX__DonorA__e001",
            description="RecipientB|contig__contam__MarkerX__DonorA__e001|contam__MarkerX__DonorA__e001",
        )

        updated = benchmark.apply_replacement_event(
            recipient_records,
            native_record_id="RecipientB|native_marker",
            contaminant_record=contaminant,
        )

        self.assertNotIn("RecipientB|native_marker", updated)
        self.assertIn("RecipientB|background", updated)
        self.assertIn("RecipientB|contig__contam__MarkerX__DonorA__e001|contam__MarkerX__DonorA__e001", updated)

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

    def test_choose_markers_for_pair_requires_marker_presence_in_both_genomes(self):
        native_map = {
            "GenomeA": {"Marker1": "GenomeA|m1", "Marker2": "GenomeA|m2"},
            "GenomeB": {"Marker1": "GenomeB|m1"},
        }

        chosen = benchmark._choose_markers_for_pair(
            used_pairs=set(),
            genomes=("GenomeA", "GenomeB"),
            markers=["Marker1", "Marker2"],
            native_map=native_map,
            n_needed=2,
            rng=benchmark.Random(42),
        )

        self.assertEqual(chosen, ["Marker1"])

    def test_donor_candidates_with_marker_filters_missing_marker_genomes(self):
        native_map = {
            "GenomeA": {"Marker1": "GenomeA|m1"},
            "GenomeB": {"Marker2": "GenomeB|m2"},
            "GenomeC": {"Marker1": "GenomeC|m1"},
        }

        donors = benchmark._donor_candidates_with_marker(
            ["GenomeA", "GenomeB", "GenomeC"],
            "Marker1",
            native_map,
        )

        self.assertEqual(donors, ["GenomeA", "GenomeC"])

    def test_normalize_assembly_accession_parses_supported_filename_forms(self):
        self.assertEqual(
            benchmark._normalize_assembly_accession("FLAV__GCA_000016645-1"),
            "GCA_000016645.1",
        )
        self.assertEqual(
            benchmark._normalize_assembly_accession("GAMMA__GCA-000147015-1"),
            "GCA_000147015.1",
        )
        self.assertEqual(
            benchmark._normalize_assembly_accession("CHLAMYDIOTA__GCF_123456789-3"),
            "GCF_123456789.3",
        )

    def test_taxonomy_scope_matches_enforces_expected_rank_rules(self):
        recipient = {
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Enterobacteriaceae",
            "genus": "Escherichia",
        }
        same_family_other_genus = {
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Enterobacteriaceae",
            "genus": "Salmonella",
        }
        same_order_other_family = {
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Vibrionaceae",
            "genus": "Vibrio",
        }
        same_class_other_order = {
            "class": "Gammaproteobacteria",
            "order_name": "Pseudomonadales",
            "family": "Pseudomonadaceae",
            "genus": "Pseudomonas",
        }

        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_family_other_genus, "genus"))
        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_order_other_family, "family"))
        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_class_other_order, "order"))
        self.assertFalse(benchmark._taxonomy_scope_matches(recipient, same_order_other_family, "genus"))
        self.assertFalse(benchmark._taxonomy_scope_matches(recipient, same_class_other_order, "family"))

    def test_taxonomic_donor_candidates_filter_by_scope_and_marker(self):
        recipient_taxonomy = {
            "RecipientA": {
                "class": "Gammaproteobacteria",
                "order_name": "Enterobacterales",
                "family": "Enterobacteriaceae",
                "genus": "Escherichia",
            }
        }
        donor_taxonomy = {
            "DonorGood": {
                "class": "Gammaproteobacteria",
                "order_name": "Enterobacterales",
                "family": "Enterobacteriaceae",
                "genus": "Salmonella",
            },
            "DonorWrongFamily": {
                "class": "Gammaproteobacteria",
                "order_name": "Enterobacterales",
                "family": "Vibrionaceae",
                "genus": "Vibrio",
            },
            "DonorMissingMarker": {
                "class": "Gammaproteobacteria",
                "order_name": "Enterobacterales",
                "family": "Enterobacteriaceae",
                "genus": "Klebsiella",
            },
        }
        donor_native_map = {
            "DonorGood": {"Marker1": "DonorGood|m1"},
            "DonorWrongFamily": {"Marker1": "DonorWrongFamily|m1"},
            "DonorMissingMarker": {"Marker2": "DonorMissingMarker|m2"},
        }

        donors = benchmark._taxonomic_donor_candidates(
            recipient_genome="RecipientA",
            marker="Marker1",
            scope="genus",
            recipient_taxonomy=recipient_taxonomy,
            donor_taxonomy=donor_taxonomy,
            donor_native_map=donor_native_map,
            truth_tree=None,
        )

        self.assertEqual(donors, ["DonorGood"])

    def test_choose_recipient_sets_honors_overlap_layout(self):
        duplicate_recipients, replacement_recipients = benchmark._choose_recipient_sets(
            [f"Genome{idx}" for idx in range(10)],
            duplicate_recipients=5,
            replacement_recipients=5,
            overlap_recipients=2,
            rng=benchmark.Random(42),
        )

        overlap = set(duplicate_recipients) & set(replacement_recipients)
        self.assertEqual(len(duplicate_recipients), 5)
        self.assertEqual(len(replacement_recipients), 5)
        self.assertEqual(len(overlap), 2)
        self.assertEqual(len(set(duplicate_recipients) | set(replacement_recipients)), 8)

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

    def test_prepare_source_subset_accepts_fna_sources(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_dir = root / "source"
            outdir = root / "subset"
            source_dir.mkdir()
            (source_dir / "GenomeA.fna").write_text(">contig1\nATGAAATTTAAATAG\n")
            (source_dir / "GenomeB.fna").write_text(">contig1\nATGAAATTTAAATAG\n")

            benchmark.prepare_source_subset(
                source_dir=source_dir,
                outdir=outdir,
                list_file=None,
                n_candidates=1,
                seed=42,
            )

            files = sorted(outdir.iterdir())
            self.assertEqual(len(files), 1)
            self.assertEqual(files[0].suffix, ".fna")

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

    def test_evaluate_benchmark_run_reports_singleton_collateral_removals(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_dir = root / "benchmark"
            scenario_dir = benchmark_dir / "scenarios" / "replacement_only"
            run_dir = root / "run"
            (scenario_dir).mkdir(parents=True)
            (run_dir / "aligned_final").mkdir(parents=True)
            (run_dir / "protTrees" / "no_duplicates" / "out").mkdir(parents=True)
            (run_dir / "protTrees" / "no_singles").mkdir(parents=True)

            reference_tree = scenario_dir / "reference_tree.nwk"
            reference_tree.write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree_final.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
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
            (run_dir / "aligned_final" / "MarkerX.faa").write_text("")
            (benchmark_dir / "benchmark_manifest.json").write_text(
                json.dumps(
                    {
                        "selected_genomes": ["GenomeA", "GenomeB", "GenomeC"],
                        "scenarios": [
                            {
                                "name": "replacement_only",
                                "reference_tree_path": str(reference_tree),
                                "reference_taxa": ["GenomeA", "GenomeB", "GenomeC"],
                            }
                        ],
                    }
                )
            )

            (run_dir / "protTrees" / "no_duplicates" / "out" / "_no_dups_MarkerX_.nw").write_text(
                "(GenomeA|x,GenomeB|x,GenomeC|x);\n"
            )
            (run_dir / "protTrees" / "no_singles" / "_no_dups_MarkerX_.nw").write_text(
                "(GenomeB|x,GenomeC|x);\n"
            )
            (run_dir / "protTrees" / "no_duplicates" / "out" / "_no_dups_MarkerY_.nw").write_text(
                "(GenomeA|y,GenomeB|y,GenomeC|y);\n"
            )
            (run_dir / "protTrees" / "no_singles" / "_no_dups_MarkerY_.nw").write_text(
                "(GenomeA|y,GenomeB|y);\n"
            )

            result = benchmark.evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name="replacement_only",
                run_dir=run_dir,
                runtime_seconds=1.0,
            )

        self.assertEqual(result["singleton_intended_removed_count"], 1)
        self.assertEqual(result["singleton_intended_removed"], "GenomeA:MarkerX")
        self.assertEqual(result["singleton_collateral_removed_count"], 1)
        self.assertEqual(result["singleton_collateral_removed"], "GenomeC:MarkerY")
        self.assertEqual(result["singleton_collateral_genome_count"], 1)
        self.assertEqual(result["singleton_collateral_genomes"], "GenomeC")

    def test_export_benchmark_tables_writes_alignment_comparison_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_dir = root / "benchmark"
            results_dir = benchmark_dir / "results"
            mafft_dir = benchmark_dir / "results_mafft_20260313"
            outdir = root / "docs_data"
            benchmark_dir.mkdir()
            results_dir.mkdir()
            mafft_dir.mkdir()

            (benchmark_dir / "benchmark_manifest.json").write_text(
                json.dumps(
                    {
                        "lineage_label": "flavo",
                        "taxonomic_scope": "genus",
                        "scenarios": [{"name": "duplicate_only"}],
                    }
                )
            )
            baseline = (
                "scenario\ttree_rf_norm\tcontaminant_markers_removed_fraction\truntime_seconds\t"
                "final_missing_taxa_count\tstatus\trun_dir\n"
                "duplicate_only\t0.0\t0.5\t10.0\t0\tok\tbaseline_run\n"
            )
            comparison = (
                "scenario\ttree_rf_norm\tcontaminant_markers_removed_fraction\truntime_seconds\t"
                "final_missing_taxa_count\tstatus\trun_dir\talignment_method\n"
                "duplicate_only\t0.1\t0.75\t12.0\t0\tok\tmafft_run\tmafft\n"
            )
            (results_dir / "summary.tsv").write_text(baseline)
            (mafft_dir / "summary.tsv").write_text(comparison)

            benchmark.export_benchmark_tables([benchmark_dir], outdir)

            comparison_df = benchmark.pd.read_csv(outdir / "benchmark_alignment_comparison.tsv", sep="\t")
            summary_df = benchmark.pd.read_csv(outdir / "benchmark_alignment_comparison_summary.tsv", sep="\t")

        self.assertEqual(len(comparison_df), 1)
        self.assertEqual(comparison_df.loc[0, "comparison_alignment_method"], "mafft")
        self.assertAlmostEqual(comparison_df.loc[0, "rf_delta_comparison_minus_baseline"], 0.1)
        self.assertEqual(comparison_df.loc[0, "rf_change_direction"], "worse")
        self.assertEqual(len(summary_df), 1)
        self.assertEqual(summary_df.loc[0, "worsened_count"], 1)
        self.assertAlmostEqual(summary_df.loc[0, "runtime_ratio_mean"], 1.2)


if __name__ == "__main__":
    unittest.main()
