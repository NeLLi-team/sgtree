import json
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sgtree import benchmark


class BenchmarkTests(unittest.TestCase):
    def test_default_cleanup_profiles_use_gcp_for_singleton_scenarios(self):
        profiles = benchmark.DEFAULT_CLEANUP_PROFILES

        self.assertFalse(profiles["duplicate_only"]["singles"])
        self.assertEqual(profiles["duplicate_only"]["singles_mode"], "delta_rf")

        for scenario in ["replacement_only", "combined", "mixed_high_level"]:
            with self.subTest(scenario=scenario):
                profile = profiles[scenario]
                self.assertTrue(profile["singles"])
                self.assertEqual(profile["singles_mode"], "gcp")
                self.assertIn("gcp", profile["name"])

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

    def test_replacement_outcome_uses_exact_opaque_record_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            aligned_dir = run_dir / "aligned_final"
            aligned_dir.mkdir()
            aligned_path = aligned_dir / "MarkerX.faa"
            contaminant_id = "GenomeA|c_a13f|g_b27e"
            native_id = "GenomeA|c_c91d|g_d48a"

            aligned_path.write_text(
                f">{contaminant_id}\nMPEPTIDE\n>GenomeB|c_e22c|g_f63b\nMOTHER\n"
            )
            self.assertEqual(
                benchmark._replacement_outcome(
                    run_dir,
                    "GenomeA",
                    "MarkerX",
                    contaminant_id,
                    native_id,
                ),
                "contaminant_retained",
            )

            aligned_path.write_text(
                f">{native_id}\nMPEPTIDE\n>GenomeB|c_e22c|g_f63b\nMOTHER\n"
            )
            self.assertEqual(
                benchmark._replacement_outcome(
                    run_dir,
                    "GenomeA",
                    "MarkerX",
                    contaminant_id,
                    native_id,
                ),
                "native_retained",
            )

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

    def test_load_taxonomy_rows_handles_pandas_string_dtype_accessions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db_path = root / "taxonomy.duckdb"
            import duckdb
            import pandas as pd

            with duckdb.connect(str(db_path)) as con:
                con.execute(
                    """
                    create table gtdb_genomes (
                        assembly_accession varchar,
                        organism_name varchar,
                        phylum varchar,
                        class varchar,
                        order_name varchar,
                        family varchar,
                        genus varchar,
                        species varchar
                    )
                    """
                )
                con.execute(
                    """
                    create table non_gtdb_genomes (
                        assembly_accession varchar,
                        organism_name varchar,
                        phylum varchar,
                        class varchar,
                        order_name varchar,
                        family varchar,
                        genus varchar,
                        species varchar
                    )
                    """
                )
                con.execute(
                    """
                    insert into gtdb_genomes values
                    ('GCA_000009265.1', 'Genome A', 'P1', 'C1', 'O1', 'F1', 'G1', 'S1')
                    """
                )

            accessions = pd.Series(["GCA_000009265.1"], dtype="string").tolist()
            rows = benchmark._load_taxonomy_rows(accessions, db_path)

        self.assertEqual(rows.loc[0, "assembly_accession"], "GCA_000009265.1")
        self.assertEqual(rows.loc[0, "taxonomy_source"], "gtdb")

    def test_taxonomy_scope_matches_enforces_expected_rank_rules(self):
        recipient = {
            "domain": "Bacteria",
            "phylum": "Pseudomonadota",
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Enterobacteriaceae",
            "genus": "Escherichia",
            "species": "Escherichia coli",
        }
        same_family_other_genus = {
            "domain": "Bacteria",
            "phylum": "Pseudomonadota",
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Enterobacteriaceae",
            "genus": "Salmonella",
            "species": "Salmonella enterica",
        }
        same_order_other_family = {
            "domain": "Bacteria",
            "phylum": "Pseudomonadota",
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Vibrionaceae",
            "genus": "Vibrio",
            "species": "Vibrio cholerae",
        }
        same_class_other_order = {
            "domain": "Bacteria",
            "phylum": "Pseudomonadota",
            "class": "Gammaproteobacteria",
            "order_name": "Pseudomonadales",
            "family": "Pseudomonadaceae",
            "genus": "Pseudomonas",
            "species": "Pseudomonas aeruginosa",
        }
        same_phylum_other_class = {
            "domain": "Bacteria",
            "phylum": "Pseudomonadota",
            "class": "Alphaproteobacteria",
            "order_name": "Rhizobiales",
            "family": "Rhizobiaceae",
            "genus": "Rhizobium",
            "species": "Rhizobium leguminosarum",
        }
        same_domain_other_phylum = {
            "domain": "Bacteria",
            "phylum": "Bacteroidota",
            "class": "Bacteroidia",
            "order_name": "Flavobacteriales",
            "family": "Flavobacteriaceae",
            "genus": "Flavobacterium",
            "species": "Flavobacterium johnsoniae",
        }

        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_family_other_genus, "genus"))
        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_order_other_family, "family"))
        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_class_other_order, "order"))
        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_phylum_other_class, "class"))
        self.assertTrue(benchmark._taxonomy_scope_matches(recipient, same_domain_other_phylum, "phylum"))
        self.assertFalse(benchmark._taxonomy_scope_matches(recipient, same_order_other_family, "genus"))
        self.assertFalse(benchmark._taxonomy_scope_matches(recipient, same_class_other_order, "family"))
        self.assertFalse(benchmark._taxonomy_scope_matches(recipient, same_domain_other_phylum, "class"))

    def test_classify_taxonomic_distance_reports_source_level(self):
        recipient = {
            "phylum": "Pseudomonadota",
            "class": "Gammaproteobacteria",
            "order_name": "Enterobacterales",
            "family": "Enterobacteriaceae",
            "genus": "Escherichia",
            "species": "Escherichia coli",
        }
        cases = [
            (
                {"phylum": "Pseudomonadota", "class": "Gammaproteobacteria", "order_name": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Escherichia", "species": "Escherichia fergusonii"},
                "species",
                "different_species_same_genus",
            ),
            (
                {"phylum": "Pseudomonadota", "class": "Gammaproteobacteria", "order_name": "Enterobacterales", "family": "Enterobacteriaceae", "genus": "Salmonella", "species": "Salmonella enterica"},
                "genus",
                "different_genus_same_family",
            ),
            (
                {"phylum": "Pseudomonadota", "class": "Gammaproteobacteria", "order_name": "Enterobacterales", "family": "Vibrionaceae", "genus": "Vibrio", "species": "Vibrio cholerae"},
                "family",
                "different_family_same_order",
            ),
            (
                {"phylum": "Pseudomonadota", "class": "Gammaproteobacteria", "order_name": "Pseudomonadales", "family": "Pseudomonadaceae", "genus": "Pseudomonas", "species": "Pseudomonas aeruginosa"},
                "order",
                "different_order_same_class",
            ),
            (
                {"phylum": "Pseudomonadota", "class": "Alphaproteobacteria", "order_name": "Rhizobiales", "family": "Rhizobiaceae", "genus": "Rhizobium", "species": "Rhizobium leguminosarum"},
                "class",
                "different_class_same_phylum",
            ),
            (
                {"phylum": "Bacteroidota", "class": "Bacteroidia", "order_name": "Flavobacteriales", "family": "Flavobacteriaceae", "genus": "Flavobacterium", "species": "Flavobacterium johnsoniae"},
                "phylum",
                "different_phylum",
            ),
        ]

        for donor, expected_level, expected_label in cases:
            with self.subTest(expected_level=expected_level):
                result = benchmark.classify_taxonomic_distance(recipient, donor)
                self.assertEqual(result["contamination_source_taxonomic_level"], expected_level)
                self.assertEqual(result["contamination_source_taxonomic_label"], expected_label)
                self.assertEqual(result[f"contamination_source_differs_{expected_level}"], "yes")

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

    def test_evaluate_benchmark_run_treats_missing_and_empty_alignments_as_unknown(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_dir = root / "benchmark"
            scenario_dir = benchmark_dir / "scenarios" / "replacement_only"
            run_dir = root / "run"
            scenario_dir.mkdir(parents=True)
            (run_dir / "aligned_final").mkdir(parents=True)

            reference_tree = scenario_dir / "reference_tree.nwk"
            reference_tree.write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree_final.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "marker_selection_rf_values.txt").write_text(
                "ProteinID MarkerGene RFdistance Status\n"
            )
            (run_dir / "aligned_final" / "MarkerEmpty.faa").write_text("")
            header = [
                "event_index",
                "scenario",
                "event_type",
                "recipient_genome",
                "marker",
                "native_record_id",
                "contaminant_record_id",
            ]
            rows = [
                [
                    "1",
                    "replacement_only",
                    "replacement",
                    "GenomeA",
                    "MarkerMissing",
                    "GenomeA|c_102a|g_203b",
                    "GenomeA|c_304c|g_405d",
                ],
                [
                    "2",
                    "replacement_only",
                    "replacement",
                    "GenomeB",
                    "MarkerEmpty",
                    "GenomeB|c_506e|g_607f",
                    "GenomeB|c_708a|g_809b",
                ],
            ]
            (scenario_dir / "events.tsv").write_text(
                "\t".join(header) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n"
            )
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

            result = benchmark.evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name="replacement_only",
                run_dir=run_dir,
                runtime_seconds=1.0,
            )

        self.assertEqual(result["replacement_unknown"], 2)
        self.assertEqual(result["replacement_contaminant_removed"], 0)
        self.assertEqual(result["replacement_marker_dropped"], 0)
        self.assertEqual(result["contaminant_markers_removed"], 0)
        self.assertEqual(result["contaminant_markers_removed_fraction"], 0.0)

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
        self.assertEqual(result["replacement_marker_dropped"], 1)
        self.assertEqual(result["replacement_contaminant_removed"], 1)

    def test_evaluate_benchmark_run_does_not_credit_marker_drop_when_recipient_is_lost(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            benchmark_dir = root / "benchmark"
            scenario_dir = benchmark_dir / "scenarios" / "replacement_only"
            run_dir = root / "run"
            scenario_dir.mkdir(parents=True)
            (run_dir / "aligned_final").mkdir(parents=True)

            reference_tree = scenario_dir / "reference_tree.nwk"
            reference_tree.write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree_final.nwk").write_text("(GenomeB,GenomeC);\n")
            (run_dir / "marker_selection_rf_values.txt").write_text(
                "ProteinID MarkerGene RFdistance Status\n"
            )
            (run_dir / "aligned_final" / "MarkerX.faa").write_text(
                ">GenomeB|c_910c|g_011d\nMPEPTIDE\n"
            )
            (scenario_dir / "events.tsv").write_text(
                "event_index\tscenario\tevent_type\trecipient_genome\tmarker\t"
                "native_record_id\tcontaminant_record_id\n"
                "1\treplacement_only\treplacement\tGenomeA\tMarkerX\t"
                "GenomeA|c_122e|g_233f\tGenomeA|c_344a|g_455b\n"
            )
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

            result = benchmark.evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name="replacement_only",
                run_dir=run_dir,
                runtime_seconds=1.0,
            )

        self.assertEqual(result["replacement_recipient_lost"], 1)
        self.assertEqual(result["replacement_recipient_genome_loss_count"], 1)
        self.assertEqual(result["replacement_marker_dropped"], 0)
        self.assertEqual(result["replacement_contaminant_removed"], 0)
        self.assertEqual(result["contaminant_markers_removed"], 0)

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

    def test_evaluate_benchmark_run_resolves_manifest_paths_relative_to_runs_root(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            backup_root = root / "backup_root"
            benchmark_dir = backup_root / "runs" / "benchmarks" / "panel_a"
            scenario_dir = benchmark_dir / "scenarios" / "replacement_only"
            run_dir = root / "run"
            scenario_dir.mkdir(parents=True)
            run_dir.mkdir()
            (run_dir / "aligned_final").mkdir()

            relative_reference_tree = Path("runs/benchmarks/panel_a/scenarios/replacement_only/reference_tree.nwk")
            reference_tree = backup_root / relative_reference_tree
            reference_tree.write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "tree_final.nwk").write_text("((GenomeA,GenomeB),GenomeC);\n")
            (run_dir / "marker_selection_rf_values.txt").write_text("ProteinID MarkerGene RFdistance Status\n")
            (run_dir / "aligned_final" / "MarkerX.faa").write_text("")
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
                                "name": "replacement_only",
                                "reference_tree_path": str(relative_reference_tree),
                                "reference_taxa": ["GenomeA", "GenomeB", "GenomeC"],
                            }
                        ],
                    }
                )
            )

            result = benchmark.evaluate_benchmark_run(
                benchmark_dir=benchmark_dir,
                scenario_name="replacement_only",
                run_dir=run_dir,
                runtime_seconds=1.0,
            )

        self.assertEqual(result["truth_tree_path"], str(reference_tree))
        self.assertTrue(result["final_taxa_match_reference"])

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


if __name__ == "__main__":
    unittest.main()
