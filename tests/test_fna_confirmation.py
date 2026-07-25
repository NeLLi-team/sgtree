import hashlib
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sgtree.benchmarks import fna_confirmation


class FnaConfirmationTests(unittest.TestCase):
    def test_fragment_bounds_stays_inside_selected_marker_barriers(self):
        rows = []
        for index in range(10):
            rows.append(
                {
                    "genome_id": "donor",
                    "contig_id": "contig",
                    "normalized_header": f"donor|contig|gene_{index}",
                    "begin": index * 100 + 1,
                    "end": index * 100 + 90,
                }
            )
        gene_calls = pd.DataFrame(rows)
        target = "donor|contig|gene_5"
        barriers = {
            "donor|contig|gene_1",
            target,
            "donor|contig|gene_9",
        }

        result = fna_confirmation._fragment_bounds(
            gene_calls,
            target,
            barriers,
        )

        self.assertIsNotNone(result)
        self.assertEqual(result["planned_non_marker_gene_count"], 6)
        self.assertNotIn("donor|contig|gene_1", result["planned_gene_headers"])
        self.assertNotIn("donor|contig|gene_9", result["planned_gene_headers"])
        self.assertIn(target, result["planned_gene_headers"])

    def test_gene_rich_event_masks_native_cds_and_appends_donor_fragment(self):
        recipient = [SeqRecord(Seq("CCCATGAAATAAGGG"), id="recipient_contig")]
        donor = [SeqRecord(Seq("TTTATGAAATAACCC"), id="donor_contig")]
        protein_hash = hashlib.sha256(b"MK").hexdigest()
        event = {
            "event_id": "panel_near",
            "role": "near",
            "native_contig_id": "recipient_contig",
            "native_begin": 4,
            "native_end": 12,
            "donor_contig_id": "donor_contig",
            "donor_begin": 4,
            "donor_end": 12,
            "donor_strand": 1,
            "donor_translation_table": 11,
            "donor_protein_sha256": protein_hash,
            "fragment_start0": 0,
            "fragment_end0": 15,
            "foreign_contig_id": "foreign_contig",
        }

        fna_confirmation._apply_event_to_records(recipient, donor, event)

        self.assertEqual(str(recipient[0].seq), "CCCNNNNNNNNNGGG")
        self.assertEqual(recipient[1].id, "foreign_contig")
        self.assertEqual(str(recipient[1].seq), str(donor[0].seq))

    def test_manifest_validator_enforces_frozen_six_panel_shape(self):
        panels = []
        for lineage in fna_confirmation.LINEAGES:
            for panel_index, seed in enumerate(
                fna_confirmation.PANEL_SEEDS,
                start=1,
            ):
                panel_id = f"{lineage}_p{panel_index}_seed{seed}"
                events = []
                for event_index, role in enumerate(
                    fna_confirmation.EVENT_ORDER
                ):
                    event_class = (
                        "gene_rich_replacement"
                        if role in {"near", "intermediate", "far"}
                        else "solo_marker_control"
                        if role == "solo"
                        else "native_contig_sentinel"
                    )
                    events.append(
                        {
                            "event_id": f"{panel_id}_{role}",
                            "role": role,
                            "event_class": event_class,
                            "recipient_genome": f"{panel_id}_r{event_index}",
                            "donor_genome": f"{panel_id}_d{event_index}",
                            "marker": f"m{event_index}",
                            "donor_cds_terminal_stop": True,
                            "donor_terminal_stop_codon": "TAA",
                            "native_protein_sha256": f"native{event_index}",
                            "donor_protein_sha256": f"donor{event_index}",
                            "donor_protein_differs_from_native": True,
                        }
                    )
                panels.append(
                    {
                        "panel_id": panel_id,
                        "lineage": lineage,
                        "genomes": [
                            f"{panel_id}_g{index}"
                            for index in range(
                                fna_confirmation.PANEL_GENOME_COUNT
                            )
                        ],
                        "gunc_clean_control_genome": f"{panel_id}_g0",
                        "markers": [
                            f"m{index}"
                            for index in range(
                                fna_confirmation.MARKER_COUNT
                            )
                        ],
                        "events": events,
                        "contexts": fna_confirmation._panel_contexts(
                            panel_id,
                            events,
                        ),
                    }
                )

        fna_confirmation._validate_frozen_manifest({"panels": panels})
        panels[0]["events"][0]["donor_protein_sha256"] = panels[0]["events"][0][
            "native_protein_sha256"
        ]
        with self.assertRaisesRegex(ValueError, "sequence-identical"):
            fna_confirmation._validate_frozen_manifest({"panels": panels})

    def test_native_sentinel_may_share_its_native_contig_with_other_markers(self):
        self.assertTrue(
            fna_confirmation._target_contig_marker_layout_pass(
                {"role": "sentinel", "marker": "m2"},
                ["m1", "m2", "m3"],
            )
        )
        self.assertFalse(
            fna_confirmation._target_contig_marker_layout_pass(
                {"role": "near", "marker": "m2"},
                ["m1", "m2"],
            )
        )

    def test_partial_donor_gene_is_rejected_before_event_freeze(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "donor.fna"
            source.write_text(">contig\nATGAAA\n")
            data = {
                "source_files": {"donor": source},
            }
            gene = {
                "normalized_header": "donor|contig|gene_1",
                "genome_id": "donor",
                "contig_id": "contig",
                "begin": 1,
                "end": 6,
                "strand": 1,
            }
            self.assertFalse(
                fna_confirmation._donor_gene_has_terminal_stop(data, gene)
            )

    def test_solo_donor_requires_stop_valid_under_codes_4_and_11(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "donor.fna"
            source.write_text(">contig\nATGAAATGA\n")
            data = {
                "source_files": {"donor": source},
            }
            gene = {
                "normalized_header": "donor|contig|gene_1",
                "genome_id": "donor",
                "contig_id": "contig",
                "begin": 1,
                "end": 9,
                "strand": 1,
            }

            self.assertTrue(
                fna_confirmation._donor_gene_has_terminal_stop(data, gene)
            )
            self.assertEqual(
                fna_confirmation._donor_gene_terminal_stop_codon(data, gene),
                "TGA",
            )
            self.assertNotIn(
                fna_confirmation._donor_gene_terminal_stop_codon(data, gene),
                {"TAA", "TAG"},
            )

    def test_sentinel_donor_requires_cds_translation_to_match_staged_protein(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "donor.fna"
            source.write_text(">contig\nTTGAAATAA\n")
            data = {
                "source_files": {"donor": source},
            }
            gene = {
                "normalized_header": "donor|contig|gene_1",
                "genome_id": "donor",
                "contig_id": "contig",
                "begin": 1,
                "end": 9,
                "strand": 1,
                "translation_table": 11,
            }

            self.assertTrue(
                fna_confirmation._donor_cds_matches_staged_protein(
                    data,
                    gene,
                    "LK",
                )
            )
            self.assertFalse(
                fna_confirmation._donor_cds_matches_staged_protein(
                    data,
                    gene,
                    "MK",
                )
            )

    def test_frozen_inference_input_checksum_is_enforced(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "markers.hmm"
            source = Path(tmpdir) / "genome.fna"
            model.write_bytes(b"model")
            source.write_bytes(b">contig\nATG\n")
            manifest = {
                "models_path": str(model),
                "models_sha256": hashlib.sha256(b"model").hexdigest(),
                "models_bytes": model.stat().st_size,
                "source_fna_artifacts": [
                    {
                        "path": str(source),
                        "sha256": hashlib.sha256(b">contig\nATG\n").hexdigest(),
                        "bytes": source.stat().st_size,
                    }
                ],
            }

            fna_confirmation._verify_inference_inputs(manifest)
            source.write_bytes(b">contig\nGTG\n")
            with self.assertRaisesRegex(ValueError, "checksum changed"):
                fna_confirmation._verify_inference_inputs(manifest)

    def test_cluster_interval_uses_panel_counts_not_context_rows(self):
        frame = pd.DataFrame(
            [
                {
                    "panel_id": "alpha_p1_seed1103",
                    "numerator": 1,
                    "denominator": 1,
                },
                {
                    "panel_id": "alpha_p1_seed1103",
                    "numerator": 0,
                    "denominator": 1,
                },
                {
                    "panel_id": "alpha_p2_seed1301",
                    "numerator": 1,
                    "denominator": 1,
                },
            ]
        )

        result = fna_confirmation._cluster_ratio(
            frame,
            numerator_column="numerator",
            denominator_column="denominator",
        )

        self.assertEqual(result["numerator"], 2)
        self.assertEqual(result["denominator"], 3)
        self.assertAlmostEqual(result["estimate"], 2 / 3)
        self.assertEqual(result["independent_unit"], "base_panel")

    def test_source_archive_is_byte_reproducible(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "source.py"
            source.write_text("value = 1\n")
            first = root / "first.tar.gz"
            second = root / "second.tar.gz"
            members = [(source, "src/source.py")]

            fna_confirmation._write_deterministic_tar_gz(first, members)
            fna_confirmation._write_deterministic_tar_gz(second, members)

            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(first.read_bytes()[4:8], b"\0\0\0\0")


if __name__ == "__main__":
    unittest.main()
