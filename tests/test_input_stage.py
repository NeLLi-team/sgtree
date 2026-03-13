import sys
import tempfile
import types
import unittest
from pathlib import Path

import pandas as pd

from sgtree.config import Config
from sgtree.fasta_normalize import normalize_and_concat_proteomes
from sgtree.input_stage import detect_input_format, gene_call_inputs, write_genome_manifest
from sgtree.search import build_working_df


class _FakeGene:
    def __init__(self, protein: str, begin: int = 1, end: int = 12, strand: int = 1, translation_table: int = 11):
        self._protein = protein
        self.begin = begin
        self.end = end
        self.strand = strand
        self.translation_table = translation_table

    def translate(self, include_stop: bool = False, strict: bool = False):
        return self._protein


class _FakeGeneFinder:
    def __init__(self, meta: bool = True):
        self.meta = meta

    def find_genes(self, sequence: bytes):
        return [_FakeGene("MPEPTIDE"), _FakeGene("MSEQ", begin=30, end=38)]


class InputStageTests(unittest.TestCase):
    def test_detect_input_format_distinguishes_fna_and_faa(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            faa_dir = tmp / "faa"
            fna_dir = tmp / "fna"
            faa_dir.mkdir()
            fna_dir.mkdir()
            (faa_dir / "GenomeA.faa").write_text(">prot1\nMPEPTIDE\n")
            (fna_dir / "GenomeA.fna").write_text(">contig_1\nATGAAATTTAAATAG\n")

            self.assertEqual(detect_input_format(str(faa_dir)), "faa")
            self.assertEqual(detect_input_format(str(fna_dir)), "fna")

    def test_normalize_and_concat_proteomes_infers_contig_from_protein_suffix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            faa_dir = tmp / "faa"
            faa_dir.mkdir()
            (faa_dir / "GenomeA.faa").write_text(
                ">contig0001_42\nMPEPTIDE\n>contig0002_7\nMSEQ\n"
            )
            out_fasta = tmp / "proteomes"
            map_path = tmp / "map.tsv"

            stats = normalize_and_concat_proteomes(str(faa_dir), str(out_fasta), str(map_path))

            self.assertEqual(stats["genomes"], 1)
            self.assertEqual(stats["contigs"], 2)
            contents = out_fasta.read_text()
            self.assertIn(">GenomeA|contig0001|contig0001_42", contents)
            self.assertIn(">GenomeA|contig0002|contig0002_7", contents)
            map_text = map_path.read_text()
            self.assertIn("contig0001", map_text)
            self.assertIn("suffix", map_text)

    def test_gene_call_inputs_emits_genome_contig_gene_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fna_dir = tmp / "fna"
            out_dir = tmp / "called"
            fna_dir.mkdir()
            (fna_dir / "GenomeA.fna").write_text(">contigAlpha\nATGAAATTTAAATAG\n")

            fake_module = types.SimpleNamespace(GeneFinder=_FakeGeneFinder)
            original = sys.modules.get("pyrodigal")
            sys.modules["pyrodigal"] = fake_module
            try:
                stats = gene_call_inputs(str(fna_dir), str(out_dir), str(tmp / "gene_calls.tsv"))
            finally:
                if original is None:
                    del sys.modules["pyrodigal"]
                else:
                    sys.modules["pyrodigal"] = original

            self.assertEqual(stats.input_format, "fna")
            called = (out_dir / "GenomeA.faa").read_text()
            self.assertIn(">GenomeA|contigAlpha|gene_000001", called)
            self.assertIn(">GenomeA|contigAlpha|gene_000002", called)

    def test_write_genome_manifest_records_fna_paths_and_sizes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fna_dir = tmp / "fna"
            staged = tmp / "staged"
            fna_dir.mkdir()
            staged.mkdir()
            (fna_dir / "GenomeA.fna").write_text(">contigAlpha\nATGAAATTTAAATAG\n>contigBeta\nATGAAATTT\n")
            manifest = tmp / "manifest.tsv"

            write_genome_manifest(
                str(fna_dir),
                input_format="fna",
                manifest_path=str(manifest),
                staged_source=str(staged),
            )

            text = manifest.read_text()
            self.assertIn("GenomeA", text)
            self.assertIn("GenomeA.fna", text)
            self.assertIn("GenomeA.faa", text)
            self.assertIn("\t2\t24\n", text)

    def test_build_working_df_adds_contig_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cfg = Config(
                genomedir=str(tmp / "input"),
                modeldir=str(tmp / "models.hmm"),
                outdir=str(tmp / "run"),
                num_cpus=1,
                percent_models=0,
                input_format="faa",
                lflt_fraction=0.0,
                aln_method="hmmalign",
                tree_method="fasttree",
                iqtree_fast=True,
                iqtree_model="LG+F+I+G4",
                hmmsearch_cutoff="cut_ga",
                hmmsearch_evalue=1e-5,
                selection_mode="coordinate",
                selection_max_rounds=5,
                selection_global_rounds=1,
                lock_references=False,
                max_sdup=-1,
                max_dupl=-1.0,
                ref=None,
                ref_concat=str(tmp / "ref_cache"),
                marker_selection=False,
                singles=False,
                singles_mode="delta_rf",
                num_nei=0,
                singles_min_rfdist=0.25,
                keep_intermediates=True,
                is_ref=False,
                start_time="now",
            )
            Path(cfg.tables_dir).mkdir(parents=True, exist_ok=True)
            finaldf = pd.DataFrame(
                [
                    {
                        0: "GenomeA|contigA|gene_000001",
                        3: "MarkerX",
                        7: "120.0",
                    }
                ]
            )

            df, df_fordups = build_working_df(cfg, finaldf)

            self.assertEqual(df.iloc[0]["genome_id"], "GenomeA")
            self.assertEqual(df.iloc[0]["contig_id"], "contigA")
            self.assertEqual(df.iloc[0]["gene_id"], "gene_000001")
            self.assertIn("GenomeA/contigA/gene_000001", df_fordups.index)


if __name__ == "__main__":
    unittest.main()
