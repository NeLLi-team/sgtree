import io
import tempfile
import tomllib
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import sgtree
from sgtree._version import DISPLAY_VERSION
from sgtree import sgtree_logging
from sgtree.config import Config


def _make_config(tmp: Path) -> Config:
    genomes = tmp / "genomes"
    genomes.mkdir()
    (genomes / "GenomeA.faa").write_text(">proteinA\nMKT\n")
    model = tmp / "models.hmm"
    model.write_text("NAME  marker_one\nNAME  marker_two\n")
    outdir = tmp / "run"
    outdir.mkdir()
    return Config(
        genomedir=str(genomes),
        modeldir=str(model),
        outdir=str(outdir),
        num_cpus=1,
        percent_models=0,
        input_format="faa",
        lflt_fraction=0.0,
        aln_method="hmmalign",
        tree_method="veryfasttree",
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
        start_time="2026-04-25 12:34:56.000000",
    )


class VersionBannerTests(unittest.TestCase):
    def test_runtime_version_matches_pyproject_metadata(self):
        pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
        metadata = tomllib.loads(pyproject.read_text())
        self.assertEqual(sgtree.__version__, metadata["project"]["version"])

    def test_print_banner_uses_package_display_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = _make_config(Path(tmpdir))
            output = io.StringIO()

            with redirect_stdout(output):
                cfg.print_banner()

            text = output.getvalue()
            self.assertIn(DISPLAY_VERSION, text)
            self.assertNotIn("Sg_Tree v.2", text)

    def test_logfile_header_uses_package_display_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = _make_config(Path(tmpdir))

            sgtree_logging.write_logfile(cfg, {"example": ("timestamp", "0.0")})

            [logfile] = Path(cfg.outdir).glob("logfile_*.txt")
            text = logfile.read_text()
            self.assertIn(DISPLAY_VERSION, text)
            self.assertNotIn("Sg_Tree v.2", text)


if __name__ == "__main__":
    unittest.main()
