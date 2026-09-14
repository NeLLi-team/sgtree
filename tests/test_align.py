import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from sgtree.align import run_alignment


class HmmalignModelPathTests(unittest.TestCase):
    def test_dotted_marker_name_maps_to_complete_hmm_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            extracted = tmp / "extracted"
            aligned = tmp / "aligned"
            extracted.mkdir()
            (extracted / "marker.v1.faa").write_text(">GenomeA|c1|p1\nMKT\n")
            cfg = SimpleNamespace(
                extracted_seqs_dir=str(extracted),
                aligned_dir=str(aligned),
                outdir=str(tmp),
                models_path=str(tmp / "models.hmm"),
                num_cpus=1,
                aln_method="hmmalign",
            )

            with (
                patch("sgtree.align._split_models"),
                patch("sgtree.align.map_processed") as map_processed,
            ):
                run_alignment(cfg)

            payload = map_processed.call_args.args[1]
            self.assertEqual(len(payload), 1)
            self.assertEqual(
                payload[0][2],
                str(tmp / "models_split" / "marker.v1.hmm"),
            )


if __name__ == "__main__":
    unittest.main()
