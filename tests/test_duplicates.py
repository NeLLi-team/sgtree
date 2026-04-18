import tempfile
import unittest
from pathlib import Path

from sgtree.duplicates import _process_file_worker


class EliminateDuplicatesTests(unittest.TestCase):
    def test_keeps_highest_scored_id_per_genome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "M1.faa"
            out_dir = tmp / "out"
            out_dir.mkdir()

            fasta.write_text(
                ">g1|p1\nAAAA\n>g1|p2\nCCCC\n>g2|p3\nGGGG\n>g1|p4\nTTTT\n"
            )

            # score_lookup key uses "|" form; g1|p1 has best score
            score_lookup = {
                "g1|p1": "g1/p1:10.0",
                "g1|p2": "g1/p2:8.0",
                "g1|p4": "g1/p4:7.0",
                "g2|p3": "g2/p3:5.0",
            }

            _process_file_worker((str(fasta), str(out_dir), score_lookup))

            kept = (out_dir / "M1.faa").read_text()
            self.assertIn(">g1|p1", kept)
            self.assertNotIn(">g1|p2", kept)
            self.assertNotIn(">g1|p4", kept)
            self.assertIn(">g2|p3", kept)

    def test_many_hits_scales_linearly(self):
        # 20 genomes x 30 copies = 600 hits. On the old O(n^2) path this
        # runs ~360k comparisons; O(n) via Counter runs ~1200. We don't
        # time-assert here — the correctness check is that each genome
        # keeps exactly one (the one with the highest score).
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "M1.faa"
            out_dir = tmp / "out"
            out_dir.mkdir()

            headers: list[str] = []
            score_lookup: dict[str, str] = {}
            for g in range(20):
                for copy in range(30):
                    key = f"g{g}|p{copy}"
                    # Ascending score so copy 29 wins.
                    score_lookup[key] = f"g{g}/p{copy}:{float(copy)}"
                    headers.append(f">{key}\nAAAA")
            fasta.write_text("\n".join(headers) + "\n")

            _process_file_worker((str(fasta), str(out_dir), score_lookup))
            kept = (out_dir / "M1.faa").read_text().splitlines()
            kept_ids = [h[1:] for h in kept if h.startswith(">")]
            self.assertEqual(len(kept_ids), 20)
            self.assertEqual(
                sorted(kept_ids),
                sorted(f"g{g}|p29" for g in range(20)),
            )


if __name__ == "__main__":
    unittest.main()
