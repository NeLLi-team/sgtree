import tempfile
import unittest
from pathlib import Path

from sgtree.marker_selection import _write_rf_values_file


def _records_bundle() -> list[list[dict]]:
    return [
        [
            {"protein_id": "gB/p1", "marker": "M2", "rf_distance": 0.31234, "status": "Kept"},
            {"protein_id": "gA/p2", "marker": "M1", "rf_distance": 0.05, "status": "Removed"},
        ],
        [],
        [
            {"protein_id": "gA/p1", "marker": "M1", "rf_distance": 0.1, "status": "Kept"},
        ],
    ]


class RfValuesDeterministicTests(unittest.TestCase):
    def test_output_is_byte_identical_across_worker_orderings(self):
        bundle = _records_bundle()
        reversed_bundle = list(reversed([list(rs) for rs in bundle]))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            p1 = tmp / "run1.txt"
            p2 = tmp / "run2.txt"
            _write_rf_values_file(str(p1), bundle)
            _write_rf_values_file(str(p2), reversed_bundle)
            self.assertEqual(p1.read_bytes(), p2.read_bytes())

    def test_lines_are_sorted_by_marker_then_protein(self):
        bundle = _records_bundle()
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "rf.txt"
            _write_rf_values_file(str(out), bundle)
            lines = out.read_text().splitlines()
            self.assertEqual(lines[0], "ProteinID MarkerGene RFdistance Status")
            data = lines[1:]
            keys = [(parts[1], parts[0]) for parts in (l.split() for l in data)]
            self.assertEqual(keys, sorted(keys))

    def test_handles_empty_and_none_worker_results(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "rf.txt"
            _write_rf_values_file(str(out), None)
            self.assertEqual(out.read_text(), "ProteinID MarkerGene RFdistance Status\n")

            _write_rf_values_file(str(out), [[]])
            self.assertEqual(out.read_text(), "ProteinID MarkerGene RFdistance Status\n")


if __name__ == "__main__":
    unittest.main()
