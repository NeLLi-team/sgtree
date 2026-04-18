import io
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from sgtree.extract import _concat_proteomes


def _write(path: Path, text: str) -> None:
    path.write_text(text)


def _records(path: Path) -> list[tuple[str, str]]:
    return [(r.id, str(r.seq)) for r in SeqIO.parse(str(path), "fasta")]


class ConcatProteomesTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.mkdtemp()
        self.tmp = Path(self._tmp)

    def tearDown(self) -> None:
        shutil.rmtree(self._tmp, ignore_errors=True)

    def test_concat_preserves_all_records_when_both_have_trailing_newline(self) -> None:
        a = self.tmp / "a.faa"
        b = self.tmp / "b.faa"
        out = self.tmp / "combined.faa"
        _write(a, ">g1|p1\nMKT\n>g1|p2\nAAA\n")
        _write(b, ">g2|p1\nCCC\n>g2|p2\nDDD\n")

        _concat_proteomes([str(a), str(b)], str(out))

        expected = set(_records(a)) | set(_records(b))
        self.assertEqual(set(_records(out)), expected)

    def test_concat_inserts_separator_when_first_file_lacks_trailing_newline(self) -> None:
        # This is the protective behavior the legacy implementation relied on:
        # if file A has no trailing newline, file B's first header must not
        # be glued to file A's last sequence.
        a = self.tmp / "a.faa"
        b = self.tmp / "b.faa"
        out = self.tmp / "combined.faa"
        _write(a, ">g1|p1\nMKT")  # no trailing newline
        _write(b, ">g2|p1\nCCC\n")

        _concat_proteomes([str(a), str(b)], str(out))

        recs = _records(out)
        self.assertEqual(sorted(recs), sorted([("g1|p1", "MKT"), ("g2|p1", "CCC")]))

    def test_concat_streams_without_loading_file_into_memory(self) -> None:
        # Proxy check: a larger file concatenates successfully.
        a = self.tmp / "big_a.faa"
        b = self.tmp / "big_b.faa"
        out = self.tmp / "combined.faa"
        buf_a = io.StringIO()
        buf_b = io.StringIO()
        for i in range(500):
            buf_a.write(f">gA|p{i}\n")
            buf_a.write("M" * 200 + "\n")
            buf_b.write(f">gB|p{i}\n")
            buf_b.write("K" * 200 + "\n")
        _write(a, buf_a.getvalue())
        _write(b, buf_b.getvalue())

        _concat_proteomes([str(a), str(b)], str(out))
        self.assertEqual(len(_records(out)), 1000)

    def test_concat_single_source_file_matches_bytes(self) -> None:
        a = self.tmp / "a.faa"
        out = self.tmp / "combined.faa"
        _write(a, ">g1|p1\nMKT\n>g1|p2\nAAA\n")

        _concat_proteomes([str(a)], str(out))

        self.assertEqual(out.read_bytes(), a.read_bytes())


if __name__ == "__main__":
    unittest.main()
