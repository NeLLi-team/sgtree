import io
import shutil
import sys
import tempfile
import unittest
import zipfile
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from sgtree._subprocess import run_check
from sgtree.cleanup import cleanup_basic, cleanup_marker_selection
from sgtree.cli import parse_args
from sgtree.input_stage import write_genome_manifest
from sgtree.reference import check_duplicate_proteomes, prepare_reference


HEADER_MAP = "source_file\tcontig_header\nassembly.fna\tcontig_1\n"
GENE_CALLS = "genome_id\tgene_id\tbegin\tend\ng1\tgene_000001\t1\t99\n"
SINGLETONS = "genome_id\tmarker\tloo_robust_z\ng1\tPF00001\t3.4\n"
NEWICK = "((a:0.1,b:0.2):0.3,c:0.4);\n"


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def _archived_member(path: Path) -> zipfile.ZipInfo:
    """Return the single member of an archive, failing if there is not exactly one."""
    with zipfile.ZipFile(path) as archive:
        infos = archive.infolist()
        assert len(infos) == 1, f"{path} holds {len(infos)} members"
        return infos[0]


def _archived_text(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        return archive.read(archive.namelist()[0]).decode()


class ReferenceArchivingTests(unittest.TestCase):
    """The reference cache must survive its own archiving step."""

    def setUp(self) -> None:
        self._tmp = tempfile.mkdtemp()
        self.tmp = Path(self._tmp)
        genomedir = self.tmp / "query"
        refdir = self.tmp / "ref"
        genomedir.mkdir()
        refdir.mkdir()
        _write(genomedir / "q1.faa", ">q1|c1|gene_000001\nMKV\n")
        _write(refdir / "r1.faa", ">r1|c1|gene_000001\nMKV\n")

        argv = [
            "sgtree", str(genomedir), str(self.tmp / "UNI56.hmm"),
            "--save_dir", str(self.tmp / "out"),
            "--ref", str(refdir),
            "--ref_concat", str(self.tmp / "cache"),
        ]
        with patch.object(sys, "argv", argv):
            self.cfg = parse_args()
        self.ref_dir = Path(self.cfg.ref_dir_path())

        with patch("sgtree.reference.run_check", side_effect=self._fake_run_check):
            with redirect_stdout(io.StringIO()):
                prepare_reference(self.cfg)

    def tearDown(self) -> None:
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _fake_run_check(self, cmd, **kwargs):
        """Stand in for the reference sgtree run; run real `zip` calls."""
        if cmd[0] == "zip":
            return run_check(cmd, **kwargs)
        self.ref_dir.mkdir(parents=True)
        _write(self.ref_dir / "log_genomes_removed.txt", "genome_a\ngenome_b\n")
        _write(self.ref_dir / "tree.nwk", NEWICK)
        _write(self.ref_dir / "ref_and_query_proteomes", ">r1|c1|gene_000001\nMKV\n")
        _write(self.ref_dir / "proteomes_header_map.tsv", HEADER_MAP)
        _write(self.ref_dir / "hits.hmmout", "# hmmsearch output\n")
        _write(self.ref_dir / "table_elim_dups", "savedname,namemodel,score_bits\n")
        (self.ref_dir / "models").mkdir()
        _write(self.ref_dir / "models" / "PF00001.hmm", "HMMER3/f\n")
        return None

    def test_archived_files_keep_their_contents(self) -> None:
        cases = {
            self.ref_dir / "temp" / "log_genomes_removed.txt": "genome_a\ngenome_b\n",
            self.ref_dir / "temp" / "tree.nwk": NEWICK,
            self.ref_dir / "temp" / "ref_and_query_proteomes": ">r1|c1|gene_000001\nMKV\n",
            self.ref_dir / "proteomes_header_map.tsv": HEADER_MAP,
        }
        for path, expected in cases.items():
            with self.subTest(archive=path.name):
                self.assertTrue(path.exists(), f"{path} disappeared")
                self.assertTrue(zipfile.is_zipfile(path), f"{path} is not an archive")
                member = _archived_member(path)
                self.assertEqual(member.file_size, len(expected.encode()))
                self.assertEqual(_archived_text(path), expected)

    def test_kept_reference_files_stay_plain(self) -> None:
        for name in ("hits.hmmout", "table_elim_dups"):
            with self.subTest(name=name):
                path = self.ref_dir / name
                self.assertFalse(zipfile.is_zipfile(path))
                self.assertTrue(path.read_text().strip())

    def test_directory_archive_is_not_archived_a_second_time(self) -> None:
        archive = self.ref_dir / "temp" / "models.zip"
        self.assertTrue(zipfile.is_zipfile(archive))
        with zipfile.ZipFile(archive) as zf:
            names = zf.namelist()
            self.assertTrue(
                any(name.endswith("models/PF00001.hmm") for name in names), names
            )
            self.assertFalse(any(name.endswith(".zip") for name in names), names)
            member = next(n for n in names if n.endswith("models/PF00001.hmm"))
            self.assertEqual(zf.read(member).decode(), "HMMER3/f\n")


class ReferenceGenomeIdTests(unittest.TestCase):
    """Reference genome ids must match the ids input_stage derives."""

    def setUp(self) -> None:
        self._tmp = tempfile.mkdtemp()
        self.tmp = Path(self._tmp)

    def tearDown(self) -> None:
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _dir_with(self, name: str, *filenames: str) -> Path:
        path = self.tmp / name
        path.mkdir()
        for filename in filenames:
            _write(path / filename, ">g|c|gene_000001\nMKV\n")
        return path

    def test_ids_keep_interior_dots_and_match_input_stage(self) -> None:
        refdir = self._dir_with("ref", "GCF_000005845.2.faa", "strain.A.faa")
        empty_query = self._dir_with("query")

        ls_ref = check_duplicate_proteomes(str(empty_query), str(refdir))
        self.assertEqual(ls_ref, ["GCF_000005845.2", "strain.A"])

        manifest = self.tmp / "genome_manifest.tsv"
        write_genome_manifest(
            str(refdir), input_format="faa", manifest_path=str(manifest)
        )
        rows = manifest.read_text().splitlines()[1:]
        self.assertEqual(ls_ref, [row.split("\t")[0] for row in rows])

    def test_shared_prefix_is_not_reported_as_a_duplicate(self) -> None:
        query = self._dir_with("query", "GCF_000005845.1.faa")
        refdir = self._dir_with("ref", "GCF_000005845.2.faa")
        self.assertEqual(
            check_duplicate_proteomes(str(query), str(refdir)), ["GCF_000005845.2"]
        )

    def test_identical_genome_in_both_directories_still_fails(self) -> None:
        query = self._dir_with("query", "GCF_000005845.2.faa")
        refdir = self._dir_with("ref", "GCF_000005845.2.faa")
        with redirect_stdout(io.StringIO()) as out:
            with self.assertRaises(SystemExit):
                check_duplicate_proteomes(str(query), str(refdir))
        self.assertIn("GCF_000005845.2", out.getvalue())


class CleanupOutputTests(unittest.TestCase):
    """Documented outputs must stay readable after a run is cleaned up."""

    def setUp(self) -> None:
        self._tmp = tempfile.mkdtemp()
        self.outdir = Path(self._tmp)

    def tearDown(self) -> None:
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _populate(self) -> None:
        _write(self.outdir / "proteomes_header_map.tsv", HEADER_MAP)
        _write(self.outdir / "gene_calls.tsv", GENE_CALLS)
        _write(self.outdir / "singleton_candidates.tsv", SINGLETONS)
        _write(self.outdir / "table_elim_dups", "savedname,namemodel,score_bits\n")

    def _assert_plain(self, name: str, expected: str) -> None:
        path = self.outdir / name
        self.assertTrue(path.exists(), f"{name} left the output directory")
        self.assertNotEqual(path.read_bytes()[:2], b"PK", f"{name} is a zip archive")
        self.assertEqual(path.read_text(), expected)

    def _assert_archived(self, name: str) -> None:
        path = self.outdir / "temp" / name
        self.assertTrue(path.exists(), f"{name} was not archived")
        self.assertTrue(zipfile.is_zipfile(path), f"{name} was not archived")

    def test_cleanup_basic_keeps_documented_outputs_readable(self) -> None:
        self._populate()
        cleanup_basic(str(self.outdir))

        self._assert_plain("proteomes_header_map.tsv", HEADER_MAP)
        self._assert_plain("gene_calls.tsv", GENE_CALLS)
        self._assert_plain("singleton_candidates.tsv", SINGLETONS)
        self._assert_archived("table_elim_dups")

    def test_user_zip_is_not_overwritten_or_moved(self) -> None:
        """A user's own .zip in the output directory is left where it is."""
        self._populate()
        (self.outdir / "protTrees").mkdir()
        _write(self.outdir / "protTrees" / "t.nw", NEWICK)
        mine = _write(self.outdir / "my_results.zip", "not really a zip\n")

        cleanup_marker_selection(str(self.outdir))

        self.assertTrue(mine.exists(), "user zip was moved into temp/")
        self.assertEqual(mine.read_text(), "not really a zip\n")
        self.assertTrue((self.outdir / "temp" / "protTrees.zip").exists())

    def test_nested_directories_survive_archiving(self) -> None:
        """The directory archive is recursive, so nested trees are not lost."""
        self._populate()
        nested = self.outdir / "protTrees" / "no_duplicates" / "out"
        nested.mkdir(parents=True)
        _write(nested / "_no_dups_COG0012_.nw", NEWICK)
        _write(self.outdir / "protTrees" / "top.txt", "top\n")

        cleanup_marker_selection(str(self.outdir))

        archive = self.outdir / "temp" / "protTrees.zip"
        self.assertTrue(archive.exists(), "protTrees was not archived")
        with zipfile.ZipFile(archive) as zf:
            stored = {name for name in zf.namelist() if not name.endswith("/")}
            self.assertIn("protTrees/no_duplicates/out/_no_dups_COG0012_.nw", stored)
            self.assertIn("protTrees/top.txt", stored)
            self.assertEqual(
                zf.read("protTrees/no_duplicates/out/_no_dups_COG0012_.nw").decode(),
                NEWICK,
            )

    def test_user_files_survive_cleanup(self) -> None:
        """Cleanup archives only what SGTree wrote, so user data is untouched."""
        self._populate()
        notes = _write(self.outdir / "my_notes.md", "keep me\n")
        data_dir = self.outdir / "my_data"
        data_dir.mkdir()
        precious = _write(data_dir / "precious.csv", "a,b\n1,2\n")

        cleanup_basic(str(self.outdir))

        self.assertEqual(notes.read_text(), "keep me\n")
        self.assertNotEqual(notes.read_bytes()[:2], b"PK")
        self.assertTrue(data_dir.is_dir(), "user directory was archived away")
        self.assertEqual(precious.read_text(), "a,b\n1,2\n")
        self.assertFalse((self.outdir / "temp" / "my_data.zip").exists())
        # the generated intermediate is still archived
        self._assert_archived("table_elim_dups")

    def test_user_files_survive_marker_selection_cleanup(self) -> None:
        self._populate()
        notes = _write(self.outdir / "my_notes.md", "keep me\n")
        cleanup_marker_selection(str(self.outdir))

        self.assertEqual(notes.read_text(), "keep me\n")
        self._assert_archived("table_elim_dups")

    def test_cleanup_marker_selection_keeps_documented_outputs_readable(self) -> None:
        self._populate()
        _write(self.outdir / "tree_round_1.nwk", NEWICK)
        _write(self.outdir / "tree.nwk", NEWICK)
        cleanup_marker_selection(str(self.outdir))

        self._assert_plain("proteomes_header_map.tsv", HEADER_MAP)
        self._assert_plain("gene_calls.tsv", GENE_CALLS)
        self._assert_plain("singleton_candidates.tsv", SINGLETONS)
        self._assert_plain("tree_round_1.nwk", NEWICK)
        self._assert_archived("table_elim_dups")
        self.assertFalse((self.outdir / "tree.nwk").exists())


if __name__ == "__main__":
    unittest.main()
