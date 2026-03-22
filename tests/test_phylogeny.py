import subprocess
import unittest
from unittest.mock import Mock
from unittest.mock import patch

from sgtree import phylogeny


class PhylogenyTests(unittest.TestCase):
    def test_fasttree_executable_can_require_veryfasttree(self):
        which = Mock(side_effect=lambda name: {"VeryFastTree": "/tmp/VeryFastTree"}.get(name))

        with patch.dict("sgtree.phylogeny.os.environ", {"SGTREE_REQUIRE_VERYFASTTREE": "1"}, clear=False):
            with patch("sgtree.phylogeny.shutil.which", which):
                executable = phylogeny._fasttree_executable()

        self.assertEqual(executable, "VeryFastTree")

    def test_fasttree_executable_raises_when_veryfasttree_required_but_missing(self):
        with patch.dict("sgtree.phylogeny.os.environ", {"SGTREE_REQUIRE_VERYFASTTREE": "1"}, clear=False):
            with patch("sgtree.phylogeny.shutil.which", return_value=None):
                with self.assertRaises(FileNotFoundError):
                    phylogeny._fasttree_executable()

    def test_run_fasttree_retries_without_threads_when_binary_rejects_flag(self):
        calls = []

        def fake_run(cmd, stdout=None, stderr=None, check=None):
            calls.append(cmd)
            if len(calls) == 1:
                raise subprocess.CalledProcessError(
                    returncode=1,
                    cmd=cmd,
                    stderr=b"Unknown or incorrect use of option -threads",
                )
            return subprocess.CompletedProcess(cmd, 0)

        with patch.object(phylogeny, "_fasttree_executable", return_value="FastTree"):
            with patch("sgtree.phylogeny.subprocess.run", side_effect=fake_run):
                phylogeny.run_fasttree("input.faa", "tree.nwk", threads=24)

        self.assertEqual(
            calls,
            [
                ["FastTree", "-threads", "24", "-quiet", "-out", "tree.nwk", "input.faa"],
                ["FastTree", "-quiet", "-out", "tree.nwk", "input.faa"],
            ],
        )

    def test_build_tree_worker_retries_without_threads_when_binary_rejects_flag(self):
        calls = []

        def fake_run(cmd, stdout=None, stderr=None, check=None):
            calls.append(cmd)
            if len(calls) == 1:
                raise subprocess.CalledProcessError(
                    returncode=1,
                    cmd=cmd,
                    stderr=b"Unknown or incorrect use of option -threads",
                )
            return subprocess.CompletedProcess(cmd, 0)

        with patch.object(phylogeny, "_fasttree_executable", return_value="FastTree"):
            with patch("sgtree.phylogeny.subprocess.run", side_effect=fake_run):
                phylogeny._build_tree_worker(("input.faa", "treeouts", "fasttree", "LG", False))

        self.assertEqual(
            calls,
            [
                ["FastTree", "-threads", "1", "-quiet", "-out", "treeouts/input.faa_tree.out", "input.faa"],
                ["FastTree", "-quiet", "-out", "treeouts/input.faa_tree.out", "input.faa"],
            ],
        )


if __name__ == "__main__":
    unittest.main()
