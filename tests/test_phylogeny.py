import subprocess
import unittest
from unittest.mock import Mock
from unittest.mock import patch

from sgtree import phylogeny


class PhylogenyTests(unittest.TestCase):
    def test_fasttree_executable_prefers_veryfasttree(self):
        which = Mock(side_effect=lambda name: {"VeryFastTree": "/tmp/VeryFastTree"}.get(name))
        with patch("sgtree.phylogeny.shutil.which", which):
            executable = phylogeny._fasttree_executable()

        self.assertEqual(executable, "VeryFastTree")

    def test_fasttree_executable_raises_when_veryfasttree_missing(self):
        with patch("sgtree.phylogeny.shutil.which", return_value=None):
            with self.assertRaises(FileNotFoundError):
                phylogeny._fasttree_executable()

    def setUp(self) -> None:
        phylogeny._THREADS_SUPPORT.clear()

    def test_run_fasttree_retries_without_threads_when_binary_rejects_flag(self):
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            if len(calls) == 1:
                # First call: run_capture(capture_output=True, text=True, check=False).
                return subprocess.CompletedProcess(
                    cmd, 1, stdout="", stderr="Unknown or incorrect use of option -threads"
                )
            # Second call: run_check(stdout=PIPE, stderr=PIPE, check=True).
            return subprocess.CompletedProcess(cmd, 0)

        with patch.object(phylogeny, "_fasttree_executable", return_value="VeryFastTree"):
            with patch("sgtree._subprocess.subprocess.run", side_effect=fake_run):
                phylogeny.run_fasttree("input.faa", "tree.nwk", threads=24)

        self.assertEqual(
            calls,
            [
                ["VeryFastTree", "-threads", "24", "-quiet", "-out", "tree.nwk", "input.faa"],
                ["VeryFastTree", "-quiet", "-out", "tree.nwk", "input.faa"],
            ],
        )

    def test_threads_support_probe_is_cached_after_first_failure(self):
        # First call: binary rejects -threads -> both threaded (capture) and
        # fallback (check) subprocess.run invocations fire.
        # Second call: cache says -threads unsupported -> only fallback fires.
        phylogeny._THREADS_SUPPORT.clear()
        calls: list[list[str]] = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            if cmd[1] == "-threads":
                return subprocess.CompletedProcess(
                    cmd, 1, stdout="", stderr="Unknown or incorrect use of option -threads"
                )
            return subprocess.CompletedProcess(cmd, 0)

        with patch.object(phylogeny, "_fasttree_executable", return_value="VeryFastTree"):
            with patch("sgtree._subprocess.subprocess.run", side_effect=fake_run):
                phylogeny.run_fasttree("input.faa", "tree.nwk", threads=4)
                call_count_after_first = len(calls)
                phylogeny.run_fasttree("input2.faa", "tree2.nwk", threads=4)

        # First call: one threaded + one fallback = 2 subprocess.run calls.
        self.assertEqual(call_count_after_first, 2)
        # Second call: cache hit -> fallback only, so +1 total.
        self.assertEqual(len(calls), 3)
        self.assertEqual(calls[2][0], "VeryFastTree")
        self.assertNotIn("-threads", calls[2])
        # Cleanup cache so other tests get a fresh probe.
        phylogeny._THREADS_SUPPORT.clear()

    def test_build_tree_worker_retries_without_threads_when_binary_rejects_flag(self):
        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            if len(calls) == 1:
                return subprocess.CompletedProcess(
                    cmd, 1, stdout="", stderr="Unknown or incorrect use of option -threads"
                )
            return subprocess.CompletedProcess(cmd, 0)

        with patch.object(phylogeny, "_fasttree_executable", return_value="VeryFastTree"):
            with patch("sgtree._subprocess.subprocess.run", side_effect=fake_run):
                phylogeny._build_tree_worker(("input.faa", "treeouts", "fasttree", "LG", False))

        self.assertEqual(
            calls,
            [
                ["VeryFastTree", "-threads", "1", "-quiet", "-out", "treeouts/input.faa_tree.out", "input.faa"],
                ["VeryFastTree", "-quiet", "-out", "treeouts/input.faa_tree.out", "input.faa"],
            ],
        )


if __name__ == "__main__":
    unittest.main()
