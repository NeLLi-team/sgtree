import unittest

from sgtree.benchmarks import gunc_low_memory


class GuncLowMemoryTests(unittest.TestCase):
    def test_diamond_command_changes_only_the_block_size(self):
        command = gunc_low_memory._diamond_command(
            "query.faa",
            4,
            "tmp",
            "database.dmnd",
            "hits.tsv",
        )

        self.assertEqual(command[:2], ["diamond", "blastp"])
        self.assertEqual(command[-2:], ["--block-size", "1.0"])
        self.assertEqual(command[command.index("--threads") + 1], "4")
        self.assertEqual(
            command[command.index("--max-target-seqs") + 1],
            "1",
        )
        self.assertEqual(command[command.index("--masking") + 1], "0")
        self.assertEqual(command[command.index("--evalue") + 1], "1")
        self.assertEqual(command[command.index("--algo") + 1], "0")


if __name__ == "__main__":
    unittest.main()
