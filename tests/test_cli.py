import sys
import unittest
from unittest.mock import patch

from sgtree.cli import parse_args


class CliTests(unittest.TestCase):
    def test_snp_defaults_to_disabled(self):
        argv = ["sgtree", "input_dir", "models.hmm"]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertFalse(cfg.snp)
        self.assertFalse(cfg.ani_cluster)

    def test_snp_requires_ani_cluster(self):
        argv = ["sgtree", "input_dir", "models.hmm", "--snp", "yes"]
        with patch.object(sys, "argv", argv):
            with self.assertRaisesRegex(ValueError, "--snp requires --ani_cluster yes"):
                parse_args()

    def test_snp_enabled_with_ani_cluster(self):
        argv = [
            "sgtree",
            "input_dir",
            "models.hmm",
            "--ani_cluster",
            "yes",
            "--snp",
            "yes",
        ]
        with patch.object(sys, "argv", argv):
            cfg = parse_args()
        self.assertTrue(cfg.ani_cluster)
        self.assertTrue(cfg.snp)


if __name__ == "__main__":
    unittest.main()
