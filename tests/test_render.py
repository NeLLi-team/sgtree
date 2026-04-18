import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sgtree import render


class RenderTests(unittest.TestCase):
    def test_midpoint_rooting_failure_logs_warning_and_continues(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            tree_path = tmp / "tree.nwk"
            color_path = tmp / "color.txt"
            out_png = tmp / "tree.png"

            tree_path.write_text("((A:0.1,B:0.1):0.2,(C:0.3,D:0.3):0.2);\n")
            color_path.write_text(
                "\n".join(
                    [
                        "DATASET_COLORSTRIP",
                        "SEPARATOR SPACE",
                        "DATASET_LABEL label1",
                        "COLOR #ff0000",
                        "COLOR_BRANCHES 0",
                        "DATA",
                    ]
                )
                + "\n"
            )

            # Force root_at_midpoint to fail; render must continue and still
            # produce the PNG, while logging the failure at WARNING.
            from Bio.Phylo.BaseTree import Tree as BioTree

            original_root = BioTree.root_at_midpoint

            def _raise(self):
                raise RuntimeError("simulated midpoint failure")

            with patch.object(BioTree, "root_at_midpoint", _raise):
                with self.assertLogs("sgtree", level="WARNING") as cm:
                    render.render_tree_file(str(tree_path), str(color_path), str(out_png))

            self.assertTrue(any("midpoint rooting failed" in m for m in cm.output))
            self.assertTrue(out_png.exists())

    def test_render_tree_file_creates_png_without_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            tree_path = tmp / "tree.nwk"
            color_path = tmp / "color.txt"
            out_png = tmp / "tree.png"

            tree_path.write_text("((A:0.1,B:0.1):0.2,C:0.3);\n")
            color_path.write_text(
                "\n".join(
                    [
                        "DATASET_COLORSTRIP",
                        "SEPARATOR SPACE",
                        "DATASET_LABEL label1",
                        "COLOR #ff0000",
                        "COLOR_BRANCHES 0",
                        "DATA",
                        "A #FF0000",
                        "B #00AA00",
                        "C #0000FF",
                    ]
                )
                + "\n"
            )

            render.render_tree_file(str(tree_path), str(color_path), str(out_png))

            self.assertTrue(out_png.exists())
            self.assertGreater(out_png.stat().st_size, 0)
            self.assertFalse((tmp / "tree.png.warning.txt").exists())


if __name__ == "__main__":
    unittest.main()
