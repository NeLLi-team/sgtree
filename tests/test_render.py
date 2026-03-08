import tempfile
import unittest
from pathlib import Path

from sgtree import render


class RenderTests(unittest.TestCase):
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
