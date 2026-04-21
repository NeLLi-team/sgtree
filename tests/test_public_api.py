"""Pin the import surface of `sgtree.marker_selection` and `sgtree.ani`.

These symbols are imported by the pipeline (`src/sgtree/cli.py`), the
benchmarks layer, and external tooling under `bin/`. The Phase 4
structural splits (module -> package) must preserve every name here.
A regression that removes or renames an export breaks downstream
callers silently at import time; this test makes that loud.
"""

from __future__ import annotations

import importlib
import unittest


_MARKER_SELECTION_PUBLIC = [
    # from bin/remove_singles.py
    "_load_score_table",
    "_tree_to_genome_level",
    "choose_tree_by_rf",
    "prune_singletons",
    # from src/sgtree/cli.py
    "run_noperm",
    "remove_singles",
    "write_cleaned_sequences",
    # from tests/test_marker_selection.py
    "resolve_marker_tree",
]

_ANI_PUBLIC = [
    # from tests/test_ani_clustering.py
    "GenomeRecord",
    "_merge_intervals",
    "choose_cluster_representative",
    # from src/sgtree/ani_clustering.py (ani_core.*)
    "run_ani_clustering",
    "build_snp_trees",
]


class PublicApiTests(unittest.TestCase):
    def test_marker_selection_exports(self):
        module = importlib.import_module("sgtree.marker_selection")
        for name in _MARKER_SELECTION_PUBLIC:
            with self.subTest(symbol=name):
                self.assertTrue(
                    hasattr(module, name),
                    f"sgtree.marker_selection.{name} missing after refactor",
                )

    def test_ani_exports(self):
        module = importlib.import_module("sgtree.ani")
        for name in _ANI_PUBLIC:
            with self.subTest(symbol=name):
                self.assertTrue(
                    hasattr(module, name),
                    f"sgtree.ani.{name} missing after refactor",
                )


if __name__ == "__main__":
    unittest.main()
