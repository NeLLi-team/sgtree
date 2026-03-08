import importlib.util
import unittest
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = Path(__file__).resolve().parents[1] / "docs" / "figures" / "generate_workflow_figure.py"
SPEC = importlib.util.spec_from_file_location("workflow_figure", MODULE_PATH)
workflow_figure = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(workflow_figure)


class WorkflowFigureTests(unittest.TestCase):
    def test_normalize_spec_fills_missing_fields_from_defaults(self):
        spec = workflow_figure.normalize_spec(
            {
                "title": "Custom Title",
                "core_boxes": [{"title": "Only One", "lines": ["Line A"]}],
                "selection_boxes": [{"title": "Marker Trees", "lines": ["Per-marker trees"]}, {"title": "", "lines": []}],
                "benchmark_lines": ["One line only"],
            }
        )

        self.assertEqual(spec["title"], "Custom Title")
        self.assertEqual(spec["core_boxes"][0]["title"], "Only One")
        self.assertEqual(spec["core_boxes"][1]["title"], "Search and Parse")
        self.assertEqual(spec["selection_boxes"][1]["title"], "RF-guided Selection")
        self.assertEqual(spec["selection_boxes"][1]["lines"], ["Copy retained if it minimizes", "species-tree RF distance"])
        self.assertEqual(spec["benchmark_lines"], ["One line only"])

    def test_choose_spec_uses_openrouter_when_available(self):
        raw = {
            "title": "OpenRouter Title",
            "subtitle": "Subtitle",
            "core_heading": "Core",
            "core_boxes": [{"title": f"Core {idx}", "lines": [f"Line {idx}"]} for idx in range(5)],
            "selection_heading": "Selection",
            "selection_boxes": [{"title": f"Sel {idx}", "lines": [f"Sel Line {idx}"]} for idx in range(4)],
            "bridge_label": "bridge",
            "benchmark_heading": "Harness",
            "benchmark_lines": ["Line 1", "Line 2"],
        }
        with patch.object(workflow_figure, "load_openrouter_spec", return_value=raw):
            spec, planner = workflow_figure.choose_spec("auto")

        self.assertEqual(planner, "openrouter")
        self.assertEqual(spec["title"], "OpenRouter Title")
        self.assertEqual(spec["selection_boxes"][3]["title"], "Sel 3")

    def test_choose_spec_falls_back_to_static(self):
        with patch.object(workflow_figure, "load_openrouter_spec", return_value=None):
            spec, planner = workflow_figure.choose_spec("auto")

        self.assertEqual(planner, "static")
        self.assertEqual(spec["title"], workflow_figure.DEFAULT_SPEC["title"])


if __name__ == "__main__":
    unittest.main()
