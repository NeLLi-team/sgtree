"""The bundled marker sets must be usable.

`gtdbbac.hmm` and `gtdbarc.hmm` once mixed HMMER3/b and HMMER3/f model blocks
in one file. pyhmmer locks onto the format of the first block, so it stopped
after the leading models and the run failed before the marker search. Both
files parse end to end now, and this test keeps every bundled set that way.
"""

import re
import unittest
from pathlib import Path

import pyhmmer


MODEL_DIR = Path(__file__).resolve().parents[1] / "resources" / "models"


class BundledMarkerSetTests(unittest.TestCase):
    def test_every_bundled_set_parses_completely(self) -> None:
        paths = sorted(MODEL_DIR.glob("*.hmm"))
        self.assertTrue(paths, f"no marker sets under {MODEL_DIR}")
        for path in paths:
            with self.subTest(marker_set=path.name):
                text = path.read_text()
                declared = len(re.findall(r"^NAME", text, re.M))
                with pyhmmer.plan7.HMMFile(path) as handle:
                    parsed = sum(1 for _ in handle)
                self.assertEqual(
                    parsed, declared,
                    f"{path.name}: parsed {parsed} of {declared} declared models",
                )

    def test_each_set_uses_one_hmm_format(self) -> None:
        for path in sorted(MODEL_DIR.glob("*.hmm")):
            with self.subTest(marker_set=path.name):
                formats = set(re.findall(r"^HMMER3/[a-z]", path.read_text(), re.M))
                self.assertEqual(
                    len(formats), 1,
                    f"{path.name} mixes model formats: {sorted(formats)}",
                )


if __name__ == "__main__":
    unittest.main()
