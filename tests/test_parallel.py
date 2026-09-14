"""Regression tests for multiprocessing fallback behavior."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sgtree.parallel import map_processed


def _record_then_maybe_fail(payload: tuple[str, int]) -> int:
    record_path, value = payload
    with Path(record_path).open("a", encoding="utf-8") as handle:
        handle.write(f"{value}\n")
    if value == 2:
        raise OSError("worker failure")
    return value


class ProcessMappingTests(unittest.TestCase):
    def test_worker_oserror_propagates_without_repeating_completed_work(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            record = Path(tmpdir) / "calls.txt"
            arguments = [(str(record), 1), (str(record), 2)]

            with self.assertRaisesRegex(OSError, "worker failure"):
                map_processed(_record_then_maybe_fail, arguments, workers=2)

            self.assertEqual(
                record.read_text(encoding="utf-8").splitlines().count("1"), 1
            )

    def test_pool_construction_error_falls_back_to_serial_execution(self) -> None:
        with (
            patch("sgtree.parallel.mp.Pool", side_effect=PermissionError("denied")),
            self.assertLogs("sgtree.parallel", level="WARNING") as logs,
        ):
            result = map_processed(abs, [-1, -2], workers=2)

        self.assertEqual(result, [1, 2])
        self.assertIn("falling back to serial execution", logs.output[0])


if __name__ == "__main__":
    unittest.main()
