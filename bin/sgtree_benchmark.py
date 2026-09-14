#!/usr/bin/env python
"""CLI wrapper for SGTree benchmark generation and evaluation."""

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"


def main() -> None:
    """Load and run the benchmark CLI from the source checkout."""
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))
    benchmarks = importlib.import_module("sgtree.benchmarks")
    benchmarks.main()


if __name__ == "__main__":
    raise SystemExit(main())
