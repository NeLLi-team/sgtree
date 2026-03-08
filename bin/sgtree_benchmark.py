#!/usr/bin/env python
"""CLI wrapper for SGTree benchmark generation and evaluation."""
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.benchmarks import main


if __name__ == "__main__":
    main()
