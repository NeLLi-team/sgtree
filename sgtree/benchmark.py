"""Backward-compatible alias for the benchmark package."""

from __future__ import annotations

import sys

from . import benchmarks as _benchmarks


sys.modules[__name__] = _benchmarks
