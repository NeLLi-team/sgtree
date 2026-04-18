"""Shared subprocess helpers used by the pipeline stages.

Two helpers intentionally, not one. See `docs/plans/active/2026-04-17-subprocess-audit.md`
for the call-site inventory. Callers have genuinely different needs:

- ``run_check`` for fire-and-check callers that only care whether the process
  succeeded. Raises ``CalledProcessError`` on non-zero exit.
- ``run_capture`` for callers that need to inspect ``stdout``/``stderr`` to make
  control-flow decisions (e.g. probing a tool's ``-threads`` capability) or
  that want captured output written to a file downstream. Returns the
  ``CompletedProcess`` without raising; caller decides.
"""

from __future__ import annotations

import logging
import subprocess
from collections.abc import Mapping
from typing import IO, Optional

logger = logging.getLogger("sgtree")


def run_check(
    cmd: list[str],
    *,
    stdout: Optional[IO | int] = None,
    stderr: Optional[IO | int] = None,
    env: Optional[Mapping[str, str]] = None,
    cwd: Optional[str] = None,
) -> None:
    """Run ``cmd`` with ``check=True``. Log and re-raise on failure.

    Callers pass ``stdout``/``stderr`` when they want to redirect or discard
    those streams. When not provided, the streams inherit the parent.
    """
    try:
        subprocess.run(
            cmd,
            stdout=stdout,
            stderr=stderr,
            env=dict(env) if env is not None else None,
            cwd=cwd,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        logger.error(
            "subprocess '%s' exited with code %s", cmd[0], exc.returncode
        )
        raise


def run_capture(
    cmd: list[str],
    *,
    env: Optional[Mapping[str, str]] = None,
    cwd: Optional[str] = None,
) -> subprocess.CompletedProcess[str]:
    """Run ``cmd`` with stdout/stderr captured as text. Does not raise.

    Returns the full ``CompletedProcess`` so the caller can branch on
    ``returncode`` or parse ``stderr`` (e.g. the VeryFastTree ``-threads``
    probe). The caller is responsible for raising if the result is not
    acceptable.
    """
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=dict(env) if env is not None else None,
        cwd=cwd,
        check=False,
    )
    logger.info(
        "subprocess '%s' completed with code %s", cmd[0], result.returncode
    )
    return result
