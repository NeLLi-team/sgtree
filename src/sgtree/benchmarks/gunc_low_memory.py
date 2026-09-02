"""Run pinned GUNC with a smaller DIAMOND sequence block."""

from __future__ import annotations

import subprocess
import sys

from gunc import external_tools
from gunc.gunc import main as gunc_main


DIAMOND_BLOCK_SIZE = "1.0"


def _diamond_command(
    input_file: str,
    threads: int,
    temp_dir: str,
    database_file: str,
    out_file: str,
) -> list[str]:
    return [
        "diamond",
        "blastp",
        "--query",
        input_file,
        "--threads",
        str(threads),
        "--max-target-seqs",
        "1",
        "--masking",
        "0",
        "--evalue",
        "1",
        "--algo",
        "0",
        "--tmpdir",
        temp_dir,
        "--db",
        database_file,
        "--out",
        out_file,
        "--quiet",
        "--block-size",
        DIAMOND_BLOCK_SIZE,
    ]


def _diamond(
    input_file: str,
    threads: int,
    temp_dir: str,
    database_file: str,
    out_file: str,
) -> None:
    try:
        subprocess.check_output(
            _diamond_command(
                input_file,
                threads,
                temp_dir,
                database_file,
                out_file,
            ),
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"Low-memory DIAMOND failed for {input_file}: {error.output}"
        ) from error


def main() -> int:
    external_tools.diamond = _diamond
    result = gunc_main()
    return int(result or 0)


if __name__ == "__main__":
    raise SystemExit(main())
