#!/usr/bin/env python
"""Render a Newick tree as PNG with color-coded leaf labels.

If rendering fails for any reason, write a placeholder PNG and warning sidecar
instead of failing the pipeline.
"""
import argparse
import base64
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.render import render_tree_file

PLACEHOLDER_PNG_B64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO7+9pUAAAAASUVORK5CYII="
)


def write_placeholder_png(out_png: str, reason: str) -> None:
    with open(out_png, "wb") as handle:
        handle.write(base64.b64decode(PLACEHOLDER_PNG_B64))
    warning_path = f"{out_png}.warning.txt"
    with open(warning_path, "w") as handle:
        handle.write(f"{reason}\n")
    print(f"warning: {reason}", file=sys.stderr)
    print(f"warning: wrote placeholder image {out_png}", file=sys.stderr)
    print(f"warning: details in {warning_path}", file=sys.stderr)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True)
    parser.add_argument("--color", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    try:
        render_tree_file(args.tree, args.color, args.out)
    except Exception as exc:
        write_placeholder_png(args.out, f"tree rendering skipped: {exc}")


if __name__ == "__main__":
    main()
