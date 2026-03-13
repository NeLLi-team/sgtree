#!/usr/bin/env python
"""Stage genome inputs and concatenate proteomes into a SGTree-safe FASTA."""

import argparse
import subprocess
import sys
from pathlib import Path

# Ensure local package imports work when invoked as a script.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from sgtree.fasta_normalize import normalize_and_concat_proteomes
from sgtree.input_stage import detect_input_format, gene_call_inputs, write_genome_manifest


def _maybe_reexec_with_pixi_python(args: argparse.Namespace) -> bool:
    pixi_python = ROOT / ".pixi" / "envs" / "default" / "bin" / "python"
    current_python = Path(sys.executable).resolve()
    if not pixi_python.exists() or current_python == pixi_python.resolve():
        return False
    cmd = [
        str(pixi_python),
        str(Path(__file__).resolve()),
        "--input",
        args.input,
        "--out",
        args.out,
        "--input-format",
        args.input_format,
    ]
    if args.map is not None:
        cmd.extend(["--map", args.map])
    if args.manifest is not None:
        cmd.extend(["--manifest", args.manifest])
    if args.gene_call_map is not None:
        cmd.extend(["--gene-call-map", args.gene_call_map])
    subprocess.run(cmd, check=True)
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input proteome directory (*.faa) or genome directory (*.fna)")
    parser.add_argument("--out", required=True, help="Normalized output FASTA")
    parser.add_argument("--map", default=None, help="Optional header mapping TSV")
    parser.add_argument("--manifest", default=None, help="Optional genome manifest TSV")
    parser.add_argument("--gene-call-map", default=None, help="Optional gene-call TSV for FNA inputs")
    parser.add_argument("--input-format", default="auto", choices=["auto", "faa", "fna"])
    args = parser.parse_args()

    input_format = args.input_format
    if input_format == "auto":
        input_format = detect_input_format(args.input)

    source_path = args.input
    staged_source = None
    if input_format == "fna":
        if args.gene_call_map is None:
            parser.error("--gene-call-map is required when staging FNA inputs")
        staged_source = str(Path(args.out).with_name("staged_proteomes"))
        try:
            gene_call_inputs(args.input, staged_source, args.gene_call_map)
        except RuntimeError as exc:
            if "pyrodigal" in str(exc) and _maybe_reexec_with_pixi_python(args):
                return
            raise
        source_path = staged_source

    if args.manifest:
        write_genome_manifest(
            args.input,
            input_format=input_format,
            manifest_path=args.manifest,
            staged_source=staged_source,
        )

    stats = normalize_and_concat_proteomes(source_path, args.out, args.map)
    print(
        "normalized proteomes: "
        f"format={input_format} genomes={stats['genomes']} records={stats['records']} "
        f"invalid_chars_replaced={stats['invalid_chars_replaced']}"
    )


if __name__ == "__main__":
    main()
