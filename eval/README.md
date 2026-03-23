# Eval

This directory is the fresh benchmark-evaluation root.

- `full/`: complete compact metadata package for the 3x50-gen and 2x100-gen benchmark families, including all retained panel/variation metadata.
- `pilot/`: one combined representative dataset per lineage family for fast inspection and writeup work.

Consistency status:

- `50gen/` is exact benchmark metadata copied from the original manuscript-panel benchmark manifests and event tables.
- `100gen/` is reconstructed from the latest full completed rerun because the original source benchmark directories were no longer present locally when this package was built.
- For `100gen/`, exact source manifests, exact scenario event tables, exact variant-specific contaminant placement, genome membership, and exact per-marker reference-genome presence are retained after regeneration.

Why this directory is only about `19M`:

- It contains TSV/JSON metadata only.
- It does not contain the bulky SGTree run trees, alignments, per-run proteome blobs, or intermediate tables that previously dominated disk usage.
- The removed run trees were tens of gigabytes because each run stored repeated `proteomes`, `ref_and_query_proteomes`, and header-map files.

The bulky SGTree run trees and intermediate benchmark outputs were removed after packaging this metadata.
