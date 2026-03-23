# Eval

This directory is the fresh benchmark-evaluation root.

- `full/`: complete compact metadata package for the 3x50-gen and 2x100-gen benchmark families, including all retained panel/variation metadata.
- `pilot/`: one combined representative dataset per lineage family for fast inspection and writeup work.

Consistency status:

- `50gen/` is exact benchmark metadata copied from the original manuscript-panel benchmark manifests and event tables.
- `100gen/` is reconstructed from the latest full completed rerun because the original source benchmark directories were no longer present locally when this package was built.
- For `100gen/`, exact source manifests, exact scenario event tables, exact variant-specific contaminant placement, genome membership, and exact per-marker reference-genome presence are retained after regeneration.

Provenance status:

- The retained dataset package under `eval/full/50gen` is sourced from the manuscript-panel benchmark manifests under the private backup root.
- The retained dataset package under `eval/full/100gen` was later upgraded to exact source metadata after regenerating the 24 proof-of-concept benchmark panels locally.
- The retained result summaries under `eval/full/summaries/` still come from the last fully completed benchmark rerun that finished on March 21, 2026.
- A fresh full SGTree rerun against the regenerated exact 100-gen source panels has not yet been completed, so dataset provenance is now exact while results provenance still points to the earlier completed rerun.

Why this directory is only about `19M`:

- It contains TSV/JSON metadata only.
- It does not contain the bulky SGTree run trees, alignments, per-run proteome blobs, or intermediate tables that previously dominated disk usage.
- The removed run trees were tens of gigabytes because each run stored repeated `proteomes`, `ref_and_query_proteomes`, and header-map files.

The bulky SGTree run trees and intermediate benchmark outputs were removed after packaging this metadata.
