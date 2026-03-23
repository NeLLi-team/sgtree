# Eval Datasets Latest

This directory keeps the compact benchmark dataset package after deleting bulky SGTree run outputs.

Contents:

- `50gen/`: exact manuscript-panel benchmark definitions for Flavo, Gamma, and Chlam replicate panels.
- `100gen/`: reconstructed Alpha and Bacteroidota proof-of-concept dataset metadata from the latest full completed rerun.
- `summaries/50gen_combined_summary.tsv`: per-run 50-gen summary table copied from the completed manuscript rerun root.
- `summaries/50gen_aggregate_summary.tsv`: aggregated 50-gen benchmark summary.
- `summaries/50gen_singleton_summary.tsv`: 50-gen singleton classification summary.
- `summaries/100gen_comparison.tsv`: final 100-gen comparison table from the latest full rerun.
- `manifest.json`: top-level package description and known limitations.
- `100gen/scripts/`: small source/recovery scripts preserved for provenance after the bulky 100-gen run trees are deleted.
- `PROVENANCE.md`: explicit source/result provenance notes for this package.

Interpretation notes:

- `introduced_markers.tsv` describes which contaminant marker was added to which recipient genome and, when recoverable, the donor genome and contaminant contig identifier.
- `genome_truth_marker_completeness.tsv` reports completeness against the benchmark truth-marker panel used for that dataset panel. These truth-marker panels are subsets of UNI56-derived markers selected for the benchmark.
- 50-gen event tables are exact copies from benchmark scenario definitions.
- 100-gen event tables are exact copies derived from regenerated source benchmark manifests and scenario event tables.
- In 100-gen panels, variant-specific contaminant contig placement is computed exactly from the source scenario events plus the documented rewrite rules for `solo_contig`, `merged_contam`, and `native_contig`.
- Exact per-marker reference-genome presence is available for both 50-gen and regenerated 100-gen panels.
- After this cleanup pass, `eval/full/` is the canonical retained benchmark-data package inside the repo; the bulky `runs/contig_variant_recipient_consensus_20260318` and regenerated `runs/contig_variant_proof_of_concept` trees are expected to be removed.
- The retained summary tables in `summaries/` are provenance-stable copies of the last completed benchmark rerun; they are not yet a fresh rerun against the regenerated exact 100-gen source package.
