# Provenance

This file explains where the retained benchmark metadata came from and what the retained summary tables do and do not represent.

## Dataset Provenance

### 50-gen

- Source root: `/home/fschulz/dev/software/sgtree_private_backup_20260316/runs/benchmarks/manuscript_panel_replicates_20260318_recipient_consensus`
- Status: exact
- Retained items:
  - benchmark manifests
  - scenario `events.tsv`
  - selected genomes and taxonomy
  - selected truth-marker panels
  - truth-marker completeness derived from scenario `reference_run/marker_count_matrix.csv`

### 100-gen

- Original bulky source root: `runs/contig_variant_proof_of_concept/benchmarks`
- Status: exact after regeneration
- Retained items:
  - regenerated benchmark manifests
  - exact scenario `events.tsv`
  - exact variant-specific contaminant placement metadata
  - selected genomes and taxonomy
  - selected truth-marker panels
  - truth-marker completeness derived from exact source scenario `reference_run/marker_count_matrix.csv`
- Source/recovery scripts retained under [100gen/scripts](/home/fschulz/dev/software/sgtree/eval/full/100gen/scripts)

## Result Provenance

- The retained summary tables in [summaries](/home/fschulz/dev/software/sgtree/eval/full/summaries) are copied from the last fully completed `recipient_consensus` benchmark rerun.
- That rerun completed on **March 21, 2026 at 11:08:44 PM PDT**.
- After that completed rerun, the exact 100-gen source benchmark panels were regenerated and the retained `eval/full/100gen` metadata was upgraded from reconstructed to exact.

## Important Distinction

- Dataset provenance is now exact for both `50gen` and `100gen`.
- Result provenance is still tied to the earlier completed rerun in `summaries/`.
- Therefore, the retained metadata package is internally consistent and exact, but a new full SGTree rerun is still required if you want result tables that are explicitly derived from the regenerated exact 100-gen source benchmark tree.
