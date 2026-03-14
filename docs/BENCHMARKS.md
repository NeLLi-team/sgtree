# SGTree Taxonomy-Aware Benchmarking

This document describes the current SGTree taxonomy-aware contamination benchmark. The benchmark notebook that regenerates the exported tables and figure is [`docs/benchmark_summary.ipynb`](benchmark_summary.ipynb), and the executed notebook is [`docs/benchmark_summary.executed.ipynb`](benchmark_summary.executed.ipynb).

This document covers the contamination benchmark suite under `runs/benchmarks/`. The separate ANI/SNP Burkholderiaceae runtime benchmark and its dataset-preparation workflow are documented in [`README.md`](/home/fschulz/dev/software/sgtree/README.md).

## Benchmark Structure

The benchmark is organized around manifest-defined panels rather than a fixed six-dataset list.

- Rank-based panels:
  - `genus`
  - `family`
  - `order`
- Mixed high-level panels:
  - `mixed_high_level`

Each benchmark directory writes a [`benchmark_manifest.json`](/home/fschulz/dev/software/sgtree/runs/benchmarks) that records:

- `lineage_label`
- `taxonomic_scope`
- selected truth genomes and truth markers
- donor lineage metadata
- per-scenario event tables
- per-scenario genome summaries

Suite directories additionally write a `suite_manifest.json` that lists the benchmark directories belonging to that lineage/scope suite.

## Current Lineage Panels

The active benchmark families are:

- `flavo`
- `gamma`
- `chlam`

The mixed high-level design reuses those three lineages with explicit donor pairs:

- Flavo + Chlam into Gamma
- Gamma + Chlam into Flavo
- Flavo + Gamma into Chlam

## Contamination Scenarios

Rank-based panels use three fixed contamination scenarios that map directly to SGTree cleanup modes:

| Scenario | Cleanup profile | Meaning |
|---|---|---|
| `duplicate_only` | `duplicate_cleanup` | duplicate-marker contamination only |
| `replacement_only` | `singles_delta_rf` | singleton-replacement contamination only |
| `combined` | `duplicate_plus_singles_delta_rf` | duplicate plus singleton-replacement contamination |

Mixed panels use:

| Scenario | Cleanup profile | Meaning |
|---|---|---|
| `mixed_high_level` | `duplicate_plus_singles_delta_rf` | class/phylum-level mixed donor contamination with duplicate and replacement recipients |

Mixed panels also encode the requested overlap layout in the manifest:

- `duplicate_recipients`
- `replacement_recipients`
- `overlap_recipients`

## Exported Data Products

The notebook and export helper now write three stable TSVs under [`docs/data`](/home/fschulz/dev/software/sgtree/docs/data):

- [`benchmark_dataset_overview.tsv`](/home/fschulz/dev/software/sgtree/docs/data/benchmark_dataset_overview.tsv)
  - one row per panel/scenario pair
  - carries `panel_id`, `panel_label`, `lineage_label`, `taxonomic_scope`, `taxonomic_scope_label`, donor metadata, scenario counts, and mixed-recipient layout fields
- [`benchmark_genome_contamination.tsv`](/home/fschulz/dev/software/sgtree/docs/data/benchmark_genome_contamination.tsv)
  - one row per affected recipient genome
  - carries panel/scenario metadata plus per-genome contamination counts
- [`benchmark_summary_all.tsv`](/home/fschulz/dev/software/sgtree/docs/data/benchmark_summary_all.tsv)
  - one row per completed panel/scenario benchmark run
  - carries RF, contaminant-removal, taxa-retention, and singleton-pruning outcomes
- [`benchmark_alignment_comparison.tsv`](/home/fschulz/dev/software/sgtree/docs/data/benchmark_alignment_comparison.tsv)
  - one row per panel/scenario comparison between the baseline `results/summary.tsv` and the latest `results_mafft_*/summary.tsv`
  - carries baseline vs comparison RF, contaminant-removal, runtime, and missing-taxa counts together with signed deltas
- [`benchmark_alignment_comparison_summary.tsv`](/home/fschulz/dev/software/sgtree/docs/data/benchmark_alignment_comparison_summary.tsv)
  - one row per benchmark panel summarizing the detailed comparison table
  - carries counts of improved/unchanged/worsened scenarios plus mean RF and runtime shifts

These TSVs are the preferred manuscript-facing source of truth because they are regenerated from manifests and `results/summary.tsv` files rather than from hardcoded dataset assumptions.

## Alignment Replay Summary

The March 13, 2026 replay added a side-by-side comparison between the pre-existing benchmark baseline under each panel's `results/summary.tsv` and a full MAFFT rerun under `results_mafft_20260313/summary.tsv`. The comparison spans the 12 completed lineage and mixed high-level panels documented above, for a total of 30 scenario-level comparisons.

At the scenario level, MAFFT was worse in `25/30` comparisons, unchanged in `2/30`, and better in `3/30`. The three MAFFT improvements were limited to `gamma/genus/replacement_only`, `chlam/order/replacement_only`, and `chlam/order/combined`. The broad pattern was therefore unfavorable to switching the default alignment mode away from profile-guided alignment.

Lineage-specific behavior was uneven:

- `flavo` regressed in all `10/10` comparisons. The worst regression was `flavo/order/combined`, which shifted from `0.0213` to `0.1915` final normalized RF. Runtime increased by roughly `1.13x` to `1.30x`.
- `gamma` was the least affected lineage. Most RF shifts were small, usually around `+0.0233`, and the only clear MAFFT improvement in the Gamma panels was `genus/replacement_only`, which improved from `0.0465` to `0.0233`. Runtime increased by about `1.05x` to `1.09x`.
- `chlam` was mixed. `order/replacement_only` and `order/combined` improved modestly, but `family` and `genus` mostly worsened, including `chlam/family/combined` from `0.0000` to `0.1277`. Runtime increased by roughly `1.42x` to `1.64x`.

Those results motivated restoring `hmmalign` as the SGTree default. MAFFT remains available as an explicit alternative (`--aln mafft` or `--aln mafft-linsi`) for cases where de novo alignment is still desired.

## Event Tables

Per-scenario `events.tsv` files carry the full traceability needed for supplement or audit use:

- `taxonomic_scope`
- `taxonomic_scope_label`
- `recipient_*` taxonomy fields
- `donor_*` taxonomy fields
- accession, organism name, record ids, and contig ids

The event schema is richer than the manuscript TSVs and is intended for detailed provenance rather than the main narrative summary.

## Figure Outputs

The notebook regenerates these figure files under [`docs/figures`](/home/fschulz/dev/software/sgtree/docs/figures):

- [`benchmark_summary.png`](/home/fschulz/dev/software/sgtree/docs/figures/benchmark_summary.png)
- [`benchmark_summary.svg`](/home/fschulz/dev/software/sgtree/docs/figures/benchmark_summary.svg)
- [`benchmark_summary_manuscript.png`](/home/fschulz/dev/software/sgtree/docs/figures/benchmark_summary_manuscript.png)
- [`benchmark_summary_manuscript.svg`](/home/fschulz/dev/software/sgtree/docs/figures/benchmark_summary_manuscript.svg)

## Validation

Core validation commands:

- `pixi run python -m unittest tests.test_benchmark tests.test_input_stage tests.test_marker_selection tests.test_benchmark_dataset tests.test_ani tests.test_ani_clustering tests.test_cli`
- `pixi run python bin/sgtree_benchmark.py export-docs --benchmarks-root runs/benchmarks --outdir docs/data`

Notebook execution example:

```bash
pixi run python -m ipykernel install --user --name sgtree-pixi --display-name "SGTree (pixi current)"
pixi run python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

src = Path("docs/benchmark_summary.ipynb")
out = Path("docs/benchmark_summary.executed.ipynb")
nb = nbformat.read(src, as_version=4)
NotebookClient(nb, timeout=1800, kernel_name="sgtree-pixi").execute(cwd=str(src.parent))
nbformat.write(nb, out)
PY
```

## Current Caveat

The benchmark catalog is now manifest-driven, so partially completed benchmark dirs can exist while long runs are still in progress. The notebook/export path only reflects completed runs that already have `results/summary.tsv`.

## ANI/SNP Runtime Benchmark

The Burkholderiaceae ANI/SNP benchmark is intentionally separate from the contamination suite above. It uses a real 50-genome assembly panel prepared with `pixi run benchmark-prepare-burkholderiaceae`, then runs SGTree with ANI clustering enabled and SNP trees optionally enabled via `--snp yes`. This benchmark measures three things that the contamination catalog does not: representative collapse after 95% ANI clustering, strain-level SNP resolution inside retained species clusters, and the new contig-level contamination guard applied before SNP alignment.

In the current verified run, the 50-genome Burkholderiaceae panel collapsed to 20 ANI representatives. One 26-member cluster produced a 524-site SNP alignment after filtering to shared cluster-core UNI56-bearing contigs that also aligned back to the representative backbone at or above 95% ANI. A second 6-member cluster passed the same contig filter but contained no variable SNP sites and was recorded as `no_snp_sites`. Each multi-member SNP cluster writes both `members.tsv` and `contig_filter.tsv`, together with the filtered contig FASTA files actually used for alignment, so the retained backbone is traceable at the contig level rather than only at the whole-genome level.
