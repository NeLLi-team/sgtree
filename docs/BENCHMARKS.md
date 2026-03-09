# SGTree Systematic Benchmarking

This document describes the current manuscript benchmark design. The benchmark track is now limited to three explicit contamination scenarios that map directly to the cleanup modes used in SGTree:

- `duplicate_only`: added extra copies of existing markers, no singleton replacement
- `replacement_only`: replacement of correct single-copy markers, no added duplicate copies
- `combined`: both added duplicate copies and singleton replacement

The benchmark notebook that generates the exported tables and figure is [`docs/benchmark_summary.ipynb`](benchmark_summary.ipynb), and the executed notebook is [`docs/benchmark_summary.executed.ipynb`](benchmark_summary.executed.ipynb).

## Dataset Families

The current benchmark set contains six datasets:

| Dataset | Family | Size | Selected genomes | Truth markers | Donor scope |
|---|---|---:|---:|---:|---|
| Flavo Small | Flavobacteriaceae | small | 18 | 10 | within-group only |
| Flavo Large | Flavobacteriaceae | large | 45 | 15 | within-group only |
| Gamma Small | Gammaproteobacteria | small | 14 | 13 | within-group only |
| Gamma Large | Gammaproteobacteria | large | 31 | 10 | within-group only |
| Cross Small | Flavo ↔ Gamma | small | 33 | 6 | cross-group only |
| Cross Large | Flavo ↔ Gamma | large | 79 | 8 | cross-group only |

The full dataset/scenario table is exported to [`docs/data/benchmark_dataset_overview.tsv`](data/benchmark_dataset_overview.tsv). The per-genome contamination table is exported to [`docs/data/benchmark_genome_contamination.tsv`](data/benchmark_genome_contamination.tsv).

## What Is Measured

For each dataset and scenario, the benchmark output now records:

- the clean truth tree built from genomes without added contaminant markers
- the initial RF of the contaminated tree before cleanup
- the final RF after applying the scenario-appropriate cleanup mode
- the RF delta (`initial_tree_rf_norm - tree_rf_norm`)
- the number of contaminant markers removed out of the number introduced
- duplicate-event counts and replacement-event counts separately

The results table is [`docs/data/benchmark_summary_all.tsv`](data/benchmark_summary_all.tsv).

## Scenario-to-Mode Mapping

The benchmark no longer sweeps over five cleanup heuristics. Instead, it evaluates only the cleanup mode that matches each contamination class:

| Scenario | Cleanup profile | Meaning |
|---|---|---|
| `duplicate_only` | `duplicate_cleanup` | duplicate-marker cleanup only |
| `replacement_only` | `singles_neighbor` | singleton-replacement cleanup with `singles_neighbor` |
| `combined` | `duplicate_plus_singles_neighbor` | duplicate cleanup plus `singles_neighbor` |

This matches the manuscript use case directly: one duplicate-only cleanup path, one singleton-replacement cleanup path, and one combined cleanup path.

## Current Results

The benchmark figure is [`docs/figures/benchmark_summary.svg`](figures/benchmark_summary.svg) and [`docs/figures/benchmark_summary.png`](figures/benchmark_summary.png). It shows, for each scenario, the initial RF, the final RF after cleanup, the RF delta, and the contaminant-removed fraction.

The main numerical results are:

| Dataset | Duplicate only: RF change / removed | Replacement only: RF change / removed | Combined: RF change / removed |
|---|---|---|---|
| Flavo Small | `0.133 → 0.067` (`4/8`) | `0.133 → 0.067` (`4/4`) | `0.200 → 0.000` (`6/12`) |
| Flavo Large | `0.048 → 0.000` (`4/8`) | `0.024 → 0.048` (`4/4`) | `0.000 → 0.024` (`4/12`) |
| Gamma Small | `0.091 → 0.091` (`6/8`) | `0.000 → 0.000` (`4/4`) | `0.182 → 0.000` (`7/12`) |
| Gamma Large | `0.036 → 0.036` (`4/8`) | `0.036 → 0.000` (`3/4`) | `0.071 → 0.000` (`11/12`) |
| Cross Small | `0.333 → 0.333` (`3/8`) | `0.367 → 0.172` (`4/4`) | `0.567 → 0.400` (`5/12`) |
| Cross Large | `0.145 → 0.000` (`8/8`) | `0.197 → 0.108` (`4/4`) | `0.145 → 0.066` (`10/12`) |

## Interpretation

Within-group added-copy contamination is now handled more effectively than in the earlier pre-fix runs. Flavo Large duplicate-only cleanup now reaches `0.048 → 0.000` while removing `4/8` added contaminant copies, and Cross Large duplicate-only cleanup reaches `0.145 → 0.000` with `8/8` duplicates removed. Duplicate-only cleanup still has little topological effect in Gamma Small and Gamma Large, and Cross Small remains difficult.

Singleton replacement is no longer hardest in the cross-clade panels. After the singleton-filter update, both cross-clade replacement-only datasets remove `4/4` contaminants and improve RF (`0.367 → 0.172` in Cross Small; `0.197 → 0.108` in Cross Large), which is consistent with the larger divergence between donor and recipient clades. The more difficult replacement-only cases are now the within-group panels, especially Flavo Large.

The combined scenario remains the most informative. Gamma Large is the strongest combined case by contaminant removal, with `11/12` contaminants removed and RF improving from `0.071` to `0.000`. Cross Large is the clearest cross-clade success case, improving from `0.145` to `0.066` while removing `10/12` contaminant markers. The weakest combined case remains Flavo Large, where the final tree is slightly worse than the initial contaminated tree despite removal of all four singleton replacements.

## Validation

The benchmark implementation was validated with the repository unit suite and the Python-versus-Nextflow regression harness:

- `pixi run test-unit`
- `pixi run test-regression-clean`

The new benchmark notebook also executes successfully top-to-bottom with:

```bash
pixi run python /home/fschulz/dev/omics-skills/skills/notebook-ai-agents-skill/scripts/execute_notebook.py \
  docs/benchmark_summary.ipynb \
  --out docs/benchmark_summary.executed.ipynb \
  --kernel sgtree \
  --timeout 1200
```
