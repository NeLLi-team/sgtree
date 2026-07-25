# Changelog

All notable changes to SGTree are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `--singles-mode loo_profile`, a report-only leave-one-marker-out
  placement-consensus scorer. It writes candidate evidence without pruning
  marker trees.
- A fixed tree benchmark includes 24 mechanism fixtures plus eight scale
  fixtures at 50 and 100 leaves. It checks clean-panel safety, far-source
  singleton detection, same-contig evidence, and RF guards.
- A fixed 12-case aligned-protein benchmark covers native-contig controls
  and held-out near-source and far-source events.

### Fixed

- Replacement-event evaluation now uses exact native and contaminant record
  IDs and treats missing or empty alignments as unknown outcomes.
- Singleton analysis runs after duplicate convergence, scores each candidate
  once, preserves cached branch support, and applies the RF guard to every
  automatic pruning mode.
- GCP low-panel fallback now runs recipient-consensus ranking and
  classification instead of applying the legacy classifier to every
  GCP candidate.
- Full GCP scoring now requires at least 5 markers; 3-4 marker panels
  use the recipient-consensus fallback because the per-genome z-score
  outlier-count gate is not reliable at that size.
- `--singles-mode` is now the preferred CLI spelling for singleton mode
  selection, while the legacy `--singles_mode` alias remains accepted.
- Runtime banner and logfile headers now report the package version
  (`SGTree 1.2.0`) instead of the stale `Sg_Tree v.2` label.
- README workflow diagram now includes the `famsa` alignment backend.
- Test output is quieter: expected negative-path stderr/log records and
  known third-party invalid-escape `SyntaxWarning`s are now captured or
  suppressed, and the trim fallback closes its parsed FASTA handle.

## [1.2.0] - 2026-04-20

### Highlights

- SGTree now uses a proper `src/sgtree` package layout with minimal
  `pyproject.toml` metadata instead of the previous flat package tree.
- The public bundled data footprint was reduced to one small
  10-genome FNA example dataset under `testgenomes/example_fna10`;
  benchmark-only source panels moved to the local ignored
  `benchmarking/` area.
- The built-in iTOL formatting and PNG tree rendering path was removed,
  and the dependency surface was reduced accordingly.
- `VeryFastTree` is now the only packaged fast tree backend, and
  `famsa` is available as a new optional de novo alignment backend.

### Added

- `pyproject.toml` with `src`-layout package discovery.
- `testgenomes/example_fna10/` as the single bundled FNA smoke-test
  panel.
- `--aln famsa` de novo alignment backend. SGTree runs FAMSA with
  `-refine_mode on` and the default guide-tree mode.

### Changed

- **Packaging**
  - The runtime package moved from `sgtree/` to `src/sgtree/`.
  - Repo-local wrappers under `bin/` now resolve imports through the
    `src` layout.
  - Pixi tasks now run the package via
    `env PYTHONPATH=src python -m sgtree ...`.

- **CLI and UX**
  - `--tree_method` now defaults to `veryfasttree`; `fasttree`
    remains accepted as a legacy alias.
  - Example task now uses `8` CPUs by default.
  - Runs now print explicit finalization progress after tree building,
    plus the absolute path to the final tree before the final timestamp
    line.

- **Data layout**
  - Public README/examples/tasks now point at the small FNA example
    dataset.
  - Benchmark preparation defaults now write local benchmark assets
    under `benchmarking/testgenomes/`.

- **Dependencies**
  - Version caps widened for `pandas`, `biopython`, `numpy`, and
    `xgboost`.
  - `VeryFastTree` replaces the old `FastTree` package fallback.

### Removed

- Built-in iTOL color-strip and heatmap output generation.
- Built-in PNG tree rendering and the related helper scripts/tests.
- Unused dependencies: `matplotlib`, `nbformat`, `nbclient`,
  `ipykernel`, `quarto`, `mummer4`, and `fasttree`.
- Obsolete root `sgtree.py` wrapper.

### Fixed

- The CLI now reports `VeryFastTree` consistently instead of showing
  the stale `fasttree` label while actually executing `VeryFastTree`.
- Long post-tree silent periods now show explicit cleanup/finalization
  progress so short example runs do not appear hung.

## [1.1.0] - 2026-04-18

### Highlights

- End-to-end pipeline throughput up ~22% on the reference example
  (~23 s → ~18 s) from streaming I/O and vectorized pandas paths.
- Marker-selection cleanup is now deterministic and parallel: the
  `marker_selection_rf_values.txt` file is written once from the parent
  process (no more concurrent-append race), and the neighbor-context
  builder is threaded across markers.
- `ani` and `marker_selection` are promoted to packages with pinned
  import surfaces, enabling a future modular split.

### Added

- `sgtree/_subprocess.py` — shared `run_check` / `run_capture` helpers.
  All modules that shell out to external binaries (`cleanup`,
  `phylogeny`, `align`, `reference`, `ani`, MCL) now use them, so
  subprocess invocations log consistently and surface stderr on
  failure.
- `sgtree/_fasta_utils.py::fasta_contig_bases_stats` consolidates the
  two previous duplicates in `input_stage.py` and `ani_clustering.py`.
- `sgtree/config.py::SingletonThresholds` frozen dataclass bundles ~30
  singleton-classifier constants with a drift test against the legacy
  module globals.
- `_cached_score_table` and `_cached_species_tree` process-local
  caches in `marker_selection`: workers parse the Newick species tree
  and load `table_elim_dups` at most once per worker lifetime instead
  of once per marker.
- Parallel per-marker neighbor-context builder
  (`_build_marker_neighbor_context` now accepts `num_cpus` and
  dispatches via a thread pool).
- New unit tests: extract (4), supermatrix (3), search (7), subprocess
  helper (6), fasta utils (3), duplicates (2), rf-values determinism
  (3), ANI thread scaling (1), phylogeny threads-probe cache (1),
  public API surface (2), singleton thresholds (1), ML-proposal
  decomposition regression (3), cached loaders (2), neighbor-context
  serial/parallel parity (2) — 120 → 127 tests total.
- `sgtree/benchmarks` gained `contig_variant_benchmark.py` coverage
  across genus/family/order/phylum scopes with cross-lineage donor
  support.

### Changed

- **Performance**
  - `extract.py`: proteome concatenation streams via `shutil.copyfileobj`
    with 1 MiB chunks instead of reading the whole file into memory.
  - `supermatrix.py`: `_fill_nan_gaps` is vectorized (preserves the
    bottom-most-width fallback for leading NaN runs); a single-pass
    `pd.concat` replaces the per-row append + CSV round-trip; the FASTA
    writer indexes the DataFrame once outside the loop.
  - `search.py`: vectorized hmmsearch parsing and dedup with
    `f"{proteinid}_{model}"` as the composite key (multi-domain
    semantics preserved).
  - `duplicates.py`: per-marker duplicate counting switched from
    `O(n²)` nested generators to a single `collections.Counter` pass.
  - `phylogeny.py`: VeryFastTree `-threads` capability probe is now
    cached in a module-level dict so the CLI no longer re-runs the
    fail-and-fall-back invocation for every marker.
  - `ani/__init__.py::compute_pairwise_ani`: minimap2 thread count
    scales to `max(1, num_cpus // min(num_cpus, n_pairs))` so small
    pair counts on large CPU budgets use wider per-call threading.
  - `marker_selection`: `score_neighbor_ml_proposals` split into
    `_prepare_ml_features`, `_fit_anomaly_models`, and
    `_rank_genome_proposals` (behavior preserved bitwise on the
    synthetic regression fixture).

- **Correctness**
  - `marker_selection`: RF-values workers return records to the parent,
    which sorts by `(marker, protein_id)` and writes the file once.
    Replaces the previous append-per-worker pattern that could
    interleave writes above `PIPE_BUF` on Linux and produced
    nondeterministic line ordering across runs.
  - `duplicates.py`: per-worker FASTA file handle is now closed
    deterministically (pre-existing leak under concurrent access
    fixed).
  - `cleanup.py`: zipping now writes to a temp path first so the
    self-zip path cannot truncate the source file if the archive write
    is interrupted.

- **Code layout**
  - `sgtree/marker_selection.py` and `sgtree/ani.py` promoted to
    packages via `git mv`; external imports and `patch.object`
    targets at the package namespace still work unchanged.
  - `tests/test_public_api.py` pins the exported surface of both
    packages so a future submodule split cannot silently drop a
    caller.

- **Error handling**
  - `logger.exception(...)` added before every swallowed or re-raised
    failure in `render.py`, `sgtree_logging.py`,
    `marker_selection/__init__.py`, and `cli.py`. Control flow is
    unchanged (every `raise` that was there before is still there).

- **Documentation**
  - README: `--singles_mode` table now lists `gcp` with the correct
    expansion ("Genome Consistency Profiling") and documents its
    low-panel fallback behavior. `topoknn` and `hybrid` are labeled as
    internal baselines. The Repository Structure section was brought
    back in sync with the actual module layout.
  - `tree_round_N.nwk` is documented as a per-round output.

### Removed

- `pixi.toml`: stale `manuscript-docx` task and unused `umap-learn`
  dependency.
- `cli.py`: duplicate `import shutil`, unreachable
  `print(sys.argv[:])` before logfile init.
- `duplicates.py::_is_duplicate` — zero callers in the repo.
- README: stale `clean-regression` reference.

### Fixed

- `cli.py`: GCP fallback log message corrected to match the classifier
  dispatch it triggers.
- `ani/__init__.py`: MCL invocation via `run_check` removed the
  invalid `text=True` kwarg that would have raised `TypeError` on any
  run that actually exercised the binary MCL path.

### Deferred

The following optional follow-ups are tracked for a later release:

- Full submodule split of `marker_selection/__init__.py` into
  `rf.py`, `scoring.py`, `classification.py`, `io.py` and of
  `ani/__init__.py` into `compute.py`, `cluster.py`, `snp_tree.py`.
  Blocked on a test restructuring because several tests patch names
  at the package namespace that internal callers would resolve via
  their own submodule globals after a split.
## [1.0.0] - 2026-03-03

Initial tagged release.
