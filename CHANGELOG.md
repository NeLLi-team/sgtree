# Changelog

All notable changes to SGTree are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.2] - 2026-09-02

### Fixed

- `pixi run benchmark-generate` and `pixi run benchmark-run` work. Both spawn
  `python -m sgtree` as a child process, and the package runs from the source
  tree rather than an installed distribution, so the child exited with
  "No module named sgtree". The child now inherits the source directory on
  `PYTHONPATH`.

### Changed

- `DEVELOPMENT.md` states that the benchmark tasks need genome panels the
  repository does not ship, so a fresh clone cannot run them.

## [1.2.1] - 2026-09-02

### Fixed

- `gtdbbac.hmm` and `gtdbarc.hmm` are usable. Both files mixed `HMMER3/b` and
  `HMMER3/f` model blocks, and the parser locks onto the format of the first
  block, so a run with either set stopped after the leading models and failed
  in the marker search. Both are rewritten in one format. Model names, model
  lengths, and emission probabilities are unchanged, and a test now parses
  every bundled marker set end to end.

### Added

- `LICENSE`: the Berkeley Lab non-commercial use licence, matching the other
  NeLLi-team repositories. `pyproject.toml` records it.

## [1.2.0] - 2026-09-02

### Added

- Review-warning margin discrimination. Review warnings now require the
  informative contig votes' mean attachment-background bitscore margin to
  reach `REVIEW_MIN_VOTE_MARGIN` (28.0, calibrated on the two truth events
  of the v2 development instrument; calibration, not validation). Donor
  genes show a sharp affinity cliff at the donor-clade boundary; native
  genes sit on a smooth divergence gradient. Measured on v2: the plain gate
  passes 42 of 91 candidates (40 false); the margin floor keeps the 2 truth
  warnings and rejects all 40. `review_tier.tsv` gains `vote_margin_mean`,
  `vote_margin_min`, and `margin_pass`; the report keeps gate-only counts
  for comparison. Action-path gates are unchanged.
- Sequence-benchmark instrument v2 (`panel_realism: native_contig_genes_v2`).
  Every native marker contig now carries the gate-minimum three non-marker
  genes, so contig audits on native records are decided by vote direction
  instead of gene scarcity. Measured consequence, pinned in tests: the
  review tier's contig gate passes 42 of 91 review candidates: 2 true and
  40 false (8 of the false in clean panels), a weakness the gene-scarce v1
  instrument hid; action decisions,
  sentinel vetoes, and all frozen structural gates are unchanged.
- Donor gene-count sweep (`--sweep-donor-genes`): traces the contig-gate
  operating curve over 0-10 donor genes per event. The gate is a step
  function at the three-query floor; donor votes agree fully at every
  count.
- Greedy voter-selection tier: a 16-marker far-source panel scored outside
  the frozen 12-case matrix (own cache directory) proves the production
  greedy path detects, gates, and removes the far-source event and spares
  the sentinel. New `loo_voter_search_mode` column in LOO rows,
  `singleton_candidates.tsv`, and `review_tier.tsv`.
- Report-only LOO review tier. `score_loo_profiles` now emits `loo_robust_z`
  (MAD-scaled conflict) and `loo_review_candidate` (dispersion-ceiling
  abstentions whose conflict is still robustly outside the voter MAD band,
  fixed threshold z >= 3). `singleton_candidates.tsv` carries both columns;
  the flag is topology-only in production and confirms nothing by itself.
  The 12-case sequence benchmark audits review candidates through the
  existing contig-vote gate and reports gate-confirmed review warnings in
  `review_tier.tsv` plus a `review_tier` report section. Action decisions,
  frozen comparison tables, and check gates are unchanged.
- `--singles-mode loo_profile`, a report-only leave-one-marker-out
  placement-consensus scorer. It writes candidate evidence without pruning
  marker trees.
- A fixed tree benchmark includes 24 mechanism fixtures plus eight scale
  fixtures at 50 and 100 leaves. It checks clean-panel safety, far-source
  singleton detection, same-contig evidence, and RF guards.
- A fixed 12-case aligned-protein benchmark covers native-contig controls
  and held-out near-source and far-source events.

### Changed

- `README.md` is rewritten for users: install, quick start, input, common runs,
  a complete option reference, outputs, and how marker selection works.
  Repository layout, tests, benchmarks, and the evidence instruments moved to
  `DEVELOPMENT.md`.
- The taxonomy database path for the benchmark generators comes from
  `SGTREE_TAXONOMY_DB`. The package shipped a hardcoded local path, and
  `--taxonomy-db` is now required when the variable is unset.

### Fixed

- Output-directory cleanup no longer deletes files SGTree did not write. The
  cleanup of a previous run is restricted to the entries SGTree generates, and
  it now runs before ANI clustering and reference preparation write into the
  output directory. Re-running with `--ani_cluster yes` into a finished output
  directory previously deleted `ani/` and then failed with `FileNotFoundError`.
- Reference-cache archiving destroyed the files it archived. Opening a path for
  writing as a zip truncated it first, so every archived member was zero bytes.
- Directory archiving is recursive. A marker-selection run archived `protTrees/`
  with two empty directory entries and then removed the originals, losing every
  per-marker tree. Archiving uses `zipfile` instead of the `zip` binary, which
  was never a declared dependency and fails on an empty directory.
- End-of-run archiving leaves user files alone, and no longer overwrites or
  relocates a zip file the user placed in the output directory.
- `proteomes_header_map.tsv`, `gene_calls.tsv`, `singleton_candidates.tsv`, and
  the per-round trees stay readable text. They were compressed in place and kept
  their original names, so opening them returned zip data.
- SNP alleles on reverse-strand alignments were wrong. SAM stores SEQ in
  reference orientation, and the parser reverse-complemented it a second time,
  so every CIGAR-derived base index pointed at the wrong position.
- `minimap2` is a declared dependency. `--snp yes` and `--ani_backend minimap2`
  failed with `FileNotFoundError` after the species tree had already been built.
- A protein that passes the threshold of several markers is assigned to the
  marker it scores best against, with a deterministic tie-break. It was assigned
  by HMM file order.
- The optional length filter checks the exit status of `grep`. A failure
  replaced the hit table with a partial file.
- Marker names that contain `_` or `.` work. Marker-selection runs raised
  `FileNotFoundError` for 4 of the 12 bundled marker sets (`gtdbbac`, `gtdbarc`,
  `mitomarkers108`, `COX123`), and on `gtdbarc` two distinct markers collapsed
  onto one name.
- `--singles_min_rfdist` applies when `--ref` is set. The reference path skipped
  the gate, so every query leaf with a positive RF gain became a proposal.
- Reference genome names that contain dots keep their protection from pruning.
- `singleton_candidates.tsv` reports what happened. A leaf the RF guard kept was
  recorded as `pruned`, and a candidate classified `clean` was recorded as
  `blocked_by_genome_budget`.
- A header with an empty ID field no longer raises `IndexError`.
- `--aln` rejects an unknown value instead of running `hmmalign` silently.
- IQ-TREE reruns on an existing prefix, so `--selection_global_rounds 2` with
  `--tree_method iqtree` completes.
- A failing external tool reports its captured stderr, not only an exit code.
- Missing `trimal` or `mcl` reports the fallback it used.
- Marker-tree scoring reads its input in a fixed order, so `gcp` and
  `neighbor_ml` scores no longer depend on directory order.


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

### Removed

- Fourteen orphaned `bin/` scripts left over from the removed Nextflow
  runtime. Their logic lives in `src/sgtree/` and nothing invoked them; the
  standalone copies had drifted behind the package implementations.
- Dormant synthetic-benchmark machinery in `sgtree.benchmarks`: the
  mixed-lineage and cross-family generators, the taxonomic suite wrapper,
  panel/seed replicate generation, multi-suite running, run aggregation, and
  the docs table exporter, plus the matching `sgtree_benchmark.py`
  subcommands and Pixi tasks. `generate`, `generate-taxonomic`,
  `prepare-burkholderiaceae`, and `run` remain. Scenario names, cleanup
  profiles, and `evaluate_benchmark_run` are unchanged, so existing panels
  and result artifacts still evaluate.
- `SingletonThresholds`/`DEFAULT_SINGLETON_THRESHOLDS` in `sgtree.config`.
  The bundle was never consumed by the classifier and had drifted from the
  live module constants it mirrored.
- Dead code: `phylogeny.run_snp_tree`, `marker_selection.prune_singletons`,
  the always-true `singleton_mode_uses_global_rf_gate` guard and the
  `effective_singleton_mode` pass-through, four write-only `Config` path
  fields, and assorted unused imports and helpers.

## [1.2.0-dev] - 2026-04-20

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
