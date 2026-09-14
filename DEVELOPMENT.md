# Development

User instructions are in [`README.md`](README.md). This page covers source layout,
validation, and local benchmarks.

## Repository layout

```text
sgtree/
  src/sgtree/
    cli.py                  argument parsing and the run sequence
    config.py               Config dataclass; every output path is derived here
    input_stage.py          input intake, validation, gene calling
    fasta_normalize.py      header and alphabet normalization
    id_schema.py            the genome|contig|gene ID scheme
    search.py               marker search with pyhmmer, hit tables, filters
    extract.py              per-marker sequence extraction
    align.py                hmmalign, MAFFT, and FAMSA alignment
    duplicates.py           duplicate-copy resolution by bitscore
    supermatrix.py          trimming and concatenation
    phylogeny.py            VeryFastTree and IQ-TREE
    reference.py            reference genome preparation and caching
    cleanup.py              end-of-run archiving
    ani/                    pairwise ANI, MCL clustering, SNP trees
    ani_clustering.py       clustering entry point and representative choice
    marker_selection/       per-marker trees, duplicate cleanup, singleton scoring
    benchmarks/             benchmark generation and the evidence instruments
  bin/                      benchmark CLI entry points
  tests/                    unittest suite
  resources/models/         bundled marker-set HMM files
  testgenomes/              the bundled 10-genome example
  runs/                     local scratch output; not tracked
```

## Running the pipeline from source

The `sgtree` task runs the package from `src/` with `PYTHONPATH=src`:

```bash
pixi run sgtree --help
```

## Tests

```bash
pixi run test-unit
```

The task uses `unittest` discovery to collect all test modules.

Run the bundled example and check that the Python files compile:

```bash
pixi run example                                  # full pipeline to tree.nwk
pixi run python -m compileall -q src bin tests
```

`pixi run example` uses the fixed directory `runs/example_basic`. Do not run two copies
at the same time. A successful run prints the absolute path to
`runs/example_basic/tree.nwk`.

## Benchmarks

The benchmark tasks require local genome panels under `benchmarking/`. A fresh clone
does not include them. `pixi run benchmark-generate` builds the synthetic contamination
benchmark from `benchmarking/testgenomes/Chloroflexi` and writes results under `runs/`.

With the panels in place:

```bash
pixi run benchmark-generate
pixi run benchmark-run
```

The two evidence instruments below need no panel and run from a fresh clone.

`pixi run benchmark-prepare-burkholderiaceae` builds a 50-genome Burkholderiaceae
panel with taxonomy sidecars. It needs a local GTDB genome DuckDB; set its path
with the `SGTREE_TAXONOMY_DB` environment variable.

## Contamination-detection evidence instruments

Two fixed instruments gate the marker-discordance code. They test engineering and safety
properties and do not estimate biological performance. Their benchmark-only filters can
label a synthetic action but do not enable production pruning.

Tree fixture screen, 24 mechanism fixtures plus 8 scale fixtures, in memory:

```bash
env PYTHONPATH=src pixi run python -m sgtree.benchmarks.loo_tree_fixtures --check
```

Sequence benchmark, 12 held-out cases with real tree inference:

```bash
env PYTHONPATH=src pixi run python -m sgtree.benchmarks.loo_sequence_benchmark \
  --outdir runs/sequence_check --threads 1 --check
```

The sequence benchmark writes `per_event_comparison.tsv` with the scorer
comparison through the shared gate, budget, and RF pipeline, and
`review_tier.tsv` with the review candidates and their contig-vote and margin
evidence. `--sweep-donor-genes` traces the contig-gate operating curve over
donor gene counts into `donor_gene_sweep.tsv`.

### Scope of the evidence

The review tier reports marker copies for inspection and does not confirm contamination.
Its margin threshold was calibrated on two truth events in the development instrument,
which is too small to establish biological sensitivity or specificity. These fixed tests
check that the scorer and its safeguards behave as intended. Validation of automatic
removal requires independent genomes and contamination events.

## Cleaning the workspace

```bash
pixi run clean-runtime            # scratch output, keeps benchmark results
pixi run clean-benchmarks
pixi run clean-reference-cache
pixi run clean-all
```

## Pipeline stages

A basic run executes these stages in order:

1. Input intake, normalization, and gene calling for assembly input
2. Marker search with `pyhmmer`
3. Hit parsing, the marker count matrix, and the inclusion filters
4. Per-marker sequence extraction
5. Alignment
6. Duplicate-copy resolution by bitscore
7. Trimming with trimAl
8. Supermatrix concatenation
9. Species-tree inference
10. Archiving of intermediates

`--marker_selection yes` adds a second phase: per-marker trimming and tree inference,
RF-guided duplicate cleanup, optional singleton analysis, then a rebuild of the trimmed
alignments, the supermatrix, and the final tree. All singleton modes except `loo_profile`
can remove markers. `loo_profile` writes evidence and passes the marker trees through
unchanged.
`--selection_global_rounds` repeats that phase against the rebuilt guide tree.
