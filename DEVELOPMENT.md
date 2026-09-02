# Development

Notes for people who work on SGTree. User documentation is in
[`README.md`](README.md).

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

Every task runs the package from `src/` with `PYTHONPATH=src`:

```bash
pixi run sgtree --help
```

## Tests

```bash
pixi run test-unit
```

The suite uses `unittest`. Use the full discovery command above; a shorter
invocation can miss test modules.

Two behavioral checks complete the gate:

```bash
pixi run example                                  # full pipeline to tree.nwk
pixi run python -m compileall -q src bin tests
```

## Benchmarks

The benchmark tasks need genome panels that the repository does not ship. They
live under `benchmarking/`, which is not tracked, so a fresh clone cannot run
them: `pixi run benchmark-generate` builds the synthetic contamination
benchmark from `benchmarking/testgenomes/Chloroflexi`. `runs/` is disposable
scratch.

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

Two fixed instruments gate the contamination-detection code. Both are engineering
and safety screens, not biological performance estimates, and both keep
production pruning disabled.

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

The review tier reports only. A review flag marks a marker copy for inspection
and confirms nothing on its own. The margin threshold that separates true from
false review warnings is calibrated on two truth events of the development
instrument. That is calibration, not validation, and it does not transfer to
empirical assemblies without a new confirmation run.

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

`--marker_selection yes` adds a second phase: per-marker trimming and tree
inference, RF-guided duplicate cleanup, optional singleton removal, then a
rebuild of the trimmed alignments, the supermatrix, and the final tree.
`--selection_global_rounds` repeats that phase against the rebuilt guide tree.
