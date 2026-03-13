# SGTree

SGTree is an end-to-end workflow for phylogenetic tree building. Use the provided sets of HMMs or provide your own HMMs to find the proteins of interest. In the default branch, duplicate copies are resolved by HMM bitscore before species-tree inference. In `--marker_selection` mode, SGTree builds per-marker trees and performs RF-guided duplicate cleanup plus optional singleton removal to improve contamination robustness.

## Setup

Install the Pixi environment:

```bash
pixi install
```

The environment is managed through `pixi.toml` and `pixi.lock`.

## Run

Primary interface (`pixi run sgtree` wrapper):

```bash
pixi run sgtree --help
```
Basic run:

```bash
pixi run sgtree \
  --genomedir <path to dir with protein faa files, one faa file per genome> \
  --modeldir <path to marker set .hmm>
```

Example run:

```bash
pixi run sgtree \
  --genomedir testgenomes/Chloroflexi \
  --modeldir resources/models/UNI56.hmm
```

Marker-selection run with references and singleton filtering:

```bash
pixi run sgtree \
  --genomedir testgenomes/Chloroflexi \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/nextflow/manual_full \
  --marker_selection true \
  --ref testgenomes/chlorref \
  --singles yes
```

Legacy baseline for direct comparisons:

```bash
pixi run sgtree \
  --genomedir testgenomes/Chloroflexi \
  --modeldir resources/models/UNI56.hmm \
  --marker_selection true \
  --ref testgenomes/chlorref \
  --selection_mode legacy
```

`pixi run sgtree` writes logs automatically to `runs/nextflow/logs/`.
For nucleotide assembly input (`*.fna`, `*.fa`, `*.fasta`) or any run with `--ani_cluster yes` or `--snp yes`, the wrapper dispatches to the Python engine automatically so gene calling, ANI clustering, and optional SNP-tree generation all happen in one code path.
Marker searches and `--aln hmmalign` are run with `pyhmmer` (HMMER-compatible search output).
Rendered PNG trees now use a headless-safe matplotlib/Biopython renderer.

Example with IQ-TREE and explicit HMM threshold mode:

```bash
pixi run sgtree \
  --genomedir testgenomes/Chloroflexi \
  --modeldir resources/models/UNI56.hmm \
  --tree_method iqtree \
  --iqtree_fast true \
  --hmmsearch_cutoff cut_ga
```

Example with `fna` input plus ANI clustering and opt-in per-cluster SNP trees:

```bash
pixi run benchmark-prepare-burkholderiaceae

pixi run sgtree \
  --genomedir testgenomes/Burkholderiaceae50 \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/burkholderiaceae_ani_benchmark \
  --num_cpus 24 \
  --ani_cluster true \
  --snp true \
  --ani_threshold 95 \
  --marker_selection true \
  --keep_intermediates true
```

Second choice (Python implementation without nextflow):

```bash
pixi run sgtree-python testgenomes/Chloroflexi resources/models/UNI56.hmm --num_cpus 8
```

Backward-compatible wrapper:

```bash
pixi run ./sgtree.py testgenomes/Chloroflexi resources/models/UNI56.hmm --num_cpus 8
```

## Settings

Core method controls:

- `--aln`: `hmmalign`, `mafft`, or `mafft-linsi` (default `hmmalign`).
- `--tree_method`: `fasttree` or `iqtree` (default `fasttree`) for both species tree and per-marker trees.
- `--iqtree_fast`: apply `-fast` when `--tree_method iqtree` (default `true`).
- `--iqtree_model`: IQ-TREE model string (default `LG+F+I+G4`).
- `--selection_mode`: `coordinate` or `legacy` (default `coordinate`) for marker-selection duplicate cleanup.
- `--selection_max_rounds`: maximum coordinate-descent rounds in `coordinate` mode (default `5`).
- `--selection_global_rounds`: rebuild the guide species tree and rerun duplicate cleanup for a small fixed number of rounds (default `1`).
- `--lock_references`: keep reference duplicate resolution score-locked instead of RF-updating them (default `false`).
- `--singles_mode`: `neighbor`, `delta_rf`, `backbone`, or `ensemble` when singleton filtering is enabled.
- `--singles_min_rfdist`: minimum marker/global RF distance required before singleton pruning activates (default `0.25`).
- `--keep_intermediates`: keep intermediate alignments/tables for debugging and benchmarking (default `false`).
- `--ani_cluster`: run pairwise ANI on the combined query+reference genome set and keep one representative per cluster for the main SGTree species tree.
- `--snp`: build cluster-level SNP trees after ANI clustering (default `false`; requires `--ani_cluster yes`). Before SNP alignment, SGTree keeps only contigs that carry shared cluster-core UNI56 markers and that still align back to the representative backbone at `>=95%` ANI.
- `--ani_threshold`: ANI cutoff used to retain graph edges before clustering (default `95`).
- `--ani_backend`: `auto`, `skani`, or `minimap2` (default `auto`; SGTree prefers `skani` when available and falls back to `minimap2` in restricted environments).
- `--ani_mcl_inflation`: MCL inflation used for ANI graph clustering (default `2.0`).
- `--snp_tree_min_cluster_size`: build a cluster-level SNP tree only when an ANI cluster has at least this many genomes (default `3`).

HMM search thresholds:

- `--hmmsearch_cutoff cut_ga`: use model gathering cutoffs (recommended for curated marker sets such as UNI56).
- `--hmmsearch_cutoff cut_tc`: use model trusted cutoffs.
- `--hmmsearch_cutoff cut_nc`: use model noise cutoffs.
- `--hmmsearch_cutoff evalue --hmmsearch_evalue <float>`: use a plain E-value threshold.

Genome inclusion/exclusion criteria:

- `--percent_models` (default `10`): minimum fraction of markers detected per genome.
- `--max_sdup` (default `-1`): maximum allowed copies of any single marker in one genome; `-1` disables.
- `--max_dupl` (default `-1`): maximum allowed fraction of markers present in multiple copies; `-1` disables.
- `--lflt` (default `0`): optional per-marker length filter (% of median hit length).
- `--num_nei` (default `0`): optional singleton-removal neighbor count override (`0` keeps auto mode).

nsgtree-style mapping:

- `minmarker` -> `--percent_models` (fraction mapped to percent).
- `maxsdup` -> `--max_sdup`.
- `maxdupl` -> `--max_dupl`.
- `hmmsearch_cutoff` -> `--hmmsearch_cutoff` and `--hmmsearch_evalue`.
- `tmethod` -> `--tree_method`.
- `iq_*` model controls -> `--iqtree_model` (and `--iqtree_fast`).
- `mafftv`/`mafft` -> `--aln mafft` or `--aln mafft-linsi` (or `--aln hmmalign`).

Practical selection guide:

- Curated marker sets (for example UNI56): start with `--hmmsearch_cutoff cut_ga`.
- Less curated/custom marker sets: start with `--hmmsearch_cutoff evalue --hmmsearch_evalue 1e-5`, then tighten if false positives appear.
- `--aln hmmalign` is the fastest stable default and keeps alignment behavior tied to each profile HMM.
- `--aln mafft-linsi` is slower but can help when marker-specific profile alignment is not desired.
- `--tree_method fasttree` is the quick default; `--tree_method iqtree --iqtree_fast true` is a practical higher-accuracy option.
- `--selection_mode coordinate` is the stronger default; `legacy` is kept for benchmark comparisons.
- `--selection_global_rounds 2` is the current practical setting for harder contamination benchmarks.
- `--singles yes` is still heuristic. On the current hard small benchmark the best topology is obtained with iterative duplicate cleanup and singleton filtering off; on the larger Flavobacteriaceae prototype, `--singles_mode neighbor` is currently the strongest singleton-aware option.
- Typical inclusion presets:
- Balanced: `--percent_models 10 --max_sdup 2 --max_dupl 0.25`
- Strict: `--percent_models 30 --max_sdup 1 --max_dupl 0.10`
- Relaxed: `--percent_models 5 --max_sdup -1 --max_dupl -1`

## Input Requirements

SGTree accepts either proteome FASTA (`*.faa`) or genome assembly FASTA (`*.fna`, `*.fa`, `*.fasta`). Inputs are normalized internally to:

```text
>IMG2684622718|2685462912
MLCAFAEEEAKIAETVGKVATELKVKKLLSDFATKEGEEHISTYNKIAMTAKAEGYADIEAMLCAFAEEEAKLQKL
```

Normalization behavior:

- Directory input (`--genomedir <dir>`): one genome/proteome per file; genome id is derived from filename stem.
- Single FASTA input (`--genomedir <file>`): if headers already contain `genome|protein`, the genome part is preserved.
- `fna` input is gene-called with `pyrodigal` before marker search; the gene-call table is written as `gene_calls.tsv`.
- Headers and IDs are sanitized to avoid delimiter collisions.
- Malformed header joins (for example `...*>next_header`) are repaired before parsing.
- Invalid amino-acid characters are replaced with `X`; `*` is removed.
- Header mapping is written as `proteomes_header_map_<input>.tsv` in `--outdir`.
- A genome manifest (`genome_manifest.tsv`) is written for every run and is reused by ANI clustering and optional SNP-tree generation.

## Output Structure

Wrapper / Python output (`--outdir` or `--save_dir`):

```text
<outdir>/
  tree.nwk
  tree_final.nwk                 # marker-selection mode
  tree_final.png                 # marker-selection mode
  marker_count_matrix.csv
  marker_count.txt               # basic mode
  marker_counts.txt              # marker-selection mode
  marker_selection_rf_values.txt # marker-selection mode
  color.txt
  log_genomes_removed.txt
  genome_manifest.tsv
  proteomes_header_map_<input>.tsv
  ani/
    ani_pairwise.tsv
    ani_clusters.tsv
    ani_representatives.tsv
    ani_kept_genomes.txt
    ani_graph.tsv
    ani_mcl_clusters.txt
  snp_trees/                     # only with --snp yes
    snp_tree_summary.tsv
    <ani_cluster_id>/
      members.tsv
      contig_filter.tsv          # retained/discarded UNI56-bearing contigs per genome
      filtered_contigs/          # marker-guided backbone contigs used for SNP alignment
      tree.nwk                   # cluster-size >= --snp_tree_min_cluster_size
      core_snps.fna              # cluster-size >= --snp_tree_min_cluster_size and variable sites present
```

Python output (`--save_dir`):

```text
<save_dir>/
  tree.nwk or tree_final.nwk
  tree_final.png                  # marker-selection mode
  marker_count_matrix.csv
  marker_selection_rf_values.txt  # marker-selection mode
  log_genomes_removed.txt
  ani/
    ani_pairwise.tsv
    ani_clusters.tsv
    ani_representatives.tsv
    ani_kept_genomes.txt
  snp_trees/                     # only with --snp yes
    snp_tree_summary.tsv
    <ani_cluster_id>/
      members.tsv
      contig_filter.tsv
      filtered_contigs/
      tree.nwk
      core_snps.fna              # only when variable SNP sites are found
  logfile_*.txt
  temp/
    *.zip
    itol/
```

## Repository Structure

```text
sgtree/
  sgtree/                 # Python package implementation
    benchmarks/           # synthetic benchmark generation/evaluation package
    ani.py                # ANI clustering + SNP-tree helpers
  sgtree.py               # backward-compatible wrapper
  main.nf                 # Nextflow entrypoint
  workflows/              # DSL2 workflow composition
  modules/                # DSL2 process modules
  bin/                    # helper scripts and launch wrappers
  tests/
    regression_parity.py  # cross-engine parity checks
  resources/
    models/               # combined marker-set HMM files
  testgenomes/            # example query/reference data
  runs/                   # runtime outputs/work/logs (.gitkeep tracked)
  pixi.toml               # reproducible environment + tasks
  nextflow.config         # runtime defaults and CPU settings
  docs/BENCHMARKS.md      # synthetic contamination benchmark design
```

## Benchmarking

Generate the default low-runtime synthetic benchmark:

```bash
pixi run benchmark-generate
```

Run the current systematic benchmark suite:

```bash
pixi run benchmark-run
```

The benchmark generator derives a clean truth panel automatically from the bundled Chloroflexi proteomes and `UNI56.hmm`, then creates duplicate, triplicate, and replacement scenarios against a fixed truth tree. Results are written under `runs/benchmarks/`. The canonical benchmark implementation now lives under `sgtree/benchmarks/`. See [docs/BENCHMARKS.md](/home/fschulz/dev/sgtree/docs/BENCHMARKS.md).

Burkholderiaceae benchmark assets:

- `pixi run benchmark-prepare-burkholderiaceae` materializes the requested 50-genome panel under [`testgenomes/Burkholderiaceae50`](/home/fschulz/dev/software/sgtree/testgenomes/Burkholderiaceae50).
- The same command writes taxonomy sidecars:
  - [`testgenomes/burkholderiaceae50.lookup`](/home/fschulz/dev/software/sgtree/testgenomes/burkholderiaceae50.lookup)
  - [`testgenomes/burkholderiaceae50_taxonomy.tsv`](/home/fschulz/dev/software/sgtree/testgenomes/burkholderiaceae50_taxonomy.tsv)
  - [`testgenomes/burkholderiaceae50_selection.tsv`](/home/fschulz/dev/software/sgtree/testgenomes/burkholderiaceae50_selection.tsv)
- The curated panel contains `20` genus/species buckets:
  - `18` singleton species from distinct genera
  - one `6`-strain species cluster
  - one `26`-strain species cluster
- Verified full run:

```bash
pixi run sgtree \
  --genomedir testgenomes/Burkholderiaceae50 \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/burkholderiaceae_ani_benchmark \
  --num_cpus 24 \
  --ani_cluster yes \
  --snp yes \
  --ani_threshold 95 \
  --marker_selection yes \
  --keep_intermediates yes
```

- Verified outputs for that run:
  - `20` ANI representatives in [`runs/burkholderiaceae_ani_benchmark/ani/ani_representatives.tsv`](/home/fschulz/dev/software/sgtree/runs/burkholderiaceae_ani_benchmark/ani/ani_representatives.tsv)
  - one `26`-member ANI cluster with a `524`-site SNP alignment after filtering to shared-core UNI56 backbone contigs
  - one `6`-member ANI cluster with no variable SNP sites, reported explicitly as `no_snp_sites`

## Workflow

```text
                            +-------------------+
                            |  Input Proteomes  |
                            |  + HMM Models     |
                            +---------+---------+
                                      |
                                      v
                             +--------+--------+
                             |    HMMSEARCH    |
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             | PARSE_HMMSEARCH |
                             | marker matrix   |
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             | EXTRACT_SEQS    |
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             | ALIGN (hmmalign/|
                             | mafft/linsi)    |
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             | ELIM_DUPLICATES |
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             |     TRIMAL      |
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             | BUILD_SUPERMATRIX|
                             +--------+--------+
                                      |
                                      v
                             +--------+--------+
                             |  TREE_BUILDER   |
                             |   tree.nwk      |
                             +--------+--------+
                                      |
                          marker_selection?
                           /            \
                        no               yes
                        |                 |
                        v                 v
                  +-----+-----+   +-------+--------+
                  | iTOL TXT  |   | per-marker     |
                  | marker_*  |   | TRIMAL+TREEBLD |
                  +-----------+   +-------+--------+
                                         |
                                         v
                                  +------+------+
                                  | RF_SELECTION|
                                  +------+------+
                                         |
                                 singles?|
                                  /      \
                               no         yes
                               |           |
                               v           v
                      +--------+---+   +---+--------+
                      | WRITE_CLEAN |   |REMOVE_     |
                      | ALIGNMENTS  |   |SINGLES     |
                      +--------+----+   +---+--------+
                               \           /
                                \         /
                                 v       v
                               +--+------+
                               |TRIMAL_FINAL
                               +--+------+
                                  |
                                  v
                             +----+-----+
                             |SUPERMATRIX|
                             +----+-----+
                                  |
                                  v
                             +----+-----+
                             |TREE_BUILDER|
                             |tree_final |
                             +----+-----+
                                  |
                                  v
                       +----------+-----------+
                       | tree_final.png       |
                       | marker_counts.txt    |
                       | marker_selection_rf  |
                       +----------------------+
```

## Repository Hygiene

Use this command for a clean runtime workspace between runs. This no longer deletes benchmark outputs:

```bash
pixi run clean-runtime
```

Use these commands for more targeted cleanup:

```bash
pixi run clean-regression
pixi run clean-benchmarks
pixi run clean-reference-cache
pixi run clean-all
```

## Authors and Contributors

| Author | Email | Date |
|---|---|---|
| Ewan Whittaker-Walker | ewanww@berkeley.edu | 05/19/2019 |
| Frederik Schulz | fschulz@lbl.gov | Since 2019 |
| Juan C. Villada | jvillada@lbl.gov | Since 2021 |
| Marianne Buscaglia | mbuscaglia@lbl.gov | Since 2022 |
