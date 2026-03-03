# SGTree

SGTree is an end-to-end workflow for phylogenetic tree building. Use the provided sets of HMMs or provide your own HMMs to find the proteins of interest. SGTree then performs gene tree to approximate species tree reconciliation to select the most likely correct copy of a protein in case of duplications (paralogs, contamination). 

## Setup

Install the Pixi environment:

```bash
pixi install
```

The environment is managed through `pixi.toml` only.

## Run

Primary interface (Nextflow):

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

`pixi run sgtree` writes logs automatically to `runs/nextflow/logs/`.
Marker searches and `--aln hmmalign` are run with `pyhmmer` (HMMER-compatible search output).

Example with IQ-TREE and explicit HMM threshold mode:

```bash
pixi run sgtree \
  --genomedir testgenomes/Chloroflexi \
  --modeldir resources/models/UNI56.hmm \
  --tree_method iqtree \
  --iqtree_fast true \
  --hmmsearch_cutoff cut_ga
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
- Typical inclusion presets:
- Balanced: `--percent_models 10 --max_sdup 2 --max_dupl 0.25`
- Strict: `--percent_models 30 --max_sdup 1 --max_dupl 0.10`
- Relaxed: `--percent_models 5 --max_sdup -1 --max_dupl -1`

## Input Requirements

Proteomes must be FASTA (`*.faa`) with headers like:

```text
>IMG2684622718|2685462912
MLCAFAEEEAKIAETVGKVATELKVKKLLSDFATKEGEEHISTYNKIAMTAKAEGYADIEAMLCAFAEEEAKLQKL
```

The text before `|` must match the genome filename stem.

## Output Structure

Nextflow output (`--outdir`):

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
```

Python output (`--save_dir`):

```text
<save_dir>/
  tree.nwk or tree_final.nwk
  tree_final.png                  # marker-selection mode
  marker_count_matrix.csv
  marker_selection_rf_values.txt  # marker-selection mode
  log_genomes_removed.txt
  logfile_*.txt
  temp/
    *.zip
    itol/
```

## Repository Structure

```text
sgtree/
  sgtree/                 # Python package implementation
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
```

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

Use this command for a clean runtime workspace between runs:

```bash
pixi run clean-runtime
```

## Authors and Contributors

| Author | Email | Date |
|---|---|---|
| Ewan Whittaker-Walker | ewanww@berkeley.edu | 05/19/2019 |
| Frederik Schulz | fschulz@lbl.gov | Since 2019 |
| Juan C. Villada | jvillada@lbl.gov | Since 2021 |
| Marianne Buscaglia | mbuscaglia@lbl.gov | Since 2022 |
