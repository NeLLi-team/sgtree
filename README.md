# SGTree

SGTree is an end-to-end workflow for phylogenetic tree building. Use the provided sets of HMMs or provide your own HMMs to find the proteins of interest. SGTree then performs gene tree to approximate species tree reconciliation to select the most likely correct copy of a protein in case of duplications (paralogs, contamination). 

## Setup

Install the Pixi environment:

```bash
pixi install
```

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
  --modeldir models/UNI56.hmm
```

Marker-selection run with references and singleton filtering:

```bash
pixi run sgtree \
  --genomedir testgenomes/Chloroflexi \
  --modeldir models/UNI56.hmm \
  --outdir runs/nextflow/manual_full \
  --marker_selection true \
  --ref testgenomes/chlorref \
  --singles yes
```

`pixi run sgtree` writes logs automatically to `runs/nextflow/logs/`.
Marker searches are run with `pyhmmer` (HMMER-compatible domain search output).

Second choice (Python implementation without nextflow):

```bash
pixi run sgtree-python testgenomes/Chloroflexi models/UNI56.hmm --num_cpus 8
```

Backward-compatible wrapper:

```bash
pixi run ./sgtree.py testgenomes/Chloroflexi models/UNI56.hmm --num_cpus 8
```

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
  models/                 # combined marker-set HMM files
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
                             |   FASTTREE      |
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
                  | marker_*  |   | TRIMAL+FASTTREE|
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
                             |FASTTREE   |
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
