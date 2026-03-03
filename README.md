# SGTree

SGTree builds species trees from conserved marker proteins in genome proteomes.  
It now supports both:

- a modular Python package (`python -m sgtree`)
- a modular Nextflow DSL2 workflow (`main.nf`, `modules/`, `workflows/`)

## Status

The codebase and outputs are functionally aligned across Python and Nextflow, including:

- basic run (`tree.nwk`)
- marker-selection run (`tree_final.nwk`)
- singleton filtering (`--singles yes`)

Regression parity is enforced by:

```bash
pixi run test-regression
```

This test checks topology parity, marker matrix parity, and RF selection status parity.
Outputs are written to `runs/regression/`.

## Setup

Recommended: Pixi (reproducible and already configured in this repo).

```bash
pixi install
```

Alternative environments are available in:

- `environment.yml`
- `linux_env.txt`
- `osx_env.txt`

## Run

### Python package path

```bash
# basic
pixi run python -m sgtree testgenomes/Chloroflexi hmms/UNI56 --num_cpus 8

# marker selection + references
pixi run python -m sgtree \
  testgenomes/Chloroflexi hmms/UNI56 \
  --num_cpus 8 \
  --marker_selection yes \
  --ref testgenomes/chlorref \
  --singles yes \
  --save_dir runs/python/manual_full
```

Backward-compatible wrapper still works:

```bash
pixi run ./sgtree.py testgenomes/Chloroflexi hmms/UNI56 --num_cpus 8
```

### Nextflow path

```bash
# basic
pixi run ./nextflow -log runs/nextflow/logs/manual_basic.log run main.nf \
  --genomedir testgenomes/Chloroflexi \
  --modeldir hmms/UNI56 \
  --outdir runs/nextflow/manual_basic

# marker selection + references + singles
pixi run ./nextflow -log runs/nextflow/logs/manual_full.log run main.nf \
  --genomedir testgenomes/Chloroflexi \
  --modeldir hmms/UNI56 \
  --outdir runs/nextflow/manual_full \
  --marker_selection true \
  --ref testgenomes/chlorref \
  --singles yes
```

Useful Nextflow resource controls:

- `--hmmsearch_cpus`
- `--align_cpus`
- `--fasttree_cpus`
- `--render_png`

## Input Requirements

Proteomes must be FASTA (`*.faa`) with headers like:

```text
>IMG2684622718|2685462912
MLCAFAEEEAKIAETVGKVATELKVKKLLSDFATKEGEEHISTYNKIAMTAKAEGYADIEAMLCAFAEEEAKLQKL
```

The part before `|` must match the genome filename stem.

## Output Structure

### Nextflow output (`--outdir`)

Basic run:

```text
<outdir>/
  tree.nwk
  marker_count_matrix.csv
  marker_count.txt
  color.txt
  log_genomes_removed.txt
```

Marker-selection run:

```text
<outdir>/
  tree.nwk
  tree_final.nwk
  tree_final.png
  marker_count_matrix.csv
  marker_counts.txt
  marker_selection_rf_values.txt
  color.txt
  log_genomes_removed.txt
```

### Python output (`--save_dir`)

Python keeps final files at top level and archives intermediates under `temp/`:

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
  sgtree.py               # backward-compatible thin wrapper
  main.nf                 # Nextflow entrypoint
  workflows/              # DSL2 workflow composition
  modules/                # DSL2 process modules
  bin/                    # helper scripts used by Nextflow modules
  tests/
    regression_parity.py  # cross-engine parity checks
  hmms/                   # marker model sets
  testgenomes/            # example query/reference data
  runs/                   # runtime outputs/work/logs (.gitkeep tracked)
  pixi.toml               # reproducible environment + tasks
  nextflow.config         # runtime defaults and CPU settings
```

## ASCII Workflow Diagram

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
                             | BUILD_SUPERMAT  |
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
                             |SUPERMAT  |
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

Runtime artifacts are expected under `runs/` by default.  
Large runtime artifacts (`work/`, `.nextflow*`, `runs/`, `results_*`) can grow quickly and are usually not versioned.  
For strict cleanup between runs:

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
