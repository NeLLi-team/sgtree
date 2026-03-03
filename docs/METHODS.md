# SGTree Methods and Implementation Documentation

This document describes the architecture and execution logic of `sgtree`, with emphasis on reproducibility, method choices, and exact software versions. The descriptions are grounded in a full instrumented run executed on March 3, 2026 (`runs/docs_audit/`), which completed with 231 successful tasks and no failures.

## Scope and Workflow

`sgtree` infers species trees from marker-gene proteins by combining profile-HMM search, per-marker sequence extraction, per-marker alignment, duplicate handling, alignment trimming, supermatrix construction, and phylogenetic inference. In the default configuration, pyhmmer is used for HMM search and profile-guided alignment, trimAl is used for gap-based column trimming, and FastTree is used for rapid approximate maximum-likelihood tree inference. Optional branches support reference-assisted duplicate resolution, RF-distance-based marker-copy selection, singleton filtering, and IQ-TREE-based phylogenetic inference.

## Repository and Execution Model

The primary implementation is a Nextflow DSL2 pipeline. The entrypoint is [`main.nf`](/home/fschulz/dev/sgtree/main.nf), workflow composition is defined in `workflows/`, and process modules are defined in `modules/`. Command-level operations are implemented in `bin/` scripts (for example, pyhmmer wrappers, parsing/filtering utilities, RF-based marker selection, and iTOL helpers). A Python package implementation is also maintained in `sgtree/` and can be invoked via `pixi run sgtree-python`, while the recommended production interface is `pixi run sgtree`, which launches Nextflow through [`bin/sgtree-nextflow.sh`](/home/fschulz/dev/sgtree/bin/sgtree-nextflow.sh).

Environment reproducibility is handled through `pixi.toml` and `pixi.lock`. Pipeline defaults are declared in [`nextflow.config`](/home/fschulz/dev/sgtree/nextflow.config), including `aln=hmmalign`, `tree_method=fasttree`, `marker_selection=false`, and `singles=false`.

## Default Branch: Process Logic

In the default branch, `CONCAT_INPUTS` merges marker models and proteomes and creates per-marker HMM files. `HMMSEARCH` runs `run_hmmsearch_pyhmmer.py` and emits HMMER-compatible domain-table output, with either model-native bit cutoffs (`cut_ga`, `cut_tc`, `cut_nc`) or an explicit E-value threshold. `PARSE_HMMSEARCH` then applies genome-level and duplication filters, writes the marker count matrix, and generates per-marker ID lists together with combined FASTA/index files.

`EXTRACT_SEQUENCES` retrieves the selected protein sequences for each marker. Alignment proceeds according to `--aln`; in the verified run, `--aln hmmalign` was used, so `ALIGN_HMMALIGN` performed profile-guided marker alignment with pyhmmer. `ELIMINATE_DUPLICATES` removes surplus copies based on the score table from parsing. `TRIMAL` applies `-gt 0.1` to remove high-gap columns. `BUILD_SUPERMATRIX` concatenates trimmed marker alignments, and `FASTTREE` infers `tree.nwk`. `WRITE_HEATMAP` produces the iTOL-compatible marker-count annotation file.

## Optional Marker Selection and Singleton Filtering

When `--marker_selection true`, `MARKER_SELECTION_WF` builds per-marker trees and applies `rf_marker_selection.py` to resolve duplicated marker copies against the current species tree using normalized Robinson-Foulds distances. For each duplicated genome-marker case, the retained copy is the one minimizing RF distance; removed copies are recorded in `marker_selection_rf_values.txt`.

When `--singles true`, `remove_singles.py` applies a neighborhood-consistency filter. For each leaf, nearest-neighbor ordering in the marker tree is compared against the pruned species tree. The neighborhood size is adapted from global RF discordance unless `--num_nei` is set. Leaves exceeding the scoring cutoff are pruned before cleaned alignments are rebuilt and the final concatenated tree (`tree_final.nwk`) is inferred.

## Software Versions (Verified March 3, 2026)

| Tool | Version | Evidence Command | Function in Pipeline |
|---|---:|---|---|
| Nextflow | 25.10.4 | `./nextflow -version` | Workflow execution and provenance |
| Python | 3.12.12 | `pixi run python -V` | Runtime for package/scripts |
| pyhmmer | 0.12.0 | `pixi run python -c ...` | Profile-HMM search and profile alignment |
| MAFFT | 7.526 | `pixi run mafft --version` | Optional MSA mode |
| trimAl | 1.5.rev1 | `pixi run trimal --version` | Alignment trimming |
| FastTree | 2.2.0 | `pixi run FastTree` | Default tree inference |
| IQ-TREE | 3.0.1 | `pixi run iqtree --version` | Optional ML tree inference |
| Biopython | 1.86 | `pixi run python -c ...` | FASTA parsing/indexing utilities |
| ETE3 | 3.1.3 | `pixi run python -c ...` | Tree operations and RF calculations |
| pandas | 2.2.3 | `pixi run python -c ...` | Tabular filtering and scoring |

## Methodological Rationale

This toolchain favors transparent and reproducible execution while preserving practical throughput for marker-rich datasets. Nextflow provides explicit process boundaries and run-level provenance [1]. pyhmmer allows profile-HMM operations to remain in Python-native workflow steps while preserving compatibility with HMMER-style methods [7,9]. MAFFT is retained as a non-profile alignment option [2]. trimAl provides automated removal of weakly supported alignment regions prior to concatenation [3]. FastTree is used as a high-throughput default [4], whereas IQ-TREE is available when a model-intensive maximum-likelihood workflow is preferred [5,6]. ETE3 underlies the RF-based tree-comparison logic used in marker selection and singleton filtering [8].

## Reproducibility Record

The verification run used the command below and generated full trace/report/timeline/dag outputs in `runs/docs_audit/`.

```bash
nextflow -log runs/docs_audit/nextflow.log run main.nf \
  --genomedir testgenomes/Chloroflexi \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/docs_audit/out \
  --hmmsearch_cpus 4 \
  --align_cpus 4 \
  --fasttree_cpus 1 \
  -with-trace runs/docs_audit/trace.tsv \
  -with-report runs/docs_audit/report.html \
  -with-timeline runs/docs_audit/timeline.html \
  -with-dag runs/docs_audit/dag.html
```

This run completed as `small_stonebraker` with `succeededCount=231` and `failedCount=0`, and produced `tree.nwk`, `marker_count_matrix.csv`, `marker_count.txt`, and `color.txt` in `runs/docs_audit/out/`. The machine-readable record for this run is [`docs/run_manifest.yaml`](/home/fschulz/dev/sgtree/docs/run_manifest.yaml).

## Reference Validation

All references listed below were validated against Crossref DOI metadata on March 3, 2026. DOI resolution, journal, year, and title-level consistency were checked for each entry.

## References

1. Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nature Biotechnology*. 2017;35(4):316-319. https://doi.org/10.1038/nbt.3820
2. Katoh K, Standley DM. MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. *Molecular Biology and Evolution*. 2013;30(4):772-780. https://doi.org/10.1093/molbev/mst010
3. Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*. 2009;25(15):1972-1973. https://doi.org/10.1093/bioinformatics/btp348
4. Price MN, Dehal PS, Arkin AP. FastTree 2 - Approximately Maximum-Likelihood Trees for Large Alignments. *PLoS ONE*. 2010;5(3):e9490. https://doi.org/10.1371/journal.pone.0009490
5. Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ. IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. *Molecular Biology and Evolution*. 2014;32(1):268-274. https://doi.org/10.1093/molbev/msu300
6. Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular Biology and Evolution*. 2020;37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015
7. Eddy SR. Accelerated Profile HMM Searches. *PLoS Computational Biology*. 2011;7(10):e1002195. https://doi.org/10.1371/journal.pcbi.1002195
8. Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. *Molecular Biology and Evolution*. 2016;33(6):1635-1638. https://doi.org/10.1093/molbev/msw046
9. Larralde M, Zeller G. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics*. 2023;39(5):btad214. https://doi.org/10.1093/bioinformatics/btad214
10. Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422-1423. https://doi.org/10.1093/bioinformatics/btp163
