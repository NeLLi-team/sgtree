# SGTree Methods and Implementation Documentation

This document describes the architecture and execution logic of `sgtree`, with emphasis on reproducibility, method choices, and exact software versions. The descriptions are grounded in a full instrumented run executed on March 3, 2026 (`runs/docs_audit/`), which completed with 231 successful tasks and no failures.

## Scope and Workflow

`sgtree` infers species trees from marker-gene proteins by combining profile-HMM search, per-marker sequence extraction, per-marker alignment, duplicate handling, alignment trimming, supermatrix construction, and phylogenetic inference. In the default configuration, pyhmmer is used for HMM search and profile-guided alignment, trimAl is used for gap-based column trimming, and FastTree is used for rapid approximate maximum-likelihood tree inference. Optional branches support reference-assisted duplicate resolution, RF-distance-based marker-copy selection, singleton filtering, ANI-based genome clustering, and cluster-specific SNP phylogenies.

A workflow schematic corresponding to the current implementation is provided in [`docs/figures/workflow_overview.svg`](/home/fschulz/dev/sgtree/docs/figures/workflow_overview.svg) and [`docs/figures/workflow_overview.png`](/home/fschulz/dev/sgtree/docs/figures/workflow_overview.png).

## Repository and Execution Model

The repository contains both a Nextflow DSL2 implementation and a Python implementation. The entrypoint is [`main.nf`](/home/fschulz/dev/sgtree/main.nf), workflow composition is defined in `workflows/`, and process modules are defined in `modules/`. Command-level operations are implemented in `bin/` scripts. A Python package implementation is maintained in `sgtree/` and can be invoked directly via `pixi run sgtree-python`. The recommended user-facing entrypoint remains `pixi run sgtree`, but the wrapper in [`bin/sgtree-nextflow.sh`](/home/fschulz/dev/sgtree/bin/sgtree-nextflow.sh) now dispatches nucleotide (`fna`) and ANI/SNP-enabled runs to the Python engine so that gene calling, ANI clustering, and SNP-tree generation all execute in the same runtime.

Environment reproducibility is handled through `pixi.toml` and `pixi.lock`. Pipeline defaults are declared in [`nextflow.config`](/home/fschulz/dev/sgtree/nextflow.config), including `aln=hmmalign`, `tree_method=fasttree`, `marker_selection=false`, `ani_cluster=false`, `snp=false`, and `singles=false`.

## Default Branch: Process Logic

In the default branch, `CONCAT_INPUTS` merges marker models and normalizes proteomes before concatenation. Input headers are rewritten to `genome|protein`, malformed sequence/header joins are repaired, and non-IUPAC amino-acid symbols are converted to `X` (with `*` removed) so downstream pyhmmer readers are resilient to formatting drift in source files. `HMMSEARCH` runs `run_hmmsearch_pyhmmer.py` and emits HMMER-compatible domain-table output, with either model-native bit cutoffs (`cut_ga`, `cut_tc`, `cut_nc`) or an explicit E-value threshold. `PARSE_HMMSEARCH` then applies genome-level and duplication filters, writes the marker count matrix, and generates per-marker ID lists together with combined FASTA/index files.

`EXTRACT_SEQUENCES` retrieves the selected protein sequences for each marker. Alignment proceeds according to `--aln`; in the verified run, `--aln hmmalign` was used, so `ALIGN_HMMALIGN` performed profile-guided marker alignment with pyhmmer. `ELIMINATE_DUPLICATES` removes surplus copies based on the score table from parsing. `TRIMAL` applies `-gt 0.1` to remove high-gap columns. `BUILD_SUPERMATRIX` concatenates trimmed marker alignments, and `FASTTREE` infers `tree.nwk`. `WRITE_HEATMAP` produces the iTOL-compatible marker-count annotation file.

When nucleotide assemblies are supplied, SGTree first gene-calls them with `pyrodigal` and writes both `gene_calls.tsv` and `genome_manifest.tsv`. If `--ani_cluster true` is set, all query and reference assemblies are compared pairwise with `skani`, edges at or above the selected ANI cutoff are retained, and the resulting graph is clustered with `mcl`. One representative assembly is retained per ANI cluster for species-tree inference, while cluster membership and representative choices are recorded under `ani/`. On the Burkholderiaceae 50-genome runtime benchmark, this reduced the working set to twenty representative genomes, consisting of one 26-member species cluster, one 6-member species cluster, and eighteen singletons.

## Optional Marker Selection and Singleton Filtering

When `--marker_selection true`, `MARKER_SELECTION_WF` builds per-marker trees and applies `rf_marker_selection.py` to resolve duplicated marker copies against the current species tree using normalized Robinson-Foulds distances. For each duplicated genome-marker case, the retained copy is the one minimizing RF distance; removed copies are recorded in `marker_selection_rf_values.txt`.

When `--singles true`, `remove_singles.py` applies a neighborhood-consistency filter. For each leaf, nearest-neighbor ordering in the marker tree is compared against the pruned species tree. The neighborhood size is adapted from global RF discordance unless `--num_nei` is set. In runs without duplicate-resolution activity, highly discordant singleton markers now escalate from the plain neighbor heuristic to an RF-aware backbone score so that obvious replacement-style outliers are not masked by the contaminated guide tree. Leaves exceeding the resulting cutoff are pruned before cleaned alignments are rebuilt and the final concatenated tree (`tree_final.nwk`) is inferred.

When `--snp true` is enabled in combination with ANI clustering, SGTree builds an additional cluster-level SNP phylogeny for each ANI cluster meeting the minimum size threshold. This step is deliberately more conservative than the species-tree path. Before whole-genome SNP alignment, each cluster is gene-called and searched again against the active UNI56 marker set. Only contigs carrying cluster-core UNI56 markers are eligible for SNP analysis, and contigs are discarded if they contain UNI56 markers that are not shared across the cluster. The remaining contigs must also align back to the representative backbone at or above 95% ANI before being retained. SGTree writes `contig_filter.tsv` and `filtered_contigs/` for each cluster so the retained backbone is auditable. SNP alignments are then built from the filtered contigs only. In the verified Burkholderiaceae run, the 26-member cluster retained one shared backbone contig per genome and yielded a 524-site SNP alignment, whereas the 6-member cluster retained six backbone contigs but contained no variable SNP sites and was therefore reported as `no_snp_sites`.

## Benchmark Frameworks

The synthetic contamination benchmark framework is implemented in [`sgtree/benchmarks/__init__.py`](/home/fschulz/dev/sgtree/sgtree/benchmarks/__init__.py), with [`sgtree/benchmark.py`](/home/fschulz/dev/sgtree/sgtree/benchmark.py) retained as a compatibility import alias. Benchmark generation first runs SGTree on a clean source proteome collection to identify a candidate truth panel, ranks marker genes by single-copy prevalence and score stability, extracts a smaller truth-marker HMM subset, and rebuilds a clean truth tree from the selected genomes and markers. Synthetic contamination is then injected into that fixed truth panel.

Three contamination scenarios are currently used. `duplicate_only` introduces only added contaminant copies of existing markers. `replacement_only` introduces only singleton-replacement contamination, in which a contaminant copy replaces the native single-copy marker. `combined` contains both added-copy and singleton-replacement contamination in the same panel. Each event is recorded explicitly in `events.tsv`, including recipient genome, donor genome, marker identity, source relation (`within_group` or `cross_group`), and the expected removal or retention status for evaluation. Each scenario also writes `genome_summary.tsv`, which records how many wrong markers were introduced per genome.

Cross-clade benchmarks are generated from two prebuilt benchmark families and are restricted to contamination exchange between the two group labels (`flavo` and `gamma`). This produces a more divergent contamination regime than the within-clade panels. The current benchmark notebook and summary tables for the six benchmark families are documented in [`docs/BENCHMARKS.md`](/home/fschulz/dev/sgtree/docs/BENCHMARKS.md).

Benchmark evaluation now uses one intended cleanup mode per scenario. `duplicate_only` is evaluated with duplicate cleanup only (`marker_selection=yes`, `singles=no`). `replacement_only` is evaluated with `singles_neighbor` activated (`marker_selection=yes`, `singles=yes`, `singles_mode=neighbor`). `combined` is evaluated with both duplicate cleanup and `singles_neighbor` active. Each dataset/scenario run records the initial contaminated-tree RF, the cleaned-tree RF, the RF delta, and the number of contaminant markers removed out of the number added. Exported benchmark tables are written to [`docs/data/benchmark_dataset_overview.tsv`](/home/fschulz/dev/sgtree/docs/data/benchmark_dataset_overview.tsv), [`docs/data/benchmark_genome_contamination.tsv`](/home/fschulz/dev/sgtree/docs/data/benchmark_genome_contamination.tsv), and [`docs/data/benchmark_summary_all.tsv`](/home/fschulz/dev/sgtree/docs/data/benchmark_summary_all.tsv).

A separate real-genome runtime benchmark now exists for the ANI/SNP workflow. `pixi run benchmark-prepare-burkholderiaceae` materializes a 50-genome Burkholderiaceae panel in which two species are represented by multi-strain clusters of size 26 and 6, and the remaining eighteen species are singletons from distinct genera. This panel is not part of the contamination benchmark catalog under `runs/benchmarks`; instead it is intended to measure ANI collapse, representative selection, and cluster-level SNP resolution on real assemblies.

## Implementation Validation

The current codebase has been validated in two ways. First, the repository unit suite completed successfully for the ANI/SNP additions with `pixi run python -m unittest tests.test_ani tests.test_ani_clustering tests.test_cli`, and the broader SGTree regression-facing unit slice also passed with `pixi run python -m unittest tests.test_benchmark tests.test_input_stage tests.test_marker_selection tests.test_benchmark_dataset`. Second, the Burkholderiaceae ANI/SNP runtime benchmark completed successfully through the public wrapper using `pixi run sgtree --genomedir testgenomes/Burkholderiaceae50 --modeldir resources/models/UNI56.hmm --outdir runs/burkholderiaceae_snp_backbone_20260313_fix --num_cpus 12 --ani_cluster yes --snp yes --ani_threshold 95 --marker_selection yes --keep_intermediates yes`.

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
