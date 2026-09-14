# SGTree

SGTree builds a species tree from genomes and a marker-set HMM file. It finds the markers,
aligns them, concatenates them, and infers the tree.

Input is a directory of protein FASTA files or genome assemblies. Output is a Newick tree
plus the tables that show how it was built. An optional marker-selection mode builds one
tree per marker and resolves duplicate copies against the species tree. Experimental
singleton modes can remove discordant marker copies; `loo_profile` reports candidates
without changing the marker set.

## Requirements

- Linux or macOS
- [Pixi](https://pixi.sh) for the environment

All other tools (HMMER through `pyhmmer`, MAFFT, FAMSA, trimAl, VeryFastTree, IQ-TREE,
skani, minimap2, MCL, Pyrodigal) come from the Pixi environment.

## Install

```bash
git clone https://github.com/NeLLi-team/sgtree.git
cd sgtree
pixi install
```

## Quick start

Run the bundled 10-genome example:

```bash
pixi run example
```

This writes `runs/example_basic/tree.nwk`. A successful run ends with a line like:

```text
Final tree: /path/to/sgtree/runs/example_basic/tree.nwk
```

Confirm that the tree exists and is not empty:

```bash
test -s runs/example_basic/tree.nwk
```

Then run your own data:

```bash
pixi run sgtree \
  --genomedir /path/to/genomes \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/my_first_tree
```

## Input

SGTree accepts protein FASTA (`.faa`) or genome assembly FASTA (`.fna`, `.fa`, `.fasta`).

- **Directory input** (usual case): one genome per file. The genome ID comes from the file
  name without its extension.
- **Single protein FASTA input**: one file can hold all genomes. Headers must already be
  in `genome|protein` form. SGTree keeps the genome part.
- **Single assembly FASTA input**: one genome, with one record per contig. To include
  several assemblies, put each genome in its own file.

Keep only genome FASTA files in the input directory, and use one sequence format throughout.
Hidden files are ignored. SGTree rejects empty files, files without FASTA headers, and
mixed protein/assembly input.
For other filename extensions, it infers the format from the first sequence line; use the
extensions above when the format is known. For directory input, genome filenames must
remain distinct after normalization. Contig IDs must also be unique within each assembly.

Assemblies are gene-called with Pyrodigal before the marker search. The gene coordinates
go to `gene_calls.tsv`.

SGTree rewrites every header to a three-field ID:

```text
>GCA_000012365|GCA_000012365_CP000086.1|gene_000001
MLCAFAEEEAKIAETVGKVATELKVKKLLSDFATKEGEEHISTYNKIAMTAKAEGYADIEA
```

The fields are genome, contig, and gene. Every header change is recorded in
`proteomes_header_map.tsv`, with one row per sequence and the columns
`source_file`, `original_header`, `normalized_header`, `genome_id`, `contig_id`,
`gene_id`, and `contig_inference`. The genome list goes to `genome_manifest.tsv`.

During normalization SGTree replaces invalid amino-acid characters with `X`, removes `*`,
repairs malformed header joins, and removes `|` from the source tokens.

### Marker sets

`resources/models/` holds the bundled sets. Use `UNI56.hmm` for bacteria and archaea.
Each file is one concatenated HMM file. You can supply your own.

| File | Markers | Scope |
|---|---|---|
| `UNI56.hmm` | 56 | universal single-copy, bacteria and archaea |
| `gtdbbac.hmm` | 120 | GTDB bacterial set |
| `gtdbarc.hmm` | 122 | GTDB archaeal set |
| `RProt16.hmm` | 16 | ribosomal proteins |
| `GVOG7.hmm`, `GVOG9.hmm` | 7, 9 | giant viruses |
| `mitomarkers108.hmm` | 108 | mitochondrial markers |
| `COX123.hmm`, `mcrABG.hmm`, `rnapol.hmm` | 3, 3, 3 | small targeted sets |
| `UNI56r.hmm`, `UNI56nr.hmm` | 30, 26 | UNI56 split into ribosomal and non-ribosomal |

## Common runs

Boolean options accept `yes`/`no`, `true`/`false`, or `1`/`0`.

**Basic tree, 24 threads:**

```bash
pixi run sgtree \
  --genomedir testgenomes/example_fna10 \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/basic \
  --num_cpus 24
```

**Marker selection with a report-only singleton diagnostic.** This builds per-marker
trees, resolves duplicate copies against the species tree, and writes candidate evidence
without removing single-copy markers:

```bash
pixi run sgtree \
  --genomedir testgenomes/example_fna10 \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/marker_selection \
  --marker_selection yes \
  --singles yes \
  --singles-mode loo_profile
```

**Species tree with IQ-TREE:**

```bash
pixi run sgtree \
  --genomedir testgenomes/example_fna10 \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/iqtree \
  --tree_method iqtree \
  --iqtree_fast yes
```

**Assemblies with ANI clustering and per-cluster SNP trees.** ANI clustering keeps one
representative genome per cluster for the species tree:

```bash
pixi run sgtree \
  --genomedir /path/to/assemblies \
  --modeldir resources/models/UNI56.hmm \
  --outdir runs/ani \
  --ani_cluster yes \
  --snp yes \
  --ani_threshold 95
```

**Add reference genomes.** References are placed in the tree with your genomes. SGTree
caches their proteins and marker-search results, then combines them with each query run:

```bash
pixi run sgtree \
  --genomedir /path/to/queries \
  --modeldir resources/models/UNI56.hmm \
  --ref /path/to/references \
  --outdir runs/with_refs
```

The positional form `pixi run sgtree <genomedir> <modeldir>` also works.

## Options

### Input and output

| Option | Default | Description |
|---|---|---|
| `--genomedir` | required | Directory of `.faa` or `.fna` files, or one concatenated FASTA |
| `--modeldir` | required | Marker-set `.hmm` file |
| `--outdir`, `--save_dir` | `runs/python/SG_<input>_<ref>_<model>_<timestamp>` | Output directory |
| `--ref` | none | Directory of reference genomes |
| `--ref_concat` | `runs/reference_cache` | Where reference proteins and marker results are cached |
| `--num_cpus` | `8` | Threads |
| `--keep_intermediates` | `no` | Keep intermediate directories instead of archiving them |

### Marker search

| Option | Default | Description |
|---|---|---|
| `--hmmsearch_cutoff` | `cut_ga` | `cut_ga`, `cut_tc`, `cut_nc`, or `evalue` |
| `--hmmsearch_evalue` | `1e-5` | E-value, used with `--hmmsearch_cutoff evalue` |
| `--percent_models` | `10` | Minimum percentage of markers a genome must carry |
| `--max_sdup` | `-1` | Maximum copies of one marker per genome; `-1` disables |
| `--max_dupl` | `-1` | Maximum fraction of markers in more than one copy; `-1` disables |
| `--lflt` | `0` | Drop hits shorter than this percentage of the median hit length |

Use `cut_ga` with curated sets such as UNI56. For custom sets start with
`--hmmsearch_cutoff evalue --hmmsearch_evalue 1e-5`, then tighten if false hits appear.

### Alignment and tree

| Option | Default | Description |
|---|---|---|
| `--aln` | `hmmalign` | `hmmalign`, `mafft`, `mafft-linsi`, or `famsa` |
| `--tree_method` | `veryfasttree` | `veryfasttree`, `iqtree`, or the alias `fasttree` |
| `--iqtree_fast` | `yes` | Add `-fast` when `--tree_method iqtree` |
| `--iqtree_model` | `LG+F+I+G4` | IQ-TREE model string |

`hmmalign` aligns each marker against its profile HMM. `mafft`, `mafft-linsi`, and `famsa`
align marker sequences de novo. Their speed and accuracy depend on the marker set and
dataset; SGTree does not benchmark them for a universal ranking.

### Marker selection

| Option | Default | Description |
|---|---|---|
| `--marker_selection` | `no` | Build per-marker trees and resolve duplicates against the species tree |
| `--selection_mode` | `coordinate` | `coordinate` or `legacy` duplicate resolution |
| `--selection_max_rounds` | `5` | Maximum coordinate-descent rounds |
| `--selection_global_rounds` | `1` | Maximum guide-tree rebuild rounds |
| `--lock_references` | `no` | Keep reference duplicate choices fixed by score |
| `--singles` | `no` | Enable singleton analysis; behavior depends on `--singles-mode` |
| `--singles-mode` | `delta_rf` | Detector; see the table below |
| `--singles_min_rfdist` | `0.25` | Minimum marker-tree to species-tree RF distance before singleton removal starts. `neighbor_clade`, `neighbor_ml`, and `gcp` do not use this gate |
| `--num_nei` | `0` | Neighborhood size; `0` selects it automatically |

### ANI clustering and SNP trees

| Option | Default | Description |
|---|---|---|
| `--ani_cluster` | `no` | Cluster genomes by ANI and keep one representative per cluster |
| `--ani_threshold` | `95` | ANI cutoff for graph edges |
| `--ani_backend` | `auto` | `auto`, `skani`, or `minimap2` |
| `--ani_mcl_inflation` | `2.0` | MCL inflation |
| `--snp` | `no` | Build a SNP tree per ANI cluster; needs `--ani_cluster yes` |
| `--snp_tree_min_cluster_size` | `3` | Minimum cluster size for a SNP tree |

Before SNP alignment SGTree keeps only the contigs that carry shared cluster-core markers
and that align back to the representative genome at 95% ANI or better.

For each genome, positions covered by more than one primary contig alignment are excluded
from the core SNP alignment. This avoids arbitrary haploid allele calls and can reduce the
number of usable sites. A cluster with no variable sites receives an unresolved star tree.

## Output

A basic run writes:

```text
<outdir>/
  tree.nwk                       species tree
  marker_count_matrix.csv        markers found per genome
  log_genomes_removed.txt        genomes dropped by the inclusion filters
  genome_manifest.tsv            genome list with source paths
  gene_calls.tsv                 gene coordinates, assembly input only
  proteomes_header_map.tsv       original header to normalized ID
  logfile_<timestamp>.txt        run log
  temp/                          archived intermediates
```

A marker-selection run writes `tree_final.nwk` as the final tree and
`marker_selection_rf_values.txt`. With more than one rebuild round it also keeps a copy
of each round as `tree_round_<N>.nwk`. With `--singles yes` it writes
`singleton_candidates.tsv`, one row per scored marker copy. The first-pass `tree.nwk`
moves into `temp/`.

Use `tree.nwk` as the result of a basic run. Use `tree_final.nwk` as the result of a
marker-selection run. In the latter case, the archived `tree.nwk` is the guide tree built
before marker selection.

`--ani_cluster yes` adds an `ani/` directory with the pairwise ANI table, the clusters,
the representatives, and the ANI graph. `--snp yes` adds
`snp_trees/snp_tree_summary.tsv`, which records the outcome for every cluster. Eligible
clusters have a subdirectory with their members and contig-filter results. A completed
SNP analysis adds `core_snps.fna` and `tree.nwk`; a cluster without variable sites gets an
unresolved star tree. Clusters below the size threshold have no subdirectory.

At the end of a run SGTree archives its intermediate directories as zip files under
`temp/`. Smaller intermediate files are compressed where they are and keep their original
names. Pass `--keep_intermediates yes` to keep everything readable for debugging.
Top-level entries whose names do not match SGTree outputs are left in place. Keep unrelated
work in a separate directory.

## Singleton diagnostics and experimental pruning

A marker copy can disagree with the other marker trees because of contamination,
horizontal transfer, hidden paralogy, taxon error, alignment error, or tree-estimation
error. The singleton modes detect topological discordance; they do not determine its
cause.

Marker selection builds one tree per marker. Where a genome has more than one copy of a
marker, SGTree keeps the copy that makes the marker tree agree best with the species tree,
measured by Robinson-Foulds distance. `--selection_mode coordinate` revisits that choice
over several rounds; `legacy` decides once.

With `--singles yes` SGTree also evaluates single-copy markers. All modes except
`loo_profile` can remove proposed leaves. Those automatic modes are experimental and have
not been validated for field contamination cleanup. Across all markers, the budget admits
at most one proposed removal per genome and retains at least one marker in that genome. For each
marker, the proposed leaves are removed together only if their removal strictly lowers
the normalized Robinson-Foulds distance to the species tree. A rejected proposal does not
free its genome's budget for a second candidate.

| `--singles-mode` | Behavior |
|---|---|
| `delta_rf` | Removes the leaf whose removal most improves marker-to-species RF. |
| `composite` | Requires RF gain, local topology mismatch, branch-length signal, and bitscore outlier signal to agree. |
| `recipient_consensus` | Scores the copy against the genome's nearest species-tree neighborhood. |
| `contig_consensus` | Starts from `composite`, then asks whether other markers on the same contig disagree. Needs reliable contig IDs. |
| `neighbor_clade` | Asks whether a copy sits where the genome's closest neighbors put theirs, without a whole-tree RF trigger. |
| `gcp` | Genome Consistency Profiling: per-genome z-scores over marker features, combined with IsolationForest and HDBSCAN. Flags at most one marker per genome. |
| `loo_profile` | Reports only. Compares each placement against the genome's other marker trees and writes the evidence. Removes nothing. |

`topoknn`, `hybrid`, and `neighbor_ml` are baselines for benchmark comparison. Do not use
them for analysis.

### What the evidence columns mean

`--singles-mode loo_profile` writes `singleton_candidates.tsv` and prunes nothing. Each row
carries the abstention reason, a MAD-scaled conflict score (`loo_robust_z`), the
voter-selection mode (`loo_voter_search_mode`), and a review flag
(`loo_review_candidate`) for a placement whose conflict is strong but stays inside the
normal spread between markers.

The review flag marks a copy for inspection; it does not confirm contamination.
Contamination from a close relative stays unresolved when its placement falls inside that
normal spread.

The fixed tree and sequence instruments test mechanism and safety properties. They do not
estimate biological performance. See
[`DEVELOPMENT.md`](DEVELOPMENT.md#contamination-detection-evidence-instruments) for the
commands and evidence boundary.

## Troubleshooting

- **`Input path does not exist`:** check the path passed to `--genomedir` or `--modeldir`.
  Use `.faa` for proteins or `.fna`, `.fa`, or `.fasta` for assemblies. The marker input
  must contain valid HMM entries.
- **Custom HMMs return no hits:** curated marker sets can use the default gathering cutoff.
  For models without gathering cutoffs, start with `--hmmsearch_cutoff evalue`.
- **Retrying an interrupted run:** repeat the command with the same `--outdir`. SGTree
  records ownership in `.sgtree-run` and clears its generated outputs before retrying.
  It refuses reserved output names in a directory it cannot identify as an SGTree run;
  use a new output directory in that case.
- **Reference inputs changed:** SGTree compares the contents of reference FASTA and marker
  HMM files with the cache metadata and rebuilds a stale or incomplete cache. A cache
  from an older version is rebuilt if its metadata and required files identify it as an
  intact SGTree cache. Otherwise, use a new cache directory with `--ref_concat`.
- **Expected intermediates are missing:** SGTree archives them by default. Add
  `--keep_intermediates yes` while debugging.
- **A small marker set produces only LOO abstentions:** `loo_profile` needs at least five
  voter trees and six shared coordinate taxa. An abstention is a recorded outcome, not a
  confirmed clean marker.

## Documentation

- [`DEVELOPMENT.md`](DEVELOPMENT.md): repository layout, tests, benchmarks, and the
  contamination-detection evidence instruments
- [`CHANGELOG.md`](CHANGELOG.md): release history

## Authors

| Author | Contact | Since |
|---|---|---|
| Ewan Whittaker-Walker | ewanww@berkeley.edu | 2019 |
| Frederik Schulz | fschulz@lbl.gov | 2019 |
| Juan C. Villada | jvillada@lbl.gov | 2021 |
| Marianne Buscaglia | mbuscaglia@lbl.gov | 2022 |

SGTree is developed at the DOE Joint Genome Institute, Lawrence Berkeley National
Laboratory.

## License

SGTree is released under a Berkeley Lab non-commercial use licence. See
[`LICENSE`](LICENSE) for the terms.

SGTree Copyright (c) 2026, The Regents of the University of California, through
Lawrence Berkeley National Laboratory. All rights reserved. For questions about
rights to use or distribute this software, contact Berkeley Lab's Intellectual
Property Office at IPO@lbl.gov.
