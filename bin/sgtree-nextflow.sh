#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'HELP'
SGTree (Nextflow runner)

Usage:
  pixi run sgtree --genomedir <dir> --modeldir <marker_set.hmm> [options]

Required arguments:
  --genomedir <dir>      Directory with input proteomes (*.faa) or assemblies (*.fna)
  --modeldir <file>      Combined marker-set HMM file (for example, resources/models/UNI56.hmm)

Common options:
  --outdir <dir>         Output directory (default: runs/nextflow/default)
  --marker_selection     true|false (default: false)
  --ref <dir>            Reference proteomes directory
  --ani_cluster          yes|no|true|false (default: false)
  --snp                  yes|no|true|false (default: false; requires --ani_cluster yes)
                         SNP trees use only shared UNI56-bearing backbone contigs with >=95% ANI
  --ani_threshold <f>    ANI cutoff used before MCL clustering (default: 95.0)
  --ani_backend <name>   auto | skani | minimap2 (default: auto)
  --ani_mcl_inflation <f>
                         MCL inflation used after ANI filtering (default: 2.0)
  --snp_tree_min_cluster_size <n>
                         Build SNP trees only for ANI clusters of at least n genomes (default: 3)
  --singles              yes|no|true|false (default: false)
  --percent_models <n>   Minimum marker coverage threshold (default: 10)
  --max_sdup <n>         Max copies allowed for one marker in a genome (-1 disables)
  --max_dupl <f>         Max fraction of duplicated markers in a genome (-1 disables)
  --num_nei <n>          Optional singleton-removal neighbor count (0=auto)
  --singles_mode         neighbor | delta_rf | backbone | ensemble (default: neighbor)
  --singles_min_rfdist   Minimum marker/global RF needed before singleton pruning (default: 0.25)
  --hmmsearch_cutoff     cut_ga | cut_tc | cut_nc | evalue (default: cut_ga)
  --hmmsearch_evalue <f> E-value used when --hmmsearch_cutoff evalue (default: 1e-5)
  --aln <method>         hmmalign | mafft | mafft-linsi (default: hmmalign)
  --tree_method <name>   fasttree | iqtree (default: fasttree)
  --iqtree_fast          yes|no|true|false (default: true)
  --iqtree_model <str>   IQ-TREE model string (default: LG+F+I+G4)
  --selection_mode       coordinate | legacy (default: coordinate)
  --selection_max_rounds Max coordinate-descent rounds (default: 5)
  --selection_global_rounds Fixed guide-tree rebuild rounds for marker selection (default: 1)
  --lock_references      yes|no|true|false (default: false)
  --hmmsearch_cpus <n>   CPUs for HMM search (default: 8)
  --align_cpus <n>       CPUs for alignment steps (default: 4)
  --fasttree_cpus <n>    CPUs for tree-building steps (default: 1)

Examples:
  pixi run sgtree --genomedir testgenomes/Chloroflexi --modeldir resources/models/UNI56.hmm

  pixi run sgtree \
    --genomedir testgenomes/Chloroflexi \
    --modeldir resources/models/UNI56.hmm \
    --outdir runs/nextflow/manual_full \
    --marker_selection true \
    --ref testgenomes/chlorref \
    --singles yes

Notes:
  - Logs are written automatically to runs/nextflow/logs/.
  - Nucleotide (`*.fna`) inputs, `--ani_cluster yes`, and `--snp yes` are routed to the Python engine automatically.
  - For the Python implementation directly, run: pixi run sgtree-python ...
HELP
}

if [[ $# -eq 0 ]]; then
  show_help
  exit 0
fi

case "${1:-}" in
  -h|--help|help)
    show_help
    exit 0
    ;;
esac

mkdir -p runs/nextflow/logs
log_file="runs/nextflow/logs/sgtree_$(date +%Y%m%d_%H%M%S).log"

echo "[sgtree] nextflow log: ${log_file}"

args=("$@")
genomedir=""
modeldir=""
refdir=""
ani_cluster="no"
snp="no"

for ((i=0; i<${#args[@]}; i++)); do
  case "${args[i]}" in
    --genomedir)
      genomedir="${args[i+1]:-}"
      i=$((i + 1))
      ;;
    --modeldir)
      modeldir="${args[i+1]:-}"
      i=$((i + 1))
      ;;
    --ref)
      refdir="${args[i+1]:-}"
      i=$((i + 1))
      ;;
    --ani_cluster)
      ani_cluster="${args[i+1]:-no}"
      i=$((i + 1))
      ;;
    --snp)
      snp="${args[i+1]:-no}"
      i=$((i + 1))
      ;;
  esac
done

is_true() {
  case "${1,,}" in
    yes|true|1) return 0 ;;
    *) return 1 ;;
  esac
}

looks_like_fna_input() {
  local path="$1"
  if [[ -z "$path" ]]; then
    return 1
  fi
  if [[ -f "$path" ]]; then
    case "${path##*.}" in
      fna|fa|fasta) return 0 ;;
    esac
    return 1
  fi
  if [[ -d "$path" ]]; then
    shopt -s nullglob
    local matches=("$path"/*.fna "$path"/*.fa "$path"/*.fasta)
    shopt -u nullglob
    if (( ${#matches[@]} > 0 )); then
      return 0
    fi
  fi
  return 1
}

if is_true "$ani_cluster" || is_true "$snp" || looks_like_fna_input "$genomedir" || looks_like_fna_input "$refdir"; then
  rest_args=()
  skip_next=0
  for ((i=0; i<${#args[@]}; i++)); do
    if (( skip_next )); then
      skip_next=0
      continue
    fi
    case "${args[i]}" in
      --genomedir|--modeldir)
        skip_next=1
        ;;
      *)
        rest_args+=("${args[i]}")
        ;;
    esac
  done
  exec python -m sgtree "$genomedir" "$modeldir" "${rest_args[@]}"
fi

exec ./nextflow -log "${log_file}" run main.nf "$@"
