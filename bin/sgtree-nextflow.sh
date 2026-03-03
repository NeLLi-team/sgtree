#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'HELP'
SGTree (Nextflow runner)

Usage:
  pixi run sgtree --genomedir <dir> --modeldir <marker_set.hmm> [options]

Required arguments:
  --genomedir <dir>      Directory with input proteomes (*.faa)
  --modeldir <file>      Combined marker-set HMM file (for example, resources/models/UNI56.hmm)

Common options:
  --outdir <dir>         Output directory (default: runs/nextflow/default)
  --marker_selection     true|false (default: false)
  --ref <dir>            Reference proteomes directory
  --singles              yes|no|true|false (default: false)
  --percent_models <n>   Minimum marker coverage threshold (default: 10)
  --max_sdup <n>         Max copies allowed for one marker in a genome (-1 disables)
  --max_dupl <f>         Max fraction of duplicated markers in a genome (-1 disables)
  --num_nei <n>          Optional singleton-removal neighbor count (0=auto)
  --hmmsearch_cutoff     cut_ga | cut_tc | cut_nc | evalue (default: cut_ga)
  --hmmsearch_evalue <f> E-value used when --hmmsearch_cutoff evalue (default: 1e-5)
  --aln <method>         hmmalign | mafft | mafft-linsi (default: hmmalign)
  --tree_method <name>   fasttree | iqtree (default: fasttree)
  --iqtree_fast          yes|no|true|false (default: true)
  --iqtree_model <str>   IQ-TREE model string (default: LG+F+I+G4)
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
  - For the Python implementation, run: pixi run sgtree-python ...
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

exec ./nextflow -log "${log_file}" run main.nf "$@"
