#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat <<'HELP'
SGTree (Nextflow runner)

Usage:
  pixi run sgtree --genomedir <dir> --modeldir <dir> [options]

Required arguments:
  --genomedir <dir>      Directory with input proteomes (*.faa)
  --modeldir <dir>       Directory with marker models (*.hmm)

Common options:
  --outdir <dir>         Output directory (default: runs/nextflow/default)
  --marker_selection     true|false (default: false)
  --ref <dir>            Reference proteomes directory
  --singles              yes|no|true|false (default: false)
  --percent_models <n>   Minimum marker coverage threshold (default: 10)
  --aln <method>         hmmalign | mafft | mafft-linsi (default: hmmalign)
  --hmmsearch_cpus <n>   CPUs for HMM search (default: 8)
  --align_cpus <n>       CPUs for alignment steps (default: 4)
  --fasttree_cpus <n>    CPUs for FastTree steps (default: 1)

Examples:
  pixi run sgtree --genomedir testgenomes/Chloroflexi --modeldir hmms/UNI56

  pixi run sgtree \
    --genomedir testgenomes/Chloroflexi \
    --modeldir hmms/UNI56 \
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
