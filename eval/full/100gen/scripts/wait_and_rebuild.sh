#!/usr/bin/env bash
set -euo pipefail

cd /home/fschulz/dev/software/sgtree

while true; do
  count="$(find runs/contig_variant_proof_of_concept/benchmarks -type f -name benchmark_manifest.json 2>/dev/null | wc -l)"
  if [ "${count}" -ge 24 ]; then
    break
  fi
  sleep 30
done

PYTHONPATH=. python runs/contig_variant_proof_of_concept/rebuild_eval_100gen_exact.py
