#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_self_all.sh [--dataset DATASET] [--locus LOCUS] [--in-root PATH] [--out-root PATH]

Defaults:
  --dataset  dataset_01
  --locus    IGH
  --in-root  data/processed
  --out-root outputs/self

Examples:
  ./run_self_all.sh
  ./run_self_all.sh --dataset dataset_02 --locus IGK
  ./run_self_all.sh --dataset dataset_01 --locus IGL --out-root ./outputs/custom
EOF
}

# -----------------------
# Defaults
# -----------------------
DATASET="dataset_01"
LOCUS="IGH"
BASE_IN="data/processed"
BASE_OUT="outputs"

# -----------------------
# Args parsing
# -----------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset|-d)
      DATASET="${2:?Missing value for --dataset}"; shift 2 ;;
    --locus|-l)
      LOCUS="${2:?Missing value for --locus}"; shift 2 ;;
    --in-root)
      BASE_IN="${2:?Missing value for --in-root}"; shift 2 ;;
    --out-root)
      BASE_OUT="${2:?Missing value for --out-root}"; shift 2 ;;
    --help|-h)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2 ;;
  esac
done

# -----------------------
# Paths
# -----------------------
PATH_IN="${BASE_IN}/${DATASET}/patchwork_output/${LOCUS}/pairwise_alignments"
PATH_OUT="${BASE_OUT}/${DATASET}/${LOCUS}/self"

# Option: if package isn't installed, add src to PYTHONPATH
export PYTHONPATH="${PWD}/src:${PYTHONPATH:-}"

mkdir -p "$PATH_OUT"

# -----------------------
# Main loop
# -----------------------
shopt -s nullglob
found_any=false
for tsv in "${PATH_IN}"/self*.tsv; do
  found_any=true
  name=$(basename "$tsv" .tsv)

  echo "----------------------------------------"
  echo "Running pipeline for: ${DATASET} | ${LOCUS} | ${name}"

  outdir="${PATH_OUT}/${name}"
  mkdir -p "$outdir"

  python3 scripts/run_self.py \
    --input-path "$PATH_IN" \
    --input-file "${name}.tsv" \
    --out-dir "$outdir"

  echo "Finished ${name} -> output saved to ${outdir}"
done

if [[ "$found_any" == false ]]; then
  echo "No files matched: ${PATH_IN}/pairwise_alignments/self*.tsv" >&2
  exit 1
fi

echo "All done for ${DATASET} (${LOCUS}). Outputs: ${PATH_OUT}"
