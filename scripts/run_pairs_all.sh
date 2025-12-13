#!/usr/bin/env bash
set -euo pipefail

# =======================
# Defaults
# =======================
DATASET="dataset_01"
LOCUS="IGH"
MAX=7

# =======================
# Args parsing
# =======================
usage() {
  cat <<'EOF'
Usage:
  ./run_pairs_all.sh [--dataset DATASET] [--locus LOCUS] [--max MAX]

Defaults:
  --dataset  dataset_01
  --locus    IGH
  --max      7

Example:
  ./run_pairs_all.sh --dataset dataset_02 --locus IGK
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset|-d)
      DATASET="${2:?Missing value for --dataset}"
      shift 2
      ;;
    --locus|-l)
      LOCUS="${2:?Missing value for --locus}"
      shift 2
      ;;
    --max)
      MAX="${2:?Missing value for --max}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

# =======================
# Main loop
# =======================
for (( n1=0; n1<MAX; n1++ )); do
  for (( n2=n1+1; n2<=MAX; n2++ )); do
    echo ">>> Running: n1=${n1}, n2=${n2}, dataset=${DATASET}, locus=${LOCUS}"

    python scripts/run_pairs.py \
      --n1 "$n1" \
      --n2 "$n2" \
      --dataset "$DATASET" \
      --locus "$LOCUS"
  done
done
