#!/usr/bin/env bash
set -euo pipefail

# Укажите путь к папке с tsv (от корня проекта)
F_PATH="data/processed/dataset_01/patchwork_output/IGH/pairwise_alignments"

# Опция: если пакет не установлен, укажите src в PYTHONPATH
export PYTHONPATH="${PWD}/src:${PYTHONPATH:-}"

# куда складывать результаты по каждому файлу (директория создаётся автоматически)
OUT_ROOT="./immloom_out"

mkdir -p "$OUT_ROOT"

# Проход по всем .tsv (не захватываем скрытые файлы)
shopt -s nullglob
for tsv in "$F_PATH"/self*.tsv; do
  name=$(basename "$tsv" .tsv)
  echo "----------------------------------------"
  echo "Running pipeline for: $name"
  outdir="${OUT_ROOT}/${name}"
  mkdir -p "$outdir"

  # вызываем ваш скрипт
  python3 scripts/run_pipeline.py --input-path "$F_PATH" --input-file "${name}.tsv" --out-dir "$outdir"

  echo "Finished $name -> output saved to $outdir"
done
echo "All done."
