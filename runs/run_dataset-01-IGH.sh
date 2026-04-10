#!/usr/bin/env bash
set -euo pipefail

DATASET="dataset-01"
LOCUS="IGH"
VENV="${HOME}/.venvs/immloom"

echo "==> Install Python 3.10 via Homebrew"
brew install python@3.10

echo "==> Create virtual environment"
python3.10 -m venv "${VENV}"

echo "==> Activate virtual environment"
# shellcheck disable=SC1090
source "${VENV}/bin/activate"

echo "==> Upgrade packaging tools"
python -m pip install --upgrade pip setuptools wheel
pip install uv

echo "==> Install Python dependencies"
uv pip install \
  gdown \
  numpy \
  pandas \
  matplotlib \
  ipykernel \
  patchworkplot \
  tqdm \
  biopython \
  networkx

echo "==> Register Jupyter kernel"
python -m ipykernel install --user --name immloom --display-name "immloom"

echo "==> Install external tools"
brew install lastz
brew install clustal-omega

echo "==> Install immloom package in editable mode"
uv pip install -e .

echo "==> Create raw subset dataset"
python preprocessing/create_dataset-01.py
python preprocessing/process_dataset-01.py

echo "==> Make patchworkplot config"
python -m immloom.pwp.make_pwp_config --dataset "${DATASET}"

echo "==> Run patchworkplot"
python -m immloom.pwp.run_patchworkplot --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Reorient sequences to forward orientation"
python -m immloom.pwp.reorient_to_forward --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Build blocks"
python -m immloom.blocks.build_blocks --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Connect blocks"
python -m immloom.blocks.connect_blocks --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Build component graph"
python -m immloom.blocks.build_graph --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Export component FASTA"
python -m immloom.components.export_component_fasta --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Run Clustal Omega"
python -m immloom.components.run_clustalo --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Build MSA graphs"
python -m immloom.components.build_msa_figures --dataset "${DATASET}" --locus "${LOCUS}"

echo "==> Done"