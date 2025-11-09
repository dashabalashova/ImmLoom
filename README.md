# ImmLoom

copy data
```
git clone https://github.com/yana-safonova/primate_t2t_igtr.git tmp_repo
mkdir -p data/raw/dataset_01/fasta
cp tmp_repo/data_primate_igtrloci_fasta/*.fasta data/raw/dataset_01/fasta/
rm -rf tmp_repo
```

venv
```
brew install python@3.10
python3.10 -m venv ~/.venvs/immloom
source ~/.venvs/immloom/bin/activate
python -m pip install --upgrade pip setuptools wheel
pip install uv
uv pip install numpy pandas matplotlib ipykernel patchworkplot
python -m ipykernel install --user --name immloom --display-name "immloom"
brew install lastz
```

run patchworkplot
```
python3 scripts/normalize_fasta_names.py
python3 scripts/create_ig_configs.py
mkdir -p data/processed/dataset_01/patchwork_output/

patchworkplot -i data/processed/dataset_01/configs/config_IGH.csv -o data/processed/dataset_01/patchwork_output/IGH --cmap PuBuGn --reverse-cmap false --lwidth 2 --lower
```

run immloom
```
uv pip install -e .

F_PATH=data/processed/dataset_01/patchwork_output/IGH/pairwise_alignments/
F_NAME=self_0-s01_h1_IGH.tsv
python3 scripts/run_pipeline.py --input-path ${F_PATH} --input-file ${F_NAME} --out-dir ./immloom_out/self_0-s01_h1_IGH

chmod +x scripts/run_all.sh
scripts/run_all.sh
```