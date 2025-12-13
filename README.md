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

## 3. run patchworkplot

code: patchworkplot
input: data/processed/${DATASET}/configs/config_${LOCUS}.csv
output: data/processed/${DATASET}/patchwork_output/${LOCUS}

```
python3 scripts/normalize_fasta_names.py
python3 scripts/create_ig_configs.py
mkdir -p data/processed/dataset_01/patchwork_output/IGH

patchworkplot -i data/processed/dataset_01/configs/config_IGH.csv -o data/processed/dataset_01/patchwork_output/IGH --cmap PuBuGn --reverse-cmap false --lwidth 2 --lower
```

## 4. run immloom self

code: scripts/run_self_all.sh
input: data/processed/${DATASET}/patchwork_output/${LOCUS}/pairwise_alignments
output: outputs/${DATASET}/${LOCUS}/self

```
uv pip install -e .

# test
F_PATH=data/processed/dataset_01/patchwork_output/IGH/pairwise_alignments
F_NAME=self_0-s01_h1_IGH.tsv
python3 scripts/run_self.py --input-path ${F_PATH} --input-file ${F_NAME} --out-dir outputs/dataset_01/IGH/self/self_0-s01_h1_IGH

# run
chmod +x scripts/run_self_all.sh
scripts/run_self_all.sh --dataset dataset_01 --locus IGH
```

## 5. run immloom pairs

```
chmod +x scripts/run_pairs_all.sh
scripts/run_pairs_all.sh --dataset dataset_01 --locus IGH
```

## 6. connected components

code: notebooks/connected_components.ipynb
input: outputs/${DATASET}/${LOCUS}/graphs
output-1: outputs/${DATASET}/${LOCUS}/graphs/G-merged.tsv
output-2: outputs/${DATASET}/${LOCUS}/components_sns
output-3: outputs/${DATASET}/${LOCUS}/components/sequences

## 7. clustalo

```
source ~/.venvs/immloom/bin/activate
python scripts/run_clustalo.py
```

## 8. alignment visualization

code: notebooks/alignments.ipynb
input: outputs/dataset_01/IGH/components/alignments
