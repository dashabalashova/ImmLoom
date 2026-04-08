# ImmLoom

## 1. set up the Python environment

```
brew install python@3.10
python3.10 -m venv ~/.venvs/immloom
source ~/.venvs/immloom/bin/activate
python -m pip install --upgrade pip setuptools wheel
pip install uv
```

## 2. download and unpack the raw dataset

```
uv pip install gdown
mkdir -p data/raw/dataset_01/
gdown -O data/raw/dataset_01/primate_igdetective_outputs.zip "https://drive.google.com/uc?id=1hPsaegHbII8sGf6toy8Hr0NBM80OtDZw"
unzip data/raw/dataset_01/primate_igdetective_outputs.zip -d data/raw/dataset_01
```

## 3. make smaller subset of samples

```
dataset="dataset_01_subset"
src="data/raw/dataset_01/igdetective_outputs"
dst="data/raw/${dataset}/igdetective_outputs"
mkdir -p "$dst"
samples=(
  AG05253_1_hap1 AG05253_1_hap2
  AG06939_13_hap1 AG06939_13_hap2
  PR00251_hap1 PR00251_hap2
  PR00366_hap1 PR00366_hap2
)
for sample in "${samples[@]}"; do
  src_dir="${src}/${sample}_igdetective"
  dst_dir="${dst}/"
  echo "Copying $sample..."
  cp -r "$src_dir" "$dst_dir"
done
```

## 4. install project dependencies

```
uv pip install numpy pandas matplotlib ipykernel patchworkplot tqdm biopython
python -m ipykernel install --user --name immloom --display-name "immloom"
brew install lastz
uv pip install -e .
```

## 5. build the locus–contig table

```
dataset="dataset_01_subset"
python3 scripts/build_locus_contig_table.py --dataset "$dataset"
```

## 6. run PatchWorkPlot for the locus

```
locus="IGH"
outdir="data/processed/${dataset}/patchwork_output/${locus}"
config="data/processed/${dataset}/config/config_${locus}.csv"
mkdir -p "$outdir"
patchworkplot \
  -i "$config" \
  -o "$outdir" \
  --cmap PuBuGn \
  --reverse-cmap false \
  --lwidth 2 \
  --lower
```

## 7. reverse the patchwork orientation

```
python3 scripts/reverse_patchwork.py --dataset "$dataset" --locus "$locus"
```

## 8. build blocks, pair tables, and component FASTAs

```
uv pip install networkx
python scripts/self_blocks.py --dataset "$dataset" --locus "$locus"
python scripts/build_pair_tables.py --dataset "$dataset" --locus "$locus"
python scripts/multi_blocks.py --dataset "$dataset" --locus "$locus"
python scripts/build_component_fastas.py --dataset "$dataset" --locus "$locus"
```

## 9. run multiple sequence alignment with Clustal Omega

```
brew install clustal-omega
mkdir -p results/${dataset}/clustalo/${locus}
for fasta in results/${dataset}/components/${locus}/fasta/component_*.fasta
do
  base=$(basename "$fasta" .fasta)
  clustalo \
    -i "$fasta" \
    -o results/${dataset}/clustalo/${locus}/${base}.aln \
    --outfmt=clu \
    --force
  echo "Done $base"
done
```

## 10. build alignment graphs

```
python scripts/build_alignment_graphs.py --dataset "$dataset" --locus "$locus"
```