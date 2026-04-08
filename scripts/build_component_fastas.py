#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build component FASTA files from blocks (1-based inclusive coordinates)"
    )
    parser.add_argument("--dataset", required=True, help="Dataset name")
    parser.add_argument("--locus", required=True, help="Locus name (e.g. IGH)")
    return parser.parse_args()


# === reverse complement ===
def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# === read fasta (single sequence) ===
def read_fasta(path: Path) -> str:
    with open(path) as f:
        lines = f.read().splitlines()
    return "".join([l for l in lines if not l.startswith(">")])


def main():
    args = parse_args()
    dataset = args.dataset
    locus = args.locus

    base_results = Path(f"results/{dataset}")
    base_fasta = Path(f"data/processed/{dataset}/fasta_reversed")
    blocks_dir = base_results / f"self/{locus}/tables"
    components_path = base_results / f"multi/{locus}/tables/components.tsv"

    out_dir = base_results / f"components/{locus}/fasta"
    out_dir.mkdir(parents=True, exist_ok=True)

    # === load components ===
    components = pd.read_csv(components_path, sep="\t")

    # === load all FASTA sequences ===
    fasta_dict = {}
    for fname in os.listdir(base_fasta):
        if fname.endswith(f"_{locus}.fasta"):
            key = fname.replace(f"_{locus}.fasta", "")
            fasta_dict[key] = read_fasta(base_fasta / fname)

    # === cache block tables ===
    blocks_cache = {}
    for fname in os.listdir(blocks_dir):
        if fname.startswith("blocks_") and fname.endswith(".tsv"):
            group = fname.replace("blocks_", "").replace(".tsv", "")
            blocks_cache[group] = pd.read_csv(blocks_dir / fname, sep="\t")

    # === main loop ===
    for comp_id, group_df in components.groupby("component_id"):
        sequences = []

        for _, row in group_df.iterrows():
            group = row["group"]
            block_id = row["block_id"]

            df_blocks = blocks_cache[group]
            block_rows = df_blocks[df_blocks["block_id"] == block_id]

            for _, b in block_rows.iterrows():
                d1, d2 = int(b["d1"]), int(b["d2"])
                forward = b["forward"]

                # === IMPORTANT: 1-based inclusive ===
                seq = fasta_dict[group][d1 - 1 : d2]

                if not forward:
                    seq = rev_comp(seq)

                header = f">{group}_block{block_id}_{d1}_{d2}_{'F' if forward else 'RC'}"
                sequences.append((header, seq))

        # === write fasta ===
        out_path = out_dir / f"component_{comp_id}.fasta"
        with open(out_path, "w") as f:
            for h, s in sequences:
                f.write(h + "\n")
                f.write(s + "\n")

        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()