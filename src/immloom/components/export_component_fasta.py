#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build component FASTA files from global block components."
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. dataset-02")
    parser.add_argument("--locus", required=True, help="Locus name, e.g. IGH")
    parser.add_argument(
        "--data-root",
        default="data/processed",
        help="Root directory containing processed datasets",
    )
    return parser.parse_args()


def rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def read_fasta(path: Path) -> str:
    with open(path) as f:
        lines = f.read().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def load_fasta_sequences(fasta_dir: Path, locus: str) -> dict[str, str]:
    fasta_dict: dict[str, str] = {}

    for path in sorted(fasta_dir.glob(f"*-{locus}.fasta")):
        sample_name = path.name.removesuffix(f"-{locus}.fasta")
        fasta_dict[sample_name] = read_fasta(path)

    return fasta_dict


def load_block_tables(blocks_dir: Path) -> dict[str, pd.DataFrame]:
    blocks_cache: dict[str, pd.DataFrame] = {}

    for path in sorted(blocks_dir.glob("blocks_*.tsv")):
        sample_name = path.stem.removeprefix("blocks_")
        blocks_cache[sample_name] = pd.read_csv(path, sep="\t")

    return blocks_cache


def main() -> None:
    args = parse_args()

    dataset = args.dataset
    locus = args.locus
    data_root = Path(args.data_root)

    fasta_dir = data_root / dataset / "fasta" / "mixed_forward"
    blocks_dir = data_root / dataset / locus / "blocks" / "tables"
    components_path = data_root / dataset / locus / "components" / "tables" / "components.tsv"
    out_dir = data_root / dataset / locus / "components" / "fasta"

    if not fasta_dir.exists():
        raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")
    if not blocks_dir.exists():
        raise FileNotFoundError(f"Blocks directory not found: {blocks_dir}")
    if not components_path.exists():
        raise FileNotFoundError(f"Components table not found: {components_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    components = pd.read_csv(components_path, sep="\t")
    if components.empty:
        print("components.tsv is empty. Nothing to export.")
        return

    fasta_dict = load_fasta_sequences(fasta_dir, locus)
    blocks_cache = load_block_tables(blocks_dir)

    for comp_id, group_df in components.groupby("component_id"):
        sequences: list[tuple[str, str]] = []

        for _, row in group_df.iterrows():
            group = row["group"]
            block_id = int(row["block_id"])

            if group not in fasta_dict:
                raise FileNotFoundError(
                    f"Missing FASTA sequence for sample '{group}' in {fasta_dir}"
                )
            if group not in blocks_cache:
                raise FileNotFoundError(
                    f"Missing block table for sample '{group}' in {blocks_dir}"
                )

            df_blocks = blocks_cache[group]
            block_rows = df_blocks[df_blocks["block_id"] == block_id]

            if block_rows.empty:
                print(
                    f"Warning: no rows found for component={comp_id}, "
                    f"sample={group}, block_id={block_id}"
                )
                continue

            full_seq = fasta_dict[group]

            for _, b in block_rows.iterrows():
                d1 = int(b["d1"])
                d2 = int(b["d2"])
                forward = bool(b["forward"])

                # coordinates are assumed 1-based inclusive
                seq = full_seq[d1 - 1 : d2]

                if not forward:
                    seq = rev_comp(seq)

                header = f">{group}_block{block_id}_{d1}_{d2}_{'F' if forward else 'RC'}"
                sequences.append((header, seq))

        out_path = out_dir / f"component_{int(comp_id)}.fasta"
        with open(out_path, "w") as f:
            for header, seq in sequences:
                f.write(header + "\n")
                f.write(seq + "\n")

        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()