#!/usr/bin/env python3
"""
Build pairwise block tables for all pairwise alignment files.

Usage:
    python build_pair_tables.py --dataset dataset_01_subset --locus IGH
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from immloom import (
    preprocess_alignment,
    filter_segments,
    split_segments,
    build_bipartite_block_graph,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build pairwise block tables from pairwise alignment files."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset name, e.g. dataset_01_subset",
    )
    parser.add_argument(
        "--locus",
        required=True,
        help="Locus name, e.g. IGH",
    )
    return parser.parse_args()


def get_alignment_files(dataset: str, locus: str) -> list[str]:
    pairwise_dir = Path("data/processed") / dataset / "patchwork_output_reversed" / locus / "pairwise_alignments"

    if not pairwise_dir.exists():
        raise FileNotFoundError(f"Pairwise alignment directory not found: {pairwise_dir}")

    files = [
        p.name
        for p in pairwise_dir.iterdir()
        if p.is_file() and not p.name.startswith("self")
    ]
    return sorted(files)


def extract_sample_names(file_name: str) -> tuple[str, str]:
    """
    Expected original logic:
        name_x = '_'.join(file_name.split('-')[1].split('_')[:2])
        name_y = file_name.split('-')[2].split('.')[0]
    """
    parts = file_name.split("-")
    if len(parts) < 3:
        raise ValueError(f"Unexpected filename format: {file_name}")

    name_x = "_".join(parts[1].split("_")[:2])
    name_y = parts[2].split(".")[0]
    return name_x, name_y


def process_one_file(file_name: str, dataset: str, locus: str) -> None:
    name_x, name_y = extract_sample_names(file_name)

    alignment_path = (
        Path("data/processed")
        / dataset
        / "patchwork_output_reversed"
        / locus
        / "pairwise_alignments"
        / file_name
    )

    dfx_path = Path("results") / dataset / "self" / locus / "tables" / f"blocks_{name_x}.tsv"
    dfy_path = Path("results") / dataset / "self" / locus / "tables" / f"blocks_{name_y}.tsv"

    out_dir = Path("results") / dataset / "pair" / locus / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{name_x}-{name_y}.tsv"

    if not dfx_path.exists():
        raise FileNotFoundError(f"Missing self-block table: {dfx_path}")
    if not dfy_path.exists():
        raise FileNotFoundError(f"Missing self-block table: {dfy_path}")

    df = preprocess_alignment(alignment_path)
    df = filter_segments(df, length=5000, pi=80.0)

    dfx = pd.read_csv(dfx_path, sep="\t")
    dfy = pd.read_csv(dfy_path, sep="\t")

    sx_list = list(dfx["d1"]) + list(dfx["d2"])
    sy_list = list(dfy["d1"]) + list(dfy["d2"])

    result_df = split_segments(df, sx_list=sx_list, sy_list=sy_list, gap=0.8)
    result_df = result_df[result_df["length"] >= 2500].reset_index(drop=True)

    if result_df.empty:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(
            out_path, index=False, sep="\t"
        )
        return

    G = build_bipartite_block_graph(dfx, dfy, result_df)

    edges = []
    for u, v, data in G.edges(data=True):
        edges.append(
            {
                "x_block_id": u[1],
                "y_block_id": v[1],
                "weight": data.get("weight"),
                "parent_ids": data.get("parent_ids"),
            }
        )

    edges_df = pd.DataFrame(edges)

    if not edges_df.empty and "weight" in edges_df.columns:
        edges_df = edges_df.sort_values("weight", ascending=False).reset_index(drop=True)

    if {"x_block_id", "y_block_id"}.issubset(edges_df.columns):
        edges_df[["x_block_id", "y_block_id"]].to_csv(out_path, index=False, sep="\t")
    else:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(
            out_path, index=False, sep="\t"
        )


def main() -> None:
    args = parse_args()

    files = get_alignment_files(args.dataset, args.locus)

    if not files:
        print("No pairwise alignment files found.")
        return

    for file_name in tqdm(files, desc="Processing pairwise alignments"):
        try:
            process_one_file(file_name, args.dataset, args.locus)
        except Exception as e:
            print(f"[ERROR] Failed for {file_name}: {e}")


if __name__ == "__main__":
    main()