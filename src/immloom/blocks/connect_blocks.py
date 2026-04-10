#!/usr/bin/env python3
"""Build pairwise block tables from pairwise alignments.

Example:
    python -m immloom.blocks.pair_blocks --dataset dataset-02 --locus IGH
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from immloom.geometry.segments import split_segments
from immloom.graph.blocks import build_bipartite_block_graph
from immloom.io.alignment import filter_segments, preprocess_alignment


PAIR_MIN_LEN = 5_000
MIN_PI = 80.0
SPLIT_GAP = 0.8
FINAL_MIN_LEN = 2_500


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build pairwise block tables from pairwise alignment files."
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. dataset-02")
    parser.add_argument("--locus", required=True, help="Locus name, e.g. IGH")
    parser.add_argument(
        "--data-root",
        default="data/processed",
        help="Root directory containing processed datasets",
    )
    parser.add_argument(
        "--output-root",
        default="results",
        help="Root directory for output results",
    )
    return parser.parse_args()


def get_alignment_files(dataset: str, locus: str, data_root: str) -> list[Path]:
    pairwise_dir = (
        Path(data_root)
        / dataset
        / locus
        / "pwp"
        / "pairwise_alignments_forward"
    )

    if not pairwise_dir.exists():
        raise FileNotFoundError(f"Pairwise alignment directory not found: {pairwise_dir}")

    files = [
        p
        for p in pairwise_dir.iterdir()
        if p.is_file() and p.name.startswith("pair")
    ]
    return sorted(files)


def extract_sample_names(file_name: str) -> tuple[str, str]:
    """
    Example:
        pair_0-s01-h1_1-s01-h2.tsv -> ("s01-h1", "s01-h2")
    """
    stem = Path(file_name).stem
    parts = stem.split("_")
    if len(parts) < 3:
        raise ValueError(f"Unexpected filename format: {file_name}")

    x_part = parts[1]  # 0-s01-h1
    y_part = parts[2]  # 1-s01-h2

    name_x = x_part.split("-", 1)[1]
    name_y = y_part.split("-", 1)[1]
    return name_x, name_y


def process_one_file(file_path: Path, dataset: str, locus: str, data_root: str, output_root: str) -> None:
    name_x, name_y = extract_sample_names(file_path.name)

    dfx_path = (
        Path(output_root)
        / dataset
        / locus
        / "blocks"
        / "tables"
        / f"blocks_{name_x}.tsv"
    )
    dfy_path = (
        Path(output_root)
        / dataset
        / locus
        / "blocks"
        / "tables"
        / f"blocks_{name_y}.tsv"
    )

    out_dir = Path(output_root) / dataset / locus / "block_pairs" / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{file_path.stem}.tsv"

    if not dfx_path.exists():
        raise FileNotFoundError(f"Missing self-block table for x sample {name_x}: {dfx_path}")
    if not dfy_path.exists():
        raise FileNotFoundError(f"Missing self-block table for y sample {name_y}: {dfy_path}")

    df = preprocess_alignment(file_path)
    df = filter_segments(df, length=PAIR_MIN_LEN, pi=MIN_PI)

    if df.empty:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(out_path, index=False, sep="\t")
        return

    dfx = pd.read_csv(dfx_path, sep="\t")
    dfy = pd.read_csv(dfy_path, sep="\t")

    if dfx.empty or dfy.empty:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(out_path, index=False, sep="\t")
        return

    sx_list = list(dfx["d1"]) + list(dfx["d2"])
    sy_list = list(dfy["d1"]) + list(dfy["d2"])

    result_df = split_segments(df, sx_list=sx_list, sy_list=sy_list, gap=SPLIT_GAP)
    if result_df.empty:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(out_path, index=False, sep="\t")
        return

    result_df = result_df[result_df["length"] >= FINAL_MIN_LEN].reset_index(drop=True)
    if result_df.empty:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(out_path, index=False, sep="\t")
        return

    graph = build_bipartite_block_graph(dfx, dfy, result_df)

    edges = []
    for u, v, data in graph.edges(data=True):
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
        edges_df = edges_df.sort_values(
            ["weight", "x_block_id", "y_block_id"],
            ascending=[False, True, True],
        ).reset_index(drop=True)

    if {"x_block_id", "y_block_id"}.issubset(edges_df.columns):
        edges_df[["x_block_id", "y_block_id"]].to_csv(out_path, index=False, sep="\t")
    else:
        pd.DataFrame(columns=["x_block_id", "y_block_id"]).to_csv(out_path, index=False, sep="\t")


def main() -> None:
    args = parse_args()

    files = get_alignment_files(args.dataset, args.locus, args.data_root)

    if not files:
        print("No pairwise alignment files found.")
        return

    for file_path in tqdm(files, desc="Processing pairwise alignments"):
        try:
            process_one_file(file_path, args.dataset, args.locus, args.data_root, args.output_root)
        except Exception as e:
            print(f"[ERROR] Failed for {file_path.name}: {e}")


if __name__ == "__main__":
    main()