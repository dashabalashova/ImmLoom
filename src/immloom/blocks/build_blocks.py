#!/usr/bin/env python3
"""Build block tables and plots from self alignments.

Example:
    python -m immloom.blocks.self_blocks --dataset dataset-01 --locus IGH
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from immloom.io.alignment import filter_segments, preprocess_alignment
from immloom.geometry.segments import split_segments
from immloom.plotting.necklace import plot_necklace_diagonal_strand


SELF_MIN_LEN = 10_000
MINUS_MIN_LEN = 5_000
MIN_PI = 80.0
SPLIT_GAP = 0.8
FINAL_MIN_LEN = 2_500


def find_diag_node(val: float, segm_diag: pd.DataFrame) -> int | None:
    hits = segm_diag[(segm_diag["d1"] <= val) & (val <= segm_diag["d2"])]
    if hits.empty:
        return None

    hits = hits.assign(diag_len=hits["d2"] - hits["d1"])
    return int(hits.sort_values("diag_len").iloc[0]["node_id"])


def split_and_filter_segments(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    for _ in range(2):
        coords = list(result.x1) + list(result.x2) + list(result.y1) + list(result.y2)
        result = split_segments(result, sx_list=coords, sy_list=coords, gap=SPLIT_GAP)
        result = result[result.length >= FINAL_MIN_LEN].reset_index(drop=True)
    return result


def build_component_graph(segm_diag: pd.DataFrame, segm_ndiag: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()

    for _, row in segm_diag.iterrows():
        graph.add_node(
            int(row["node_id"]),
            x1=row["x1"],
            x2=row["x2"],
            y1=row["y1"],
            y2=row["y2"],
            d1=row["d1"],
            d2=row["d2"],
            length=row["length"],
        )

    for idx, row in segm_ndiag.iterrows():
        x_mid = 0.5 * (row["x1"] + row["x2"])
        y_mid = 0.5 * (row["y1"] + row["y2"])

        node_x = find_diag_node(x_mid, segm_diag)
        node_y = find_diag_node(y_mid, segm_diag)

        if node_x is None or node_y is None or node_x == node_y:
            continue

        graph.add_edge(
            int(node_x),
            int(node_y),
            segm_id=idx,
            x1=row["x1"],
            x2=row["x2"],
            y1=row["y1"],
            y2=row["y2"],
            strand=row.get("strand"),
            length=row["length"],
        )

    return graph


def assign_blocks(segm_diag: pd.DataFrame, graph: nx.Graph) -> pd.DataFrame:
    components = list(nx.connected_components(graph))
    node_to_block: dict[int, int] = {}

    for block_id, comp in enumerate(components):
        for node in comp:
            node_to_block[int(node)] = block_id

    segm_diag = segm_diag.copy()
    segm_diag["block_id"] = segm_diag["node_id"].map(node_to_block)

    comp_sizes = {i: len(comp) for i, comp in enumerate(components)}
    non_singletons = segm_diag[segm_diag["block_id"].map(comp_sizes).fillna(0) > 1].copy()
    non_singletons["block_id"] = non_singletons["block_id"].factorize()[0] + 1
    return non_singletons.reset_index(drop=True)


def assign_minus_strand_nodes(
    segm_diag_non_singletons: pd.DataFrame,
    segm_ndiag: pd.DataFrame,
    full_diag: pd.DataFrame,
    sample_name: str,
) -> list[int]:
    node_to_block = dict(
        zip(segm_diag_non_singletons["node_id"], segm_diag_non_singletons["block_id"])
    )

    segm_ndiag_block = segm_ndiag.copy()
    segm_ndiag_block["x_mid"] = 0.5 * (segm_ndiag_block["x1"] + segm_ndiag_block["x2"])
    segm_ndiag_block["y_mid"] = 0.5 * (segm_ndiag_block["y1"] + segm_ndiag_block["y2"])
    segm_ndiag_block["block_id"] = (
        segm_ndiag_block["x_mid"].apply(lambda v: find_diag_node(v, full_diag)).map(node_to_block)
    )

    minus_blocks = sorted(
        set(segm_ndiag_block.loc[segm_ndiag_block["strand"] == "-", "block_id"].dropna())
    )

    node_id_minus_lst: list[int] = []

    for block_id in minus_blocks:
        q1 = segm_ndiag_block[segm_ndiag_block["block_id"] == block_id].copy()
        q2 = segm_diag_non_singletons[segm_diag_non_singletons["block_id"] == block_id].copy()

        q2 = q2.reset_index(drop=True)
        q2["node_id_local"] = q2.index
        q2["d1"] = q2[["d1", "d2"]].min(axis=1)
        q2["d2"] = q2[["d1", "d2"]].max(axis=1)

        def find_q2_node(val: float) -> int | None:
            hits = q2[(q2["d1"] <= val) & (val <= q2["d2"])]
            if hits.empty:
                return None
            hits = hits.assign(diag_len=hits["d2"] - hits["d1"])
            return int(hits.sort_values("diag_len").iloc[0]["node_id_local"])

        g_block = nx.Graph()
        for _, row in q2.iterrows():
            g_block.add_node(
                int(row["node_id_local"]),
                node_id=int(row["node_id"]),
                block_id=int(row["block_id"]),
                x1=row["x1"],
                x2=row["x2"],
                y1=row["y1"],
                y2=row["y2"],
                d1=row["d1"],
                d2=row["d2"],
                length=row["length"],
            )

        q1["node_x"] = q1["x_mid"].apply(find_q2_node)
        q1["node_y"] = q1["y_mid"].apply(find_q2_node)

        for idx, row in q1.iterrows():
            node_x = row["node_x"]
            node_y = row["node_y"]
            if node_x is None or node_y is None or node_x == node_y:
                continue
            if row.strand == '-':
                g_block.add_edge(
                    int(node_x),
                    int(node_y),
                    segm_id=idx,
                    x1=row["x1"],
                    x2=row["x2"],
                    y1=row["y1"],
                    y2=row["y2"],
                    x_mid=row["x_mid"],
                    y_mid=row["y_mid"],
                    length=row["length"],
                )

        if not nx.is_bipartite(g_block):
            print(f"Warning: graph is not bipartite for block {block_id} (sample={sample_name})")
            continue

        part1, part2 = nx.bipartite.sets(g_block)
        part1 = sorted(part1)
        part2 = sorted(part2)

        if len(part1) < len(part2):
            smaller_part = part1
        elif len(part2) < len(part1):
            smaller_part = part2
        else:
            smaller_part = part1 if [str(x) for x in part1] < [str(x) for x in part2] else part2

        minus_node_ids = (
            q2.loc[q2["node_id_local"].isin(smaller_part), "node_id"]
            .astype(int)
            .tolist()
        )
        node_id_minus_lst.extend(minus_node_ids)

    return node_id_minus_lst


def process_file(file_path: Path, output_root: Path) -> None:
    df = preprocess_alignment(file_path)

    df_minus = filter_segments(
        df[df["strand"] == "-"].reset_index(drop=True),
        length=MINUS_MIN_LEN,
        pi=MIN_PI,
    )
    df_plus = filter_segments(
        df[df["strand"] != "-"].reset_index(drop=True),
        length=SELF_MIN_LEN,
        pi=MIN_PI,
    )
    df = pd.concat([df_plus, df_minus], ignore_index=True)

    if df.empty:
        print(f"Skip {file_path.name}: no segments after filtering")
        return

    result_df = split_and_filter_segments(df)
    if result_df.empty:
        print(f"Skip {file_path.name}: no segments after splitting/filtering")
        return

    segm_diag = result_df[np.isclose(result_df["x1"], result_df["y1"])].copy().reset_index(drop=True)
    segm_ndiag = result_df[~np.isclose(result_df["x1"], result_df["y1"])].copy().reset_index(drop=True)

    if segm_diag.empty:
        print(f"Skip {file_path.name}: no diagonal segments")
        return

    segm_diag["node_id"] = segm_diag.index.astype(int)
    segm_diag["d1"] = segm_diag[["x1", "x2"]].min(axis=1)
    segm_diag["d2"] = segm_diag[["x1", "x2"]].max(axis=1)

    graph = build_component_graph(segm_diag, segm_ndiag)
    segm_diag_non_singletons = assign_blocks(segm_diag, graph)

    if segm_diag_non_singletons.empty:
        print(f"No non-singleton components found for {file_path.name}")
        return

    sample_name = file_path.stem.split("-", 1)[1] if "-" in file_path.stem else file_path.stem
    minus_nodes = assign_minus_strand_nodes(segm_diag_non_singletons, segm_ndiag, segm_diag, sample_name)
    minus_nodes_set = set(minus_nodes)

    segm_diag_non_singletons = segm_diag_non_singletons.copy()
    segm_diag_non_singletons["strand"] = np.where(
        segm_diag_non_singletons["node_id"].isin(minus_nodes_set),
        "-",
        "+",
    )
    segm_diag_non_singletons["forward"] = ~segm_diag_non_singletons["node_id"].isin(minus_nodes_set)
    segm_diag_non_singletons["d1"] = np.ceil(segm_diag_non_singletons["d1"]).astype(int)
    segm_diag_non_singletons["d2"] = np.floor(segm_diag_non_singletons["d2"]).astype(int)

    tables_dir = output_root / "tables"
    figures_dir = output_root / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 8))
    plot_necklace_diagonal_strand(
        ax,
        segm_diag_non_singletons,
        title=f"Blocks from self alignments\n{sample_name}",
    )
    plt.tight_layout()
    fig.savefig(figures_dir / f"blocks_{sample_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    segm_diag_non_singletons[["d1", "d2", "block_id", "forward"]].to_csv(
        tables_dir / f"blocks_{sample_name}.tsv",
        index=False,
        sep="\t",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build block tables and plots from self pairwise alignments."
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. dataset-01")
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


def main() -> None:
    args = parse_args()

    alignments_dir = (
        Path(args.data_root)
        / args.dataset
        / args.locus
        / "pwp"
        / "pairwise_alignments_forward"
    )
    output_root = Path(args.output_root) / args.dataset / args.locus / "blocks"

    if not alignments_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {alignments_dir}")

    files = sorted(
        p for p in alignments_dir.iterdir()
        if p.is_file() and p.name.startswith("self")
    )
    if not files:
        raise FileNotFoundError(f"No self* alignment files found in: {alignments_dir}")

    output_root.mkdir(parents=True, exist_ok=True)

    for file_path in tqdm(files, desc="Processing self alignments"):
        process_file(file_path, output_root)


if __name__ == "__main__":
    main()