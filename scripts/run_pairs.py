#!/usr/bin/env python3
"""Generate graph edges/nodes TSV from immloom segment/block outputs.

Usage:
    python generate_graph_from_segments.py --n1 0 --n2 1 \
        --data-dir ../data/processed/dataset_01/patchwork_output/IGH/pairwise_alignments \
        --immloom-out ../immloom_out --out-dir graphs

Defaults match the arguments you provided.
"""
from pathlib import Path
import argparse
import os
import pandas as pd
import networkx as nx
from typing import Optional

# Project-specific helpers (from immloom package)
from immloom import (
    segm_preprocess,
    segm_filter,
    split_segments_sdir,
    create_and_plot_graph,
)


def find_block_id(row: pd.Series, blocks_df: pd.DataFrame, prefix: str) -> Optional[int]:
    """Find block_id in blocks_df that fully contains the segment described by row.

    Returns the first matching block_id or pd.NA if none found.
    """
    mask = (blocks_df['d1'] <= row[f'{prefix}1']) & (blocks_df['d2'] >= row[f'{prefix}2'])
    if mask.any():
        return blocks_df.loc[mask, 'block_id'].iat[0]
    return pd.NA


def build_graph_and_save(
    n1: int,
    n2: int,
    data_dir: Path,
    immloom_out: Path,
    out_dir: Path,
    pi: float = 80.0,
    min_length: int = 2000,
) -> None:
    data_dir = Path(data_dir)
    immloom_out = Path(immloom_out)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # list of sample ids based on filenames in data_dir that start with 'self'
    s_lst = sorted([x[5:-4] for x in os.listdir(data_dir) if x[:4] == 'self'])

    if n1 < 0 or n1 >= len(s_lst) or n2 < 0 or n2 >= len(s_lst):
        raise IndexError(f"n1/n2 indices out of range: {n1}, {n2} (len samples = {len(s_lst)})")

    # choose pair or self filename
    if n1 != n2:
        f_name = f'pair_{s_lst[n1]}_{s_lst[n2]}.tsv'
    else:
        f_name = f'self_{s_lst[n1]}.tsv'

    # preprocess and filter segments
    df = segm_preprocess(data_dir, f_name)
    df1_all = segm_filter(df, pi=pi, length=10000)
    df1 = df1_all[df1_all.strand == '+']

    # read immloom block files for each side
    block_x_path = immloom_out / f'self_{s_lst[n1]}' / 'block.tsv'
    block_y_path = immloom_out / f'self_{s_lst[n2]}' / 'block.tsv'

    dfx = pd.read_csv(block_x_path, sep='\t')
    dfy = pd.read_csv(block_y_path, sep='\t')

    # split segments by block boundaries on x side
    split_points = sorted(list(set(list(dfx.d1) + list(dfx.d2))))
    df2 = split_segments_sdir(df1, split_points=split_points)

    # then split by block boundaries on y side (horiz direction)
    split_points = sorted(list(set(list(dfy.d1) + list(dfy.d2))))
    df3 = split_segments_sdir(df2, split_points=split_points, sdir='horiz')

    # mark whether segment is inside blocks on each side
    df3['in_x'] = df3.apply(
        lambda r: ((dfx['d1'] <= r['x1']) & (dfx['d2'] >= r['x2'])).any(),
        axis=1,
    )
    df3['in_y'] = df3.apply(
        lambda r: ((dfy['d1'] <= r['y1']) & (dfy['d2'] >= r['y2'])).any(),
        axis=1,
    )

    # filter by minimal length
    df3 = df3[df3.length >= min_length]

    # compute block ids for segments that fall entirely within a block
    df4 = df3.copy()
    df4['block_id_x'] = df4.apply(lambda r: find_block_id(r, dfx, 'x'), axis=1)
    df4['block_id_y'] = df4.apply(lambda r: find_block_id(r, dfy, 'y'), axis=1)

    # create graph from blocks and segments; set plot=False to avoid plotting here
    G = create_and_plot_graph(dfx, dfy, df4, n1, n2, plot=False)

    # build edges dataframe from block ids if available, otherwise from networkx graph
    if 'block_id_x' in df4.columns and 'block_id_y' in df4.columns:
        edges_df = (
            df4
            .dropna(subset=['block_id_x', 'block_id_y'])
            .assign(
                source=lambda d: d['block_id_x'].apply(lambda x: f"{n1}_{int(x)}"),
                target=lambda d: d['block_id_y'].apply(lambda x: f"{n2}_{int(x)}"),
            )
            .groupby(['source', 'target'])
            .size()
            .reset_index(name='count')
        )
    else:
        edges_df = nx.to_pandas_edgelist(G).rename(columns={'source': 'source', 'target': 'target'})
        if 'weight' not in edges_df.columns and 'count' not in edges_df.columns:
            edges_df['count'] = 1

    edge_path = out_dir / f"G-{n1}-{n2}.tsv"
    edges_df.to_csv(edge_path, sep='\t', index=False)
    print(f"Saved edges -> {edge_path}")
    print(edges_df.head().to_string())

    # nodes and degrees
    nodes = list(G.nodes())
    nodes_df = pd.DataFrame({'node': nodes})

    def side_of(node: str) -> int:
        s = str(node)
        if s.startswith(f"{n1}_"):
            return n1
        if s.startswith(f"{n2}_"):
            return n2
        return -1

    nodes_df['side'] = nodes_df['node'].apply(side_of)
    deg = dict(G.degree())
    nodes_df['degree'] = nodes_df['node'].map(deg).fillna(0).astype(int)

    nodes_path = out_dir / f"G-{n1}-{n2}-nodes.tsv"
    nodes_df.to_csv(nodes_path, sep='\t', index=False)
    print(f"Saved nodes -> {nodes_path}")
    print(nodes_df.head().to_string())


def parse_args():
    p = argparse.ArgumentParser(description="Create graph TSVs from immloom segmentation/block outputs")
    p.add_argument('--n1', type=int, default=0, help='index of left sample (default: 0)')
    p.add_argument('--n2', type=int, default=1, help='index of right sample (default: 1)')
    p.add_argument('--data-dir', type=Path, default=Path('data/processed/dataset_01/patchwork_output/IGH/pairwise_alignments'), help='directory with pairwise alignment TSV files')
    p.add_argument('--immloom-out', type=Path, default=Path('immloom_out'), help='directory with immloom outputs (self_*/block.tsv)')
    p.add_argument('--out-dir', type=Path, default=Path('immloom_out/graphs'), help='output directory for graph TSVs')
    p.add_argument('--pi', type=float, default=80.0, help='percent identity threshold passed to segm_filter')
    p.add_argument('--min-length', type=int, default=2000, help='minimum segment length to keep')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    build_graph_and_save(
        n1=args.n1,
        n2=args.n2,
        data_dir=args.data_dir,
        immloom_out=args.immloom_out,
        out_dir=args.out_dir,
        pi=args.pi,
        min_length=args.min_length,
    )
