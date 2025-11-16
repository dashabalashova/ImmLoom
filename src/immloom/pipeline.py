#!/usr/bin/env python3
"""run_pipeline.py

Run the ImmLoom pipeline (the notebook workflow) from the command line using the compact package.

Usage example:
    python run_pipeline.py --input-path /path/to --input-file sample.tsv --out-dir ./out --no-show

The script will:
  - read the input TSV with segm_preprocess (expects columns like start1,end1,start2+,end2+,strand2,id%)
  - filter segments
  - merge nearby segments
  - split segments
  - find projections and reversed intervals
  - reflect segments
  - save intermediate CSVs and (optionally) plots

Note: this script adds the package parent (/mnt/data/src) to sys.path so the local package can be imported.
"""

import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import networkx as nx
import sys
import pandas as pd
from pathlib import Path
import random


def run_pipeline(input_path: Path, 
                 input_file: str, 
                 out_dir: Path, 
                 #pi_min = 80.0, 
                 #segm_length_min = 3000, 
                 #dist_max = 3000, 
                 #block_length_min = 3000, 
                 #inv_plot=False
                 init_segm_length_min: float,
                 pi_min: float, 
                 segm_length_min: int, 
                 dist_max: int, 
                 block_length_min: int, 
                 inv_plot: bool):
    
    repo_src = Path(__file__).resolve().parents[1]
    if str(repo_src) not in sys.path:
        sys.path.insert(0, str(repo_src))
    from immloom import segm_preprocess, segm_filter, merge_segments, \
    split_segments, find_proj, find_reversed_intervals, reflect_segments, \
    plot_segments, find_overlap, plot_necklace_diagonal, segm_symmetric, \
    find_containing_projids

    ANIMAL_EMOJIS = ["🐵","🐶","🐱","🦊","🐼","🐯","🦁","🐸","🐨","🐰","🦄","🐷","🐮","🐔","🐧","🐢","🦉","🦋","🐙","🦖"]
    emoji = random.choice(ANIMAL_EMOJIS)

    input_path = Path(input_path)
    input_file = input_file
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"{emoji}  Run: {input_file}\n")
    
    df = segm_preprocess(input_path, input_file)
    print(f"{emoji}  Number of segments: {df.shape[0]}\n")

    df_filtered = segm_filter(df, length=init_segm_length_min, pi=pi_min)
    print(f"{emoji}  Number of segments after filtering: {df_filtered.shape[0]}\n")
    df_symmetric = segm_symmetric(df_filtered)
    df_symmetric = df_symmetric[df_symmetric.x1 >= df_symmetric.y1]

    df_merged = merge_segments(df_symmetric, dist_max=dist_max)
    dfn = df_merged[df_merged.strand=='-']
    split_points = sorted(list(set(list(dfn.x1)+list(dfn.x2)+list(dfn.y1)+list(dfn.y2))))
    df_split = split_segments(df_merged, split_points=split_points, length_min=segm_length_min)
    
    if df_split[df_split.strand=='-'].shape[0]>0:
        df6, df7 = find_proj(df_split)
        split_points = list(df7.d1)+list(df7.d2)
        df8 = split_segments(df_split, split_points = split_points, length_min=segm_length_min)
        res = df8.apply(find_containing_projids, axis=1, df7=df7)
        df8_with_hits = pd.concat([df8.reset_index(drop=True), res], axis=1)
        df8_in = df8_with_hits[df8_with_hits['which'] != 'none'].copy()
        df8_lst = list(df8_in[df8_in.which!='both'].proj_x) + list(df8_in[df8_in.which!='both'].proj_y)
        df8_lst = [x for sub in df8_lst for x in sub]
        G, inv_indices = find_reversed_intervals(df6, df7, target_nodes=df8_lst)
        df9 = reflect_segments(df8, df7, inv_indices, G, plot=inv_plot)
        assert not (df9['x1'] > df9['x2']).any() and not (df9['y1'] > df9['y2']).any(), \
            "Found rows with x1>x2 or y1>y2"
        df_inv = df7.loc[inv_indices, :][['d1', 'd2']].reset_index(drop=True)
        df_inv_path = out_dir / "dsegm_inv.tsv"
        df_inv.to_csv(df_inv_path, index=False, sep='\t')
        print(f"{emoji}  Saved inverted diagonal segments -> {df_inv_path} (rows: {len(df_inv)})\n")
    else:
        df9 = df_split

    df9_nd = df9[(df9.x1!=df9.y1) | (df9.x2!=df9.y2)]
    df9_d = df9[(df9.x1==df9.y1) & (df9.x2==df9.y2)]
    df10_nd = split_segments(df9_nd, length_min=segm_length_min)
    split_points = sorted(list(set(list(df10_nd.x1) + list(df10_nd.x2) + list(df10_nd.y1) + list(df10_nd.y2))))
    df10_d = split_segments(df9_d, split_points=split_points, length_min=segm_length_min)
    df10 = pd.concat([df10_nd, df10_d])
    dsegm1 = df10[(df10.x1==df10.y1) & (df10.x2==df10.y2)].reset_index(drop=True)
    dsegm1 = dsegm1[dsegm1.length>=block_length_min]
    dsegm1['dsegm_id'] = dsegm1.index
    dsegm1['d1'] = dsegm1['x1']
    dsegm1['d2'] = dsegm1['x2']
    dsegm1 = dsegm1.drop(columns=['x1', 'x2', 'y1', 'y2'])
    df11 = df10[(df10.x1!=df10.y1) | (df10.x2!=df10.y2)].reset_index(drop=True)

    df11.loc[:,'proj_id_x'] = df11.apply(lambda r: find_overlap(r, dsegm1, axis='x'), axis=1)
    df11.loc[:,'proj_id_y'] = df11.apply(lambda r: find_overlap(r, dsegm1, axis='y'), axis=1)
    df11.to_csv('test', sep='\t')
    edges = (df11[['proj_id_x', 'proj_id_y']].values.tolist())
    edges = [tuple(edge) for edge in edges if ((edge[0] != edge[1]) & (edge[0] != -1) & (edge[1] != -1))]
    G = nx.Graph()
    G.add_edges_from(edges)
    comp_dict = defaultdict(int)
    for comp_id, comp in enumerate(nx.connected_components(G), start=1):
        for node in comp:
            comp_dict[int(node)] = comp_id
    print(f"{emoji}  Number of blocks: {comp_id}\n")
    dsegm1['block_id'] = dsegm1['dsegm_id'].astype(int).apply(lambda x: comp_dict[x])

    df_block = dsegm1[dsegm1.block_id!=0][['block_id', 'd1', 'd2']].reset_index(drop=True)
    df_block_path = out_dir / "block.tsv"
    df_block.to_csv(df_block_path, index=False, sep='\t')
    print(f"{emoji}  Saved blocks -> {df_block_path} (rows: {len(df_block)})\n")
    
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plot_necklace_diagonal(ax, dsegm1[dsegm1.block_id!=0], title='', info_text='')
    plt.tight_layout()
    fig_path = out_dir / "block.png"
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"{emoji}  Saved block figure -> {fig_path}")
    return None
