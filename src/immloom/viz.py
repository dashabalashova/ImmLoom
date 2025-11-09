"""viz.py

Visualization helpers for ImmLoom.
"""
from typing import Dict, Any
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import seaborn as sns, random
import pandas as pd

__all__ = ["plot_segments", "plot_necklace_diagonal"]


def plot_segments(*dfs):
    """
    Plot one or more segment DataFrames side-by-side. Color by strand.
    """
    xs, ys = [], []
    for df in dfs:
        xs.extend(df.x1.tolist())
        xs.extend(df.x2.tolist())
        ys.extend(df.y1.tolist())
        ys.extend(df.y2.tolist())
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    n = len(dfs)
    fig, axes = plt.subplots(1, n, sharex=True, sharey=True, figsize=(16, 16/max(1,n)))
    if n == 1:
        axes = [axes]
    for ax, df in zip(axes, dfs):
        for _, row in df.iterrows():
            color = 'black' if row['strand'] == '+' else 'red'
            ax.plot([row['x1'], row['x2']], [row['y1'], row['y2']], linewidth=2, color=color)
        ax.set_aspect('equal', 'box')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_max, y_min)
    plt.tight_layout()
    plt.show()

def plot_necklace_diagonal(ax, df: pd.DataFrame, info_text: str = None, title: str = None, colors: Dict[Any, Any] = None,
                           scale: float = 1.37, aspect: float = 0.8, rotate_deg: float = 45,
                           palette_name: str = 'Spectral', draw_lines: bool = True,
                           edge_color: str = 'red', edge_width: float = 0.1):
    blocks = sorted(df['block_id'].unique())
    if colors is None:
        palette = sns.color_palette(palette_name, len(blocks))
        random.seed(0)
        random.shuffle(palette)
        colors = {b: palette[i] for i, b in enumerate(blocks)}
    ax.set_title(title or '')
    for _, r in df.sort_values(['block_id', 'd1']).iterrows():
        b = r['block_id']
        x_center = (r['d1'] + r['d2']) / 2.0
        y_center = x_center
        width = max((r['d2'] - r['d1']) * scale, 0.1)
        height = width * aspect
        diamond = np.array([
            [0.0,  height/2.0],
            [width/2.0, 0.0],
            [0.0, -height/2.0],
            [-width/2.0, 0.0]])
        theta = np.deg2rad(rotate_deg)
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta),  np.cos(theta)]])
        rot = diamond.dot(R.T)
        rot[:, 0] += x_center
        rot[:, 1] += y_center
        patch = Polygon(rot, closed=True, edgecolor=edge_color, facecolor=colors[b],
                        alpha=0.9, linewidth=edge_width)
        ax.add_patch(patch)
        if draw_lines:
            ax.plot([r['d1'], r['d2']], [r['d1'], r['d2']], color=colors[b], alpha=0.25, linewidth=1)
    ax.set_aspect('equal', adjustable='datalim')
    xmin = df['d1'].min() - (df['d2'] - df['d1']).max()
    xmax = df['d2'].max() + (df['d2'] - df['d1']).max()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmax, xmin)
    handles = [mpatches.Patch(color=colors[b], label=str(b)) for b in blocks]
    ax.legend(handles=handles, title='block_id', bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.set_aspect('equal', 'box')
    return colors
