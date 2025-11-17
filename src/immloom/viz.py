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
from pathlib import Path

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

def plot_segments_c(*dfs, seed=0):
    """
    Plot one or more segment DataFrames side-by-side. Color by strand.
    For strand '+' each segment gets a random color; for strand '-' color is red.
    Optional: seed for reproducible random colors.
    """
    if seed is not None:
        random.seed(seed)

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
            if row['strand'] == '+':
                # random RGB tuple, each component in [0,1)
                color = tuple(random.random() for _ in range(3))
            else:
                color = 'red'
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

def plot_segments_color(df: pd.DataFrame, *,
                        save_path: Path | str | None = None,
                        show: bool = False,
                        title: str | None = None,
                        linewidth: float = 2.0,
                        alpha: float = 0.9,
                        figsize=(10, 8)):
    """
    Plot segments colored by membership flags 'in_x' and 'in_y'.

    Color scheme:
      - red    : neither in_x nor in_y
      - blue   : in_x only
      - yellow : in_y only
      - green  : both in_x and in_y
    """
    if df is None or len(df) == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_title(title or "No segments to plot")
        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, bbox_inches="tight", dpi=150)
            plt.close(fig)
        elif show:
            plt.show()
        return fig, ax

    required = ['x1', 'x2', 'y1', 'y2']
    for c in required:
        if c not in df.columns:
            raise ValueError(f"DataFrame must contain column '{c}'")

    df_plot = df.dropna(subset=required).copy()

    # Normalize flags to booleans if present
    def to_bool_col(df_local, col):
        if col not in df_local.columns:
            return pd.Series([False] * len(df_local), index=df_local.index)
        return df_local[col].astype(bool)

    in_x = to_bool_col(df_plot, 'in_x')
    in_y = to_bool_col(df_plot, 'in_y')

    both_mask = in_x & in_y
    only_x_mask = in_x & ~in_y
    only_y_mask = in_y & ~in_x
    none_mask = ~in_x & ~in_y

    # extents for plotting limits (safe)
    x_vals = df_plot[['x1', 'x2']].values
    y_vals = df_plot[['y1', 'y2']].values
    x_min = int(np.nanmin(x_vals))
    x_max = int(np.nanmax(x_vals))
    y_min = int(np.nanmin(y_vals))
    y_max = int(np.nanmax(y_vals))

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # plot order: none (red) → only_x (blue) → only_y (yellow) → both (green)
    color_specs = [
        (none_mask,  'red',    'none (neither)'),
        (only_x_mask,'blue',   'only x'),
        (only_y_mask,'yellow', 'only y'),
        (both_mask,  'green',  'both x & y'),
    ]

    for mask, color, label in color_specs:
        df_subset = df_plot[mask]
        for _, row in df_subset.iterrows():
            ax.plot([row['x1'], row['x2']], [row['y1'], row['y2']],
                    linewidth=linewidth, color=color, alpha=alpha)
    # Build legend (one handle per category)
    from matplotlib.lines import Line2D
    legend_elems = [Line2D([0], [0], color=c, lw=3, label=l) for (_, c, l) in color_specs]
    ax.legend(handles=legend_elems, loc='best')

    ax.set_title(title or '')
    ax.set_aspect('equal', 'box')
    xpad = max(1, (x_max - x_min) * 0.02)
    ypad = max(1, (y_max - y_min) * 0.02)
    ax.set_xlim(x_min - xpad, x_max + xpad)
    ax.set_ylim(y_max + ypad, y_min - ypad)  # invert y by swapping if needed

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    plt.tight_layout()

    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight", dpi=150)
        plt.close(fig)
    else:
        if show:
            plt.show()

    return fig, ax
