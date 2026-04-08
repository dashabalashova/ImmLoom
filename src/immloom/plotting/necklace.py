import random
from typing import Any, Dict

import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse, Polygon


def _build_block_colors(df: pd.DataFrame, colors: Dict[Any, Any] | None, palette_name: str):
    blocks = sorted(df["block_id"].dropna().unique())
    if colors is None:
        palette = sns.color_palette(palette_name, max(len(blocks), 1))
        random.seed(0)
        random.shuffle(palette)
        colors = {block_id: palette[i] for i, block_id in enumerate(blocks)}
    return blocks, colors


def plot_necklace_diagonal(
    ax,
    df: pd.DataFrame,
    info_text: str | None = None,
    title: str | None = None,
    colors: Dict[Any, Any] | None = None,
    scale: float = 1.37,
    aspect: float = 0.8,
    rotate_deg: float = 45,
    palette_name: str = "Spectral",
    draw_lines: bool = True,
    edge_color: str = "red",
    edge_width: float = 0.1,
    show_legend: bool = True,
):
    df = df.copy()
    blocks, colors = _build_block_colors(df, colors, palette_name)
    ax.set_title(title or "")

    for _, row in df.sort_values(["block_id", "d1"]).iterrows():
        block_id = row["block_id"]
        x_center = (row["d1"] + row["d2"]) / 2.0
        y_center = x_center
        width = max((row["d2"] - row["d1"]) * scale, 0.1)
        height = width * aspect

        diamond = np.array([
            [0.0, height / 2.0],
            [width / 2.0, 0.0],
            [0.0, -height / 2.0],
            [-width / 2.0, 0.0],
        ])
        theta = np.deg2rad(rotate_deg)
        rotation = np.array(
            [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
        )
        rotated = diamond @ rotation.T
        rotated[:, 0] += x_center
        rotated[:, 1] += y_center

        ax.add_patch(
            Polygon(
                rotated,
                closed=True,
                edgecolor=edge_color,
                facecolor=colors[block_id],
                alpha=0.9,
                linewidth=edge_width,
            )
        )

        if draw_lines:
            ax.plot(
                [row["d1"], row["d2"]],
                [row["d1"], row["d2"]],
                color=colors[block_id],
                alpha=0.25,
                linewidth=1,
            )

    if len(df) > 0:
        span = (df["d2"] - df["d1"]).max()
        xmin = df["d1"].min() - span
        xmax = df["d2"].max() + span
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(xmax, xmin)

    ax.set_aspect("equal", "box")

    if show_legend and len(blocks) > 0:
        handles = [mpatches.Patch(color=colors[b], label=str(b)) for b in blocks]
        ax.legend(handles=handles, title="block_id", bbox_to_anchor=(1.02, 1), loc="upper left")

    if info_text is not None:
        ax.text(1.02, 0.5, info_text, transform=ax.transAxes, va="center", ha="left", fontsize=10)

    return colors


def plot_necklace_diagonal_strand(
    ax,
    df: pd.DataFrame,
    info_text: str | None = None,
    title: str | None = None,
    colors: Dict[Any, Any] | None = None,
    scale: float = 1.37,
    aspect: float = 0.8,
    rotate_deg: float = 45,
    palette_name: str = "Spectral",
    draw_lines: bool = True,
    edge_color: str = "red",
    edge_width: float = 0.1,
    show_legend: bool = True,
    oval_scale_major: float = 25.0,
    oval_lw: float = 0.2,
    oval_angle: float = 135,
):
    df = df.copy()
    blocks, colors = _build_block_colors(df, colors, palette_name)
    ax.set_title(title or "")
    minus_ovals = []

    for _, row in df.sort_values(["block_id", "d1"]).iterrows():
        block_id = row["block_id"]
        x_center = (row["d1"] + row["d2"]) / 2.0
        y_center = x_center
        base_len = max(row["d2"] - row["d1"], 0.1)
        width = base_len * scale
        height = width * aspect

        diamond = np.array([
            [0.0, height / 2.0],
            [width / 2.0, 0.0],
            [0.0, -height / 2.0],
            [-width / 2.0, 0.0],
        ])
        theta = np.deg2rad(rotate_deg)
        rotation = np.array(
            [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
        )
        rotated = diamond @ rotation.T
        rotated[:, 0] += x_center
        rotated[:, 1] += y_center

        ax.add_patch(
            Polygon(
                rotated,
                closed=True,
                edgecolor=edge_color,
                facecolor=colors[block_id],
                alpha=0.9,
                linewidth=edge_width,
                zorder=2,
            )
        )

        if draw_lines:
            ax.plot(
                [row["d1"], row["d2"]],
                [row["d1"], row["d2"]],
                color=colors[block_id],
                alpha=0.25,
                linewidth=1,
                zorder=1,
            )

        if row.get("strand", "+") == "-":
            minus_ovals.append({"x": x_center, "y": y_center, "len": base_len})

    if len(df) > 0:
        span = (df["d2"] - df["d1"]).max()
        xmin = df["d1"].min() - span
        xmax = df["d2"].max() + span
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(xmax, xmin)

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xmin_ax, xmax_ax = min(x0, x1), max(x0, x1)
    ymin_ax, ymax_ax = min(y0, y1), max(y0, y1)

    theta = np.deg2rad(oval_angle)
    c = abs(np.cos(theta))
    s = abs(np.sin(theta))

    for oval in minus_ovals:
        xc, yc = oval["x"], oval["y"]
        minor = oval["len"]
        major = oval["len"] * oval_scale_major

        b = minor / 2.0
        a = major / 2.0

        avail_x = min(xc - xmin_ax, xmax_ax - xc)
        avail_y = min(yc - ymin_ax, ymax_ax - yc)
        if avail_x <= 0 or avail_y <= 0:
            continue

        bx = b * s
        by = b * c
        max_a_x = (avail_x - bx) / c if c > 1e-12 else np.inf
        max_a_y = (avail_y - by) / s if s > 1e-12 else np.inf
        max_a = min(max_a_x, max_a_y)

        a = min(a, max_a * 0.95)
        if a <= 0:
            continue

        ellipse = Ellipse(
            (xc, yc),
            width=2 * a,
            height=2 * b,
            angle=oval_angle,
            facecolor="none",
            edgecolor="red",
            linewidth=oval_lw,
            zorder=5,
        )
        ellipse.set_clip_path(ax.patch)
        ax.add_patch(ellipse)

    ax.set_aspect("equal", "box")

    if show_legend and len(blocks) > 0:
        handles = [mpatches.Patch(color=colors[b], label=str(b)) for b in blocks]
        if minus_ovals:
            handles.append(mpatches.Patch(facecolor="none", edgecolor="red", label="reversed"))
        ax.legend(handles=handles, title="block_id", bbox_to_anchor=(1.02, 1), loc="upper left")

    if info_text is not None:
        ax.text(1.02, 0.5, info_text, transform=ax.transAxes, va="center", ha="left", fontsize=10)

    return colors
