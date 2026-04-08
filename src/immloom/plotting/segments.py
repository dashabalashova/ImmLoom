import matplotlib.pyplot as plt
import pandas as pd


def plot_segments(ax, df: pd.DataFrame, title: str | None = None):
    for _, row in df.iterrows():
        ax.plot([row["x1"], row["x2"]], [row["y1"], row["y2"]], linewidth=2)

    ax.set_title(title)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.invert_yaxis()


def plot_segments_with_splits(
    df: pd.DataFrame,
    result_df: pd.DataFrame,
    sx_list,
    sy_list,
    figsize=(14, 6),
):
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    ax = axes[0]
    for _, row in df.iterrows():
        ax.plot([row["x1"], row["x2"]], [row["y1"], row["y2"]], linewidth=2)
    for sx in sx_list:
        ax.axvline(sx, linestyle="--", alpha=0.6)
    for sy in sy_list:
        ax.axhline(sy, linestyle="--", alpha=0.6)
    ax.set_title("Исходные отрезки + линии разбиения")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.invert_yaxis()

    ax = axes[1]
    for _, row in result_df.iterrows():
        ax.plot([row["x1"], row["x2"]], [row["y1"], row["y2"]], linewidth=2)
    for sx in sx_list:
        ax.axvline(sx, linestyle="--", alpha=0.6)
    for sy in sy_list:
        ax.axhline(sy, linestyle="--", alpha=0.6)
    ax.set_title("Подотрезки после разбиения")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.show()
