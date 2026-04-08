import numpy as np
import pandas as pd


def segment_length(x1: float, y1: float, x2: float, y2: float) -> float:
    return float(np.hypot(x2 - x1, y2 - y1))


def point_on_segment(x1: float, y1: float, x2: float, y2: float, t: float):
    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    return x, y


def get_split_params_for_segment(
    row: pd.Series,
    sx_list,
    sy_list,
    tol: float = 1e-9,
):
    """Return sorted t in [0, 1] where the segment intersects split lines."""
    x1, x2, y1, y2 = row["x1"], row["x2"], row["y1"], row["y2"]
    dx = x2 - x1
    dy = y2 - y1

    ts = [0.0, 1.0]

    if abs(dx) > tol:
        for sx in sx_list:
            t = (sx - x1) / dx
            if tol < t < 1.0 - tol:
                ts.append(t)

    if abs(dy) > tol:
        for sy in sy_list:
            t = (sy - y1) / dy
            if tol < t < 1.0 - tol:
                ts.append(t)

    ts = sorted(ts)
    unique_ts = []
    for t in ts:
        if not unique_ts or abs(t - unique_ts[-1]) > tol:
            unique_ts.append(t)

    return unique_ts


def split_one_segment(
    row: pd.Series,
    sx_list,
    sy_list,
    gap: float = 1.0,
    min_length: float = 1e-6,
):
    """Split one segment into subsegments."""
    x1, x2, y1, y2 = row["x1"], row["x2"], row["y1"], row["y2"]
    strand = row["strand"]

    total_len = segment_length(x1, y1, x2, y2)
    if total_len == 0:
        return []

    ts = get_split_params_for_segment(row, sx_list, sy_list)
    gap_t = gap / total_len

    pieces = []
    n_intervals = len(ts) - 1
    for i in range(n_intervals):
        ta = ts[i]
        tb = ts[i + 1]

        t_start = ta + (gap_t / 2 if i > 0 else 0.0)
        t_end = tb - (gap_t / 2 if i < n_intervals - 1 else 0.0)

        if t_end <= t_start:
            continue

        px1, py1 = point_on_segment(x1, y1, x2, y2, t_start)
        px2, py2 = point_on_segment(x1, y1, x2, y2, t_end)
        piece_length = segment_length(px1, py1, px2, py2)

        if piece_length < min_length:
            continue

        pieces.append(
            {
                "parent_id": row.name,
                "x1": px1,
                "x2": px2,
                "y1": py1,
                "y2": py2,
                "strand": strand,
                "length": piece_length,
            }
        )

    return pieces


def split_segments(
    df: pd.DataFrame,
    sx_list,
    sy_list,
    gap: float = 1.0,
    min_length: float = 1e-6,
) -> pd.DataFrame:
    """Split all segments from a DataFrame."""
    all_pieces = []

    for _, row in df.iterrows():
        pieces = split_one_segment(
            row,
            sx_list=sx_list,
            sy_list=sy_list,
            gap=gap,
            min_length=min_length,
        )
        all_pieces.extend(pieces)

    result = pd.DataFrame(all_pieces)
    if not result.empty:
        result["strand"] = np.where(result["y2"] >= result["y1"], "+", "-")
    return result
