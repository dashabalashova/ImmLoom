from pathlib import Path

import numpy as np
import pandas as pd


def add_length(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["length"] = np.hypot(df["x2"] - df["x1"], df["y2"] - df["y1"]).astype(int)
    return df


def preprocess_alignment(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    required_columns = ["start1", "end1", "start2+", "end2+", "strand2"]
    if not set(required_columns).issubset(df.columns):
        raise KeyError(f"Missing columns: {required_columns}")

    numeric_columns = ["start1", "end1", "start2+", "end2+"]
    df[numeric_columns] = (
        df[numeric_columns].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    )
    df["pi"] = (
        df.get("id%", "0").astype(str).str.rstrip("%").replace("", "0").astype(float)
    )

    df = df.rename(
        columns={
            "start1": "x1",
            "end1": "x2",
            "start2+": "y1",
            "end2+": "y2",
            "strand2": "strand",
        }
    )[["x1", "x2", "y1", "y2", "strand", "pi"]]

    minus_mask = df["strand"] == "-"
    df.loc[minus_mask, ["y1", "y2"]] = df.loc[minus_mask, ["y2", "y1"]].to_numpy()

    return add_length(df)


def filter_segments(
    df: pd.DataFrame,
    length: int = 5000,
    pi: float = 80.0,
) -> pd.DataFrame:
    return (
        df[(df["length"] >= length) & (df["pi"] >= pi)]
        .drop(columns="pi")
        .reset_index(drop=True)
    )
