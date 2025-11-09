"""io_utils.py

Preprocessing utilities for ImmLoom.
"""
from pathlib import Path
from math import sqrt
import pandas as pd
from typing import Tuple

__all__ = ["segm_preprocess", "segm_filter"]


def segm_preprocess(f_path: Path, f_name: str) -> pd.DataFrame:
    """
    Read a TSV-like file and normalize to columns:
    ['x1','x2','y1','y2','strand','length','pi']

    Notes:
      - expects columns 'start1','end1','start2+','end2+','strand2','id%'
        (these come from lastz-like output that the user provided).
      - 'id%' values like '80%' are converted to float 80.0.
    """
    p = Path(f_path) / f_name
    df = pd.read_csv(p, sep='\t', dtype=str)
    for col in ['start1','end1','start2+','end2+']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
        else:
            raise KeyError(f"Expected column {col} in input file")
    df['length'] = df.apply(lambda x: int(sqrt((int(x['end1'])-int(x['start1']))**2 + \
                                               (int(x['end2+'])-int(x['start2+']))**2)), axis=1)
    if 'id%' in df.columns:
        df['pi'] = df['id%'].astype(str).str.rstrip('%').replace('', '0').astype(float)
    else:
        df['pi'] = 0.0
    df = df[['start1', 'end1', 'start2+', 'end2+', 'strand2', 'length', 'pi']].copy()
    df.columns = ['x1', 'x2', 'y1', 'y2', 'strand', 'length', 'pi']
    mask = df['strand'] == '-'
    if mask.any():
        df.loc[mask, ['y1','y2']] = df.loc[mask, ['y2','y1']].to_numpy()
    return df


def segm_filter(df: pd.DataFrame, length: int = 3000, pi: float = 80.0) -> pd.DataFrame:
    """
    Filter segments by minimum length and percent identity (pi).
    Returns DataFrame without the 'pi' column.
    """
    df2 = df[(df.length >= length) & (df.pi >= pi)].copy()
    return df2.drop(columns=['pi'])
