"""segments.py

Core operations on segments: symmetry, merging, splitting, projections and reflections.
"""
from math import sqrt
from typing import List, Tuple, Dict, Set, Optional
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

__all__ = [
    "segm_symmetric", "merge_segments", "split_segments",
    "find_containing_projids", "find_overlap", "find_proj",
    "find_reversed_intervals", "reflect_segments"]


def segm_symmetric(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only symmetric segments: a segment is symmetric if its mirrored counterpart exists.
    For '+' strand, mirror is swapping x<->y. For '-' strand, mirror depends on stored orientation.
    """
    segments = set(zip(zip(df.x1, df.y1), zip(df.x2, df.y2)))
    mirrored = df.apply(
        lambda row: ((row.y1, row.x1), (row.y2, row.x2)) if row.strand == '+'
        else ((row.y2, row.x2), (row.y1, row.x1)), axis=1)
    return df[mirrored.isin(segments)].reset_index(drop=True)

def merge_segments(df: pd.DataFrame, dist_max: float = 0.0) -> pd.DataFrame:
    """
    Merge nearby segments on the same strand. Two segments are linked if the gap between
    end of one and start of the other is < dist_max (Euclidean distance in (x,y) space).
    """
    merged = []
    merged_idx = []
    df = df.reset_index(drop=True)
    df['idx_0'] = df.index
    for strand in ['+', '-']:
        dfs = df[df.strand == strand].sort_values(['x1', 'y1']).reset_index(drop=True)
        G = nx.Graph()
        for i in range(len(dfs)):
            G.add_node(i)
        for u in range(len(dfs)):
            for v in range(u+1, len(dfs)):
                a = dfs.iloc[u]
                b = dfs.iloc[v]
                if sqrt((b.x1 - a.x2)**2 + (b.y1 - a.y2)**2) < dist_max:
                    G.add_edge(u, v)
        for comp in nx.connected_components(G):
            comp = sorted(comp)
            if len(comp) >= 2:
                segs = dfs.iloc[list(comp)]
                if strand == '+':
                    x1, y1 = segs[['x1', 'y1']].min()
                    x2, y2 = segs[['x2', 'y2']].max()
                else:
                    x1, y2 = segs[['x1', 'y2']].min()
                    x2, y1 = segs[['x2', 'y1']].max()
                merged.append({
                    'x1': int(x1), 'x2': int(x2),
                    'y1': int(y1), 'y2': int(y2),
                    'strand': strand})
                merged_idx += list(segs['idx_0'])
    merged_df = pd.DataFrame(merged)
    not_merged_df = df.drop(index=merged_idx, errors='ignore')
    df1 = pd.concat([merged_df, not_merged_df], ignore_index=True).sort_values(['x1', 'y1', 'x2', 'y2'])
    df1['length'] = df1.apply(lambda x: int(sqrt((x.x2 - x.x1)**2 + (x.y2 - x.y1)**2)), axis=1)
    return df1.reset_index(drop=True)

def split_segments(df: pd.DataFrame, split_points: Optional[List[int]] = None, 
                   length_min: Optional[int] = None) -> pd.DataFrame:
    """
    Split segments by given split_points (list of coordinates). If split_points is None,
    build them from all segment endpoints, merging very close points (<5) into one.
    """
    if split_points is None:
        split_points = sorted(list(set(list(df.x1) + list(df.x2) + list(df.y1) + list(df.y2))))
    pts = split_points
    res = []
    i = 0
    n = len(pts)
    while i < n:
        j = i
        while j + 1 < n and pts[j+1] - pts[j] < 5:
            j += 1
        if i == j:
            res.append(int(pts[i]))
        else:
            res.append((pts[i] + pts[j]) // 2)
        i = j + 1
    split_points = res
    segm_lst = []
    for _, row in df.iterrows():
        x1, y1, x2, y2, strand = row.x1, row.y1, row.x2, row.y2, row.strand
        points = [(x1, y1), (x2, y2)]
        if strand == '+':
            for xt in split_points:
                if x1 < xt < x2:
                    yt = int(y1 + (y2 - y1) / (x2 - x1) * (xt - x1))
                    points.append((xt, yt))
                    points.append((xt+1, yt+1))
            for yt in split_points:
                if y1 < yt < y2:
                    xt = int(x1 + (x2 - x1) / (y2 - y1) * (yt - y1))
                    points.append((xt, yt))
                    points.append((xt+1, yt+1))
        else:
            for xt in split_points:
                if x1 < xt < x2:
                    yt = int(y1 + (y2 - y1) / (x2 - x1) * (xt - x1))
                    points.append((xt, yt+1))
                    points.append((xt+1, yt))
            for yt in split_points:
                if y2 < yt < y1:
                    xt = int(x1 + (x2 - x1) / (y2 - y1) * (yt - y1))
                    points.append((xt, yt+1))
                    points.append((xt+1, yt))
        points = sorted(points, key=lambda p: p[0])
        for (xa, ya), (xb, yb) in zip(points[:-1], points[1:]):
            segm_lst.append({'x1': int(xa), 'y1': int(ya), 'x2': int(xb), 'y2': int(yb), 'strand': strand})
    df1 = pd.DataFrame(segm_lst)
    if df1.empty:
        return df1
    df1['length'] = df1.apply(lambda x: int(sqrt((x.x2 - x.x1)**2 + (x.y2 - x.y1)**2)), axis=1)
    if length_min is not None:
        df1 = df1[df1.length >= length_min]
    return df1[['x1', 'x2', 'y1', 'y2', 'strand', 'length']].reset_index(drop=True)

def split_segments_sdir(df: pd.DataFrame, split_points: Optional[List[int]] = None, 
                   length_min: Optional[int] = None, sdir='vert') -> pd.DataFrame:
    if split_points is None:
        split_points = sorted(list(set(list(df.x1) + list(df.x2) + list(df.y1) + list(df.y2))))

    segm_lst = []
    for _, row in df.iterrows():
        x1, y1, x2, y2, strand = row.x1, row.y1, row.x2, row.y2, row.strand
        points = [(x1, y1), (x2, y2)]
        if strand == '+':
            if sdir=='vert':
                for xt in split_points:
                    if x1 < xt < x2:
                        yt = int(y1 + (y2 - y1) / (x2 - x1) * (xt - x1))
                        points.append((xt, yt))
                        points.append((xt+1, yt+1))
            else:
                for yt in split_points:
                    if y1 < yt < y2:
                        xt = int(x1 + (x2 - x1) / (y2 - y1) * (yt - y1))
                        points.append((xt, yt))
                        points.append((xt+1, yt+1))

        points = sorted(points, key=lambda p: p[0])
        for (xa, ya), (xb, yb) in zip(points[:-1], points[1:]):
            segm_lst.append({'x1': int(xa), 'y1': int(ya), 'x2': int(xb), 'y2': int(yb), 'strand': strand})
    df1 = pd.DataFrame(segm_lst)
    if df1.empty:
        return df1
    df1['length'] = df1.apply(lambda x: int(sqrt((x.x2 - x.x1)**2 + (x.y2 - x.y1)**2)), axis=1)
    if length_min is not None:
        df1 = df1[df1.length >= length_min]
    return df1[['x1', 'x2', 'y1', 'y2', 'strand', 'length']].reset_index(drop=True)

def find_containing_projids(row: pd.Series, df7: pd.DataFrame) -> pd.Series:
    """Return proj_id lists that fully contain the x- and y- intervals of a row."""
    a, b = row['x1'], row['x2']
    c, d = row['y1'], row['y2']
    mask_x = (df7['d1'] <= a) & (df7['d2'] >= b)
    mask_y = (df7['d1'] <= c) & (df7['d2'] >= d)
    proj_x = df7.loc[mask_x, 'proj_id'].tolist()
    proj_y = df7.loc[mask_y, 'proj_id'].tolist()
    if proj_x and proj_y:
        which = 'both'
    elif proj_x:
        which = 'x'
    elif proj_y:
        which = 'y'
    else:
        which = 'none'
    return pd.Series({'proj_x': proj_x, 'proj_y': proj_y, 'which': which})

def find_overlap(row: pd.Series, dsegm1: pd.DataFrame, over_min: float = 0.8, axis: str = 'x') -> int:
    """
    Find first overlapping dsegm_id in dsegm1 with overlap >= over_min fraction.
    Returns dsegm_id or -1.
    """
    if axis == 'x':
        a, b = row.x1, row.x2
    else:
        a, b = row.y1, row.y2
    if a > b:
        a, b = b, a
    for _, row2 in dsegm1.iterrows():
        c, d = row2.d1, row2.d2
        if c > d:
            c, d = d, c
        lo = max(a, c)
        hi = min(b, d)
        inter_len = max(0, hi - lo)
        if inter_len >= max(over_min * (b - a), over_min * (d - c)):
            return row2.dsegm_id
    return -1

def find_proj(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    From df (expects '-' strand for reversed segments) build df6 and df7:
      - df6: original '-' segments with assigned proj_id_x/proj_id_y where relevant
      - df7: merged intervals (d1,d2) representing projection intervals with proj_id index
    """
    def assign_projection(row, df_intervals, axis='x'):
        if axis == 'x':
            t1, t2 = row['x1'], row['x2']
        else:
            t1, t2 = row['y2'], row['y1']
        for _, proj in df_intervals.iterrows():
            if not (t2 < proj.d1 or t1 > proj.d2):
                return proj.proj_id
        return None

    df6 = df[df.strand == '-'].copy().reset_index(drop=True)
    intervals = df6[['x1', 'x2']].values.tolist() + df6[['y2', 'y1']].values.tolist()
    intervals.sort(key=lambda x: (x[0], x[1]))
    intervals_merged = []
    for start, end in intervals:
        if not intervals_merged or intervals_merged[-1][1] < start:
            intervals_merged.append([start, end])
        else:
            intervals_merged[-1][1] = max(intervals_merged[-1][1], end)
    df7 = pd.DataFrame(intervals_merged, columns=['d1', 'd2']).reset_index().rename(columns={'index': 'proj_id'})
    df6['proj_id_x'] = df6.apply(lambda r: assign_projection(r, df7, axis='x'), axis=1)
    df6['proj_id_y'] = df6.apply(lambda r: assign_projection(r, df7, axis='y'), axis=1)
    return df6, df7


def find_reversed_intervals(df6: pd.DataFrame, df7: pd.DataFrame, target_nodes: Optional[List[int]] = None) -> Tuple[nx.Graph, List[int]]:
    """
    Build a graph of proj_id_x -- proj_id_y edges from df6. Then bipartite-color it and choose
    which side to 'invert' based on weight (interval lengths) and optional target_nodes constraint.
    Returns the graph and list of proj_id indices to invert (red set in original logic).
    """
    if target_nodes is None:
        target_nodes = []
    G = nx.Graph()
    edges = df6[['proj_id_x', 'proj_id_y']].values.tolist()
    edges = [tuple(edge) for edge in edges if edge[0] != edge[1] and (edge[0] is not None) and (edge[1] is not None)]
    G.add_edges_from(edges)
    if not nx.is_bipartite(G):
        raise ValueError("Graph is not bipartite.")
    coloring = nx.bipartite.color(G)
    node_weights: Dict[int, float] = {}
    df7_idx = df7.set_index('proj_id', drop=False)
    for node in G.nodes():
        row = df7_idx.loc[node]
        d1 = float(row['d1'])
        d2 = float(row['d2'])
        w = abs(d2 - d1)
        node_weights[node] = w
    inv_indices: List[int] = []
    colors: Dict[int, str] = {}
    for comp in nx.connected_components(G):
        comp = set(comp)
        side0 = {n for n in comp if coloring.get(n, 0) == 0}
        side1 = comp - side0
        w0 = sum(node_weights.get(n, 0.0) for n in side0)
        w1 = sum(node_weights.get(n, 0.0) for n in side1)
        if w0 <= w1:
            red = side0
            green = side1
        else:
            red = side1
            green = side0
        r_flg = 0
        g_flg = 0
        for target_node in target_nodes:
            if target_node in comp:
                if target_node in red:
                    r_flg = 1
                if target_node in green:
                    g_flg = 1
        if r_flg + g_flg == 2:
            raise AssertionError(f"Node incompatibility: {sorted(target_nodes)}")
        if r_flg == 1:
            red, green = green, red
        for n in red:
            colors[n] = 'red'
        for n in green:
            colors[n] = 'green'
        inv_indices += sorted(red)
    return G, inv_indices

def reflect_segments(df: pd.DataFrame, df7: pd.DataFrame, 
                     inv_indices: List[int], G: nx.Graph, plot: bool = False) -> pd.DataFrame:
    """
    Reflect segments that fall inside intervals listed in df7 (inv_indices).
    Returns a new DataFrame with all segments reflected for proj intervals in inv_indices.
    """
    import matplotlib.pyplot as plt

    df_reflected_x = df.copy().reset_index(drop=True)
    reflected_x_idx = set()
    df_segments = df7[df7.proj_id.isin(inv_indices)].reset_index(drop=True)
    for i in range(df_segments.shape[0]):
        x1_1, x2_1 = df_segments.iloc[i][["d1", "d2"]]
        mirror_x = (x1_1 + x2_1) / 2
        df_within = df_reflected_x[(df_reflected_x['x1'] >= x1_1) & (df_reflected_x['x2'] <= x2_1)].copy()
        reflected_x_idx.update(df_within.index.tolist())
        df_within['x1'] = mirror_x - (df_within['x1'] - mirror_x)
        df_within['x2'] = mirror_x - (df_within['x2'] - mirror_x)
        df_non = df_reflected_x.drop(index=df_within.index)
        df_reflected_x = pd.concat([df_non, df_within])
    mask = (df_reflected_x['x1'] > df_reflected_x['x2']) & (df_reflected_x['y1'] > df_reflected_x['y2'])
    if mask.any():
        df_reflected_x.loc[mask, ['x1', 'x2']] = df_reflected_x.loc[mask, ['x2', 'x1']].values
        df_reflected_x.loc[mask, ['y1', 'y2']] = df_reflected_x.loc[mask, ['y2', 'y1']].values
    df_reflected_xy = df_reflected_x.copy()
    reflected_y_idx = set()
    for i in range(df_segments.shape[0]):
        x1_1, x2_1 = df_segments.iloc[i][["d1", "d2"]]
        mirror_y = (x1_1 + x2_1) / 2
        df_within = df_reflected_xy[(df_reflected_xy['y1'] >= x1_1) & (df_reflected_xy['y2'] <= x2_1)].copy()
        reflected_y_idx.update(df_within.index.tolist())
        df_within['y1'] = mirror_y - (df_within['y1'] - mirror_y)
        df_within['y2'] = mirror_y - (df_within['y2'] - mirror_y)
        df_non = df_reflected_xy.drop(index=df_within.index)
        df_reflected_xy = pd.concat([df_non, df_within])
    mask = (df_reflected_xy['x1'] > df_reflected_xy['x2']) & (df_reflected_xy['y1'] > df_reflected_xy['y2'])
    if mask.any():
        df_reflected_xy.loc[mask, ['x1', 'x2']] = df_reflected_xy.loc[mask, ['x2', 'x1']].values
        df_reflected_xy.loc[mask, ['y1', 'y2']] = df_reflected_xy.loc[mask, ['y2', 'y1']].values
    if plot:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 2, figsize=(16, 16))
        pos = nx.spring_layout(G)
        colors = []
        for n in G.nodes():
            if n in inv_indices:
                colors.append('red')
            else:
                colors.append('green')
        nx.draw(G, pos, node_color=colors, with_labels=True, ax=axs[0, 0], node_size=200, font_size=12)
        axs[0, 0].set_title("Projection graph (red=inverted)")

        def _plot_segments(ax, df_all, df_highlight_idx, title):
            for i, row in df_all.iterrows():
                color = 'green' if i in df_highlight_idx else 'black'
                linewidth = 4 if i in df_highlight_idx else 2
                ax.plot([row["x1"], row["x2"]], [row["y1"], row["y2"]], color=color, linewidth=linewidth)
            for _, seg in df_segments.iterrows():
                d1, d2 = seg["d1"], seg["d2"]
                mirror = (d1 + d2) / 2
                ax.axvline(d1, color='red', linestyle='--', linewidth=1)
                ax.axvline(d2, color='red', linestyle='--', linewidth=1)
                ax.axhline(d1, color='red', linestyle='--', linewidth=1)
                ax.axhline(d2, color='red', linestyle='--', linewidth=1)
                ax.axvline(mirror, color='red', linestyle=':', linewidth=1.5)
                ax.axhline(mirror, color='red', linestyle=':', linewidth=1.5)
            ax.set_title(title)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_aspect("equal")
            ax.invert_yaxis()
        _plot_segments(axs[0, 1], df, set(), "Original")
        _plot_segments(axs[1, 0], df_reflected_x, reflected_x_idx, "Reflected in X (highlighted)")
        _plot_segments(axs[1, 1], df_reflected_xy, reflected_y_idx, "Reflected in X & Y (highlighted)")
        plt.tight_layout()
        plt.show()

    df_reflected_xy = df_reflected_xy[['x1', 'x2', 'y1', 'y2', 'strand']].copy()
    df_reflected_xy['strand'] = '+'
    df_reflected_xy[['x1', 'x2', 'y1', 'y2']] = df_reflected_xy[['x1', 'x2', 'y1', 'y2']].astype(int)
    return df_reflected_xy

def create_and_plot_graph(dfx, dfy, df4, n1='n1', n2='n2', figsize=(8,6), plot=False):
    if 'block_id_x' not in df4.columns:
        df4['block_id_x'] = df4.apply(lambda r: find_block_id(r, dfx, 'x'), axis=1)
    if 'block_id_y' not in df4.columns:
        df4['block_id_y'] = df4.apply(lambda r: find_block_id(r, dfy, 'y'), axis=1)
    def prefixed_bid(prefix, bid):
        return f"{prefix}_{int(bid)}" if pd.notna(bid) else None
    edges = []
    for _, r in df4.iterrows():
        bx = prefixed_bid(n1, r.get('block_id_x'))
        by = prefixed_bid(n2, r.get('block_id_y'))
        if bx is not None and by is not None:
            edges.append((bx, by))
    G = nx.Graph()
    G.add_edges_from(edges)
    if plot:
        plt.figure(figsize=figsize)
        if len(G) == 0:
            plt.text(0.5, 0.5, "Graph is empty (no edges between prefixed block ids)", ha='center', va='center')
            plt.axis('off')
        else:
            left_nodes = [n for n in G.nodes if str(n).startswith(f"{n1}_")]
            right_nodes = [n for n in G.nodes if str(n).startswith(f"{n2}_")]
            if left_nodes and right_nodes:
                pos = {}
                for i, node in enumerate(left_nodes):
                    pos[node] = (-1, (i - len(left_nodes)/2))
                for i, node in enumerate(right_nodes):
                    pos[node] = (1, (i - len(right_nodes)/2))
                combined = nx.spring_layout(G, pos=pos, fixed=list(pos.keys()), seed=3)
                pos = combined
            else:
                pos = nx.spring_layout(G, seed=4)
            nx.draw(G, pos, with_labels=True, node_size=700, font_size=9)
        plt.title(f"Edges between {n1}_block_id and {n2}_block_id (from df4)")
        plt.show()
    return G