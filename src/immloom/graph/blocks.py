import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


def interval_overlap(a1, a2, b1, b2):
    """Check whether two closed intervals overlap."""
    return max(a1, b1) <= min(a2, b2)


def build_bipartite_block_graph(
    dfx: pd.DataFrame,
    dfy: pd.DataFrame,
    result_df: pd.DataFrame,
):
    """Build a bipartite graph between x-blocks and y-blocks."""
    graph = nx.Graph()

    for _, row in dfx.iterrows():
        node = ("x", row["block_id"])
        graph.add_node(
            node,
            bipartite=0,
            axis="x",
            block_id=row["block_id"],
            d1=row["d1"],
            d2=row["d2"],
        )

    for _, row in dfy.iterrows():
        node = ("y", row["block_id"])
        graph.add_node(
            node,
            bipartite=1,
            axis="y",
            block_id=row["block_id"],
            d1=row["d1"],
            d2=row["d2"],
        )

    dfx_ = dfx.copy()
    dfy_ = dfy.copy()

    for _, seg in result_df.iterrows():
        sx1, sx2 = sorted((seg["x1"], seg["x2"]))
        sy1, sy2 = sorted((seg["y1"], seg["y2"]))

        x_hits = dfx_[np.maximum(dfx_["d1"], sx1) <= np.minimum(dfx_["d2"], sx2)]
        y_hits = dfy_[np.maximum(dfy_["d1"], sy1) <= np.minimum(dfy_["d2"], sy2)]

        if x_hits.empty or y_hits.empty:
            continue

        parent_id = seg["parent_id"] if "parent_id" in seg.index else None

        for _, xb in x_hits.iterrows():
            for _, yb in y_hits.iterrows():
                u = ("x", xb["block_id"])
                v = ("y", yb["block_id"])

                if graph.has_edge(u, v):
                    graph[u][v]["weight"] += 1
                    if parent_id is not None:
                        graph[u][v]["parent_ids"].append(parent_id)
                else:
                    graph.add_edge(
                        u,
                        v,
                        weight=1,
                        parent_ids=[] if parent_id is None else [parent_id],
                    )

    return graph


def draw_bipartite_graph(graph, figsize=(10, 8)):
    x_nodes = [n for n, d in graph.nodes(data=True) if d["axis"] == "x"]
    y_nodes = [n for n, d in graph.nodes(data=True) if d["axis"] == "y"]

    pos = {}
    for i, node in enumerate(sorted(x_nodes, key=lambda z: z[1])):
        pos[node] = (0, -i)
    for i, node in enumerate(sorted(y_nodes, key=lambda z: z[1])):
        pos[node] = (1, -i)

    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(graph, pos, nodelist=x_nodes, node_size=800)
    nx.draw_networkx_nodes(graph, pos, nodelist=y_nodes, node_size=800)

    widths = [1 + 0.5 * graph[u][v]["weight"] for u, v in graph.edges()]
    nx.draw_networkx_edges(graph, pos, width=widths, alpha=0.7)

    labels = {n: f"{n[0]}:{n[1]}" for n in graph.nodes()}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=10)

    edge_labels = {(u, v): graph[u][v]["weight"] for u, v in graph.edges()}
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_size=9)

    plt.axis("off")
    plt.tight_layout()
    plt.show()
