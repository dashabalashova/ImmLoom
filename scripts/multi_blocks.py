#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import networkx as nx
import pandas as pd
from networkx.algorithms import community


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a multi-sample block graph from pairwise TSV tables, save connected "
            "components, refine oversized components by community detection, and export graph visualizations."
        )
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. dataset_01_subset")
    parser.add_argument("--locus", required=True, help="Locus name, e.g. IGH")
    parser.add_argument(
        "--results-root",
        default="results",
        help="Root directory for results (default: results)",
    )
    parser.add_argument(
        "--use-multigraph",
        action="store_true",
        help="Use nx.MultiGraph instead of nx.Graph",
    )
    parser.add_argument(
        "--min-component-size",
        type=int,
        default=2,
        help="Only draw components with at least this many nodes (default: 2)",
    )
    parser.add_argument(
        "--layout-seed",
        type=int,
        default=42,
        help="Seed for spring layout used in component plots (default: 42)",
    )
    parser.add_argument(
        "--max-nodes-per-group",
        type=float,
        default=3.0,
        help=(
            "If len(component) / number_of_unique_groups is greater than this threshold, "
            "the component is refined with greedy modularity communities (default: 3.0)"
        ),
    )
    return parser.parse_args()


def build_paths(dataset: str, locus: str, results_root: str) -> dict[str, Path]:
    root = Path(results_root)
    pair_dir = root / dataset / "pair" / locus / "tables"
    out_tables = root / dataset / "multi" / locus / "tables"
    out_graphs = root / dataset / "multi" / locus / "graphs"
    return {
        "pair_dir": pair_dir,
        "out_tables": out_tables,
        "out_graphs": out_graphs,
    }


def build_block_graph(pair_dir: Path, use_multi: bool = False) -> tuple[nx.Graph, list[Path]]:
    files = sorted(pair_dir.glob("*.tsv"))
    graph: nx.Graph = nx.MultiGraph() if use_multi else nx.Graph()

    for file_path in files:
        try:
            group_x, group_y = file_path.stem.split("-", maxsplit=1)
        except ValueError as exc:
            raise ValueError(
                f"Unexpected filename format: {file_path.name}. Expected 'groupX-groupY.tsv'."
            ) from exc

        df = pd.read_csv(file_path, sep="\t")
        required_columns = {"x_block_id", "y_block_id"}
        missing = required_columns - set(df.columns)
        if missing:
            raise KeyError(f"{file_path.name}: missing required columns: {sorted(missing)}")

        df = df.dropna(subset=["x_block_id", "y_block_id"]).copy()
        if df.empty:
            continue

        df["x_block_id"] = df["x_block_id"].astype(int)
        df["y_block_id"] = df["y_block_id"].astype(int)

        for _, row in df.iterrows():
            node_x = (group_x, int(row["x_block_id"]))
            node_y = (group_y, int(row["y_block_id"]))

            graph.add_node(node_x, group=group_x, block_id=int(row["x_block_id"]))
            graph.add_node(node_y, group=group_y, block_id=int(row["y_block_id"]))
            graph.add_edge(node_x, node_y, source_file=file_path.name)

    return graph, files


def should_refine_component(
    graph: nx.Graph,
    component: set[tuple[str, int]],
    max_nodes_per_group: float,
) -> bool:
    unique_groups = {graph.nodes[node]["group"] for node in component}
    if not unique_groups:
        return False
    ratio = len(component) / len(unique_groups)
    return ratio > max_nodes_per_group


def refine_components(
    graph: nx.Graph,
    max_nodes_per_group: float,
) -> list[set[tuple[str, int]]]:
    """
    Start from connected components.
    If a component has too many nodes per unique group, split it using
    greedy modularity community detection.
    """
    refined_components: list[set[tuple[str, int]]] = []

    for component in nx.connected_components(graph):
        component = set(component)

        if not should_refine_component(graph, component, max_nodes_per_group):
            refined_components.append(component)
            continue

        subgraph = graph.subgraph(component).copy()

        # Convert to simple Graph for community detection.
        simple_subgraph = nx.Graph()
        simple_subgraph.add_nodes_from(subgraph.nodes(data=True))
        simple_subgraph.add_edges_from(subgraph.edges())

        communities = list(community.greedy_modularity_communities(simple_subgraph))

        if len(communities) <= 1:
            refined_components.append(component)
            continue

        refined_components.extend(set(c) for c in communities)

    refined_components = sorted(refined_components, key=len, reverse=True)
    return refined_components


def components_df_from_list(
    graph: nx.Graph,
    components: list[set[tuple[str, int]]],
) -> pd.DataFrame:
    rows: list[dict[str, int | str]] = []

    for component_id, component in enumerate(components, start=1):
        for node in component:
            rows.append(
                {
                    "component_id": component_id,
                    "group": graph.nodes[node]["group"],
                    "block_id": graph.nodes[node]["block_id"],
                }
            )

    return pd.DataFrame(rows)


def make_global_layout(graph: nx.Graph) -> tuple[dict[tuple[str, int], tuple[float, float]], dict[str, int]]:
    group_order = sorted({attrs["group"] for _, attrs in graph.nodes(data=True)})
    group_to_x = {group: idx for idx, group in enumerate(group_order)}

    pos: dict[tuple[str, int], tuple[float, float]] = {}
    for group in group_order:
        nodes_in_group = [node for node, attrs in graph.nodes(data=True) if attrs["group"] == group]
        nodes_in_group = sorted(nodes_in_group, key=lambda node: node[1])

        for j, node in enumerate(nodes_in_group):
            pos[node] = (group_to_x[group], -j)

    return pos, group_to_x


def save_global_plot(graph: nx.Graph, out_path: Path) -> None:
    if graph.number_of_nodes() == 0:
        return

    pos, group_to_x = make_global_layout(graph)
    labels = {node: f"{node[1]}" for node in graph.nodes()}

    fig, ax = plt.subplots(figsize=(14, 10))
    nx.draw(
        graph,
        pos,
        labels=labels,
        with_labels=True,
        node_size=500,
        font_size=8,
        width=1.0,
        ax=ax,
    )

    for group, x in group_to_x.items():
        ax.text(x, 1, group, ha="center", va="bottom", fontsize=12)

    ax.set_title("Blocks grouped by sample/haplotype")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_group_colormap(groups: Iterable[str]) -> dict[str, tuple[float, float, float, float]]:
    groups = sorted(set(groups))
    cmap = mpl.colormaps.get_cmap("tab20")
    denominator = max(1, len(groups) - 1)
    return {group: cmap(i / denominator) for i, group in enumerate(groups)}


def save_component_plots(
    graph: nx.Graph,
    components: list[set[tuple[str, int]]],
    out_dir: Path,
    min_component_size: int = 2,
    layout_seed: int = 42,
) -> None:
    if graph.number_of_nodes() == 0 or not components:
        return

    components = [component for component in components if len(component) >= min_component_size]
    if not components:
        return

    all_groups = [graph.nodes[node]["group"] for component in components for node in component]
    group_to_color = build_group_colormap(all_groups)

    components = sorted(components, key=len, reverse=True)

    for i, component in enumerate(components, start=1):
        subgraph = graph.subgraph(component)
        pos = nx.spring_layout(subgraph, seed=layout_seed)

        node_colors = [
            group_to_color.get(graph.nodes[node]["group"], (0.8, 0.8, 0.8, 1.0))
            for node in subgraph.nodes
        ]
        labels = {
            node: f"{graph.nodes[node]['group']}\n{graph.nodes[node]['block_id']}"
            for node in subgraph.nodes
        }

        fig, ax = plt.subplots(figsize=(7, 7))
        nx.draw(
            subgraph,
            pos,
            node_color=node_colors,
            labels=labels,
            with_labels=True,
            node_size=500,
            font_size=8,
            ax=ax,
        )

        used_groups = sorted({graph.nodes[node]["group"] for node in subgraph.nodes})
        legend_elements = [Patch(facecolor=group_to_color[group], label=group) for group in used_groups]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_title(f"Component {i} (|V|={len(component)})")
        ax.axis("off")

        fig.tight_layout()
        fig.savefig(out_dir / f"component_{i:03d}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    args = parse_args()
    paths = build_paths(args.dataset, args.locus, args.results_root)

    pair_dir = paths["pair_dir"]
    out_tables = paths["out_tables"]
    out_graphs = paths["out_graphs"]

    if not pair_dir.exists():
        raise FileNotFoundError(f"Pair tables directory does not exist: {pair_dir}")

    out_tables.mkdir(parents=True, exist_ok=True)
    out_graphs.mkdir(parents=True, exist_ok=True)

    graph, files = build_block_graph(pair_dir, use_multi=args.use_multigraph)

    print(f"Pair directory: {pair_dir}")
    print(f"Found files: {len(files)}")
    print(f"Nodes: {graph.number_of_nodes()}")
    print(f"Edges: {graph.number_of_edges()}")

    if graph.number_of_nodes() == 0:
        print("Graph is empty. Nothing to save.")
        return

    original_components = list(nx.connected_components(graph))
    original_component_sizes = sorted((len(c) for c in original_components), reverse=True)

    refined_components = refine_components(
        graph,
        max_nodes_per_group=args.max_nodes_per_group,
    )
    refined_component_sizes = sorted((len(c) for c in refined_components), reverse=True)

    components_df = components_df_from_list(graph, refined_components)
    components_path = out_tables / "components.tsv"
    components_df.to_csv(components_path, index=False, sep="\t")

    save_global_plot(graph, out_graphs / "blocks_grouped_by_sample.png")
    save_component_plots(
        graph,
        refined_components,
        out_graphs,
        min_component_size=args.min_component_size,
        layout_seed=args.layout_seed,
    )

    print(f"Original connected components: {len(original_components)}")
    print(f"Largest original component sizes: {original_component_sizes[:10]}")
    print(f"Refined components: {len(refined_components)}")
    print(f"Largest refined component sizes: {refined_component_sizes[:10]}")
    print(f"Saved table: {components_path}")
    print(f"Saved graphs to: {out_graphs}")


if __name__ == "__main__":
    main()