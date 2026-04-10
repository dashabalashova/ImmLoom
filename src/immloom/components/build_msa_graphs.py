#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build pairwise-identity and coverage plots for component alignments."
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. dataset-01")
    parser.add_argument("--locus", required=True, help="Locus name, e.g. IGH")
    parser.add_argument(
        "--data-root",
        default="data/processed",
        help="Root directory containing processed datasets",
    )
    return parser.parse_args()


def non_gap_segments(seq: str) -> list[tuple[int, int]]:
    """Return list of (start, length) segments where seq is not a gap."""
    segments: list[tuple[int, int]] = []
    in_segment = False
    start = 0

    for i, ch in enumerate(seq):
        if ch != "-" and not in_segment:
            start = i
            in_segment = True
        elif ch == "-" and in_segment:
            segments.append((start, i - start))
            in_segment = False

    if in_segment:
        segments.append((start, len(seq) - start))

    return segments


def compute_pairwise_identity(alignment) -> pd.DataFrame:
    """
    Compute pairwise identity matrix with gaps included.

    Identity definition:
    - match if symbols are equal, including gap-gap
    - denominator is full alignment length
    """
    records = list(alignment)
    names = [rec.id for rec in records]
    seqs = [np.array(list(str(rec.seq))) for rec in records]

    n_seq = len(seqs)
    aln_len = alignment.get_alignment_length()

    identity = np.zeros((n_seq, n_seq), dtype=float)

    for i in range(n_seq):
        for j in range(n_seq):
            s1 = seqs[i]
            s2 = seqs[j]
            matches = (s1 == s2).sum()
            identity[i, j] = 100.0 * matches / aln_len if aln_len > 0 else 0.0

    return pd.DataFrame(identity, index=names, columns=names)


def plot_pairwise_identity(identity_df: pd.DataFrame, out_path: Path, component_id: str) -> None:
    n_seq = len(identity_df)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(identity_df.values, aspect="auto", vmin=0, vmax=100)

    ax.set_xticks(range(n_seq))
    ax.set_yticks(range(n_seq))
    ax.set_xticklabels(identity_df.columns, rotation=90, fontsize=8)
    ax.set_yticklabels(identity_df.index, fontsize=8)
    ax.set_title(f"Component {component_id}: pairwise identity (%) with gaps included")

    fig.colorbar(im, ax=ax, label="% identity")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_coverage(alignment, out_path: Path, component_id: str) -> None:
    records = list(alignment)
    n_seq = len(records)
    aln_len = alignment.get_alignment_length()

    fig_h = max(4, 0.4 * n_seq)
    fig, ax = plt.subplots(figsize=(16, fig_h))

    y_positions = []
    labels = []

    for idx, record in enumerate(records):
        y = n_seq - 1 - idx
        seq = str(record.seq)
        segments = non_gap_segments(seq)

        if segments:
            ax.broken_barh(segments, (y - 0.35, 0.7))

        y_positions.append(y)
        labels.append(record.id)

    ax.set_xlim(0, aln_len)
    ax.set_ylim(-1, n_seq)
    ax.set_xlabel("Alignment position")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_title(f"Component {component_id}: MSA coverage by sequence")
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    alignments_dir = (
        Path(args.data_root)
        / args.dataset
        / args.locus
        / "components_msa"
    )
    graphs_dir = (
        Path(args.data_root)
        / args.dataset
        / args.locus
        / "components_msa"
        / "graphs"
    )
    tables_dir = (
        Path(args.data_root)
        / args.dataset
        / args.locus
        / "components_msa"
        / "tables"
    )

    graphs_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    aln_files = sorted(alignments_dir.glob("component_*.aln"))
    if not aln_files:
        raise FileNotFoundError(f"No alignment files found in: {alignments_dir}")

    for aln_file in aln_files:
        component_id = aln_file.stem.replace("component_", "")
        print(f"Processing component {component_id}: {aln_file}")

        alignment = AlignIO.read(aln_file, "clustal")
        identity_df = compute_pairwise_identity(alignment)

        identity_csv = tables_dir / f"component_{component_id}_pairwise_identity.csv"
        identity_png = graphs_dir / f"component_{component_id}_pairwise_identity.png"
        coverage_png = graphs_dir / f"component_{component_id}_coverage.png"

        identity_df.round(2).to_csv(identity_csv)
        plot_pairwise_identity(identity_df, identity_png, component_id)
        plot_coverage(alignment, coverage_png, component_id)

        print(f"  saved table:  {identity_csv}")
        print(f"  saved figure: {identity_png}")
        print(f"  saved figure: {coverage_png}")


if __name__ == "__main__":
    main()