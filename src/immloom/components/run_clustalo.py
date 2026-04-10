#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Clustal Omega on component FASTA files."
    )
    parser.add_argument("--dataset", required=True, help="Dataset name, e.g. dataset-01")
    parser.add_argument("--locus", required=True, help="Locus name, e.g. IGH")
    parser.add_argument(
        "--data-root",
        default="data/processed",
        help="Root directory containing processed datasets",
    )
    parser.add_argument(
        "--outfmt",
        default="clu",
        choices=["clu", "fa", "fasta", "msf", "phy", "selex", "st"],
        help="Clustal Omega output format",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if shutil.which("clustalo") is None:
        raise RuntimeError(
            "clustalo not found in PATH. Install it first, for example:\n"
            "brew install clustal-omega"
        )

    input_dir = (
        Path(args.data_root)
        / args.dataset
        / args.locus
        / "components"
        / "fasta"
    )
    output_dir = (
        Path(args.data_root)
        / args.dataset
        / args.locus
        / "components"
        / "msa"
    )

    if not input_dir.exists():
        raise FileNotFoundError(f"Input FASTA directory not found: {input_dir}")

    fasta_files = sorted(input_dir.glob("component_*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(f"No component FASTA files found in: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    for fasta_path in fasta_files:
        base = fasta_path.stem
        out_path = output_dir / f"{base}.aln"

        cmd = [
            "clustalo",
            "-i", str(fasta_path),
            "-o", str(out_path),
            "--outfmt", args.outfmt,
            "--force",
        ]

        subprocess.run(cmd, check=True)
        print(f"Done {base}")

    print(f"Saved alignments to: {output_dir}")


if __name__ == "__main__":
    main()