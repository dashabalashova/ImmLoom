#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def reverse_complement_fasta(input_file: Path, output_file: Path) -> None:
    records = []
    for record in SeqIO.parse(str(input_file), "fasta"):
        record.seq = record.seq.reverse_complement()
        records.append(record)
    SeqIO.write(records, str(output_file), "fasta")


def reverse_coordinates(df: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    """
    Reverse coordinates within one sequence axis using the maximum end coordinate
    observed in the table, same as in old.py.
    """
    max_end = df[end_col].max()
    df[f"{start_col}_new"] = max_end + 1 - df[start_col]
    df[f"{end_col}_new"] = max_end + 1 - df[end_col]
    df[start_col] = df[f"{start_col}_new"]
    df[end_col] = df[f"{end_col}_new"]
    df = df.drop(columns=[f"{start_col}_new", f"{end_col}_new"])
    return df


def flip_strand(df: pd.DataFrame, strand_col: str = "strand2") -> pd.DataFrame:
    """
    Flip strand signs in pairwise alignments, same idea as in old.py.
    Unknown values are preserved.
    """
    if strand_col in df.columns:
        df[strand_col] = df[strand_col].map({"+": "-", "-": "+"}).fillna(df[strand_col])
    return df


def get_orientation_flags(pairwise_dir: Path) -> dict[int, bool]:
    """
    Determine sample orientation from pair_0* files.

    Convention:
    - True  -> forward (+)
    - False -> reverse (-)

    Sample 0 is assumed to be forward, same as in old.py.
    """
    pair0_files = sorted(
        p for p in pairwise_dir.iterdir()
        if p.is_file() and p.name.startswith("pair_0")
    )

    orientation = {0: True}

    for n, filepath in enumerate(pair0_files):
        df = pd.read_csv(filepath, sep="\t")

        plus_sum = df.loc[df["strand2"] == "+", "length1"].sum()
        minus_sum = df.loc[df["strand2"] == "-", "length1"].sum()

        orientation[n + 1] = bool(plus_sum > minus_sum)

    print(f"Orientation flags: {orientation}  # True = forward (+), False = reverse (-)")
    return orientation


def parse_pair_indices(filename: str) -> tuple[int, int]:
    """
    Example:
      pair_0-s01-h1_1-s01-h2.tsv -> (0, 1)
    """
    stem = Path(filename).stem
    parts = stem.split("_")
    if len(parts) < 3:
        raise ValueError(f"Unexpected pair filename format: {filename}")

    x_idx = int(parts[1].split("-")[0])
    y_idx = int(parts[2].split("-")[0])
    return x_idx, y_idx


def parse_self_index_and_sample(filename: str) -> tuple[int, str]:
    """
    Example:
      self_0-s01-h1.tsv -> (0, 's01-h1')
    """
    stem = Path(filename).stem
    parts = stem.split("_")
    if len(parts) != 2:
        raise ValueError(f"Unexpected self filename format: {filename}")

    idx_and_sample = parts[1]          # 0-s01-h1
    idx = int(idx_and_sample.split("-")[0])
    sample_name = idx_and_sample.split("-", 1)[1]
    return idx, sample_name


def process_pair_files(input_dir: Path, output_dir: Path, orientation: dict[int, bool]) -> None:
    pair_files = sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.name.startswith("pair")
    )

    for filepath in pair_files:
        df = pd.read_csv(filepath, sep="\t")
        x_idx, y_idx = parse_pair_indices(filepath.name)

        x_forward = orientation[x_idx]
        y_forward = orientation[y_idx]

        if x_forward is False:
            df = flip_strand(df)
            df = reverse_coordinates(df, "start1", "end1")
            df["start1"], df["end1"] = df["end1"], df["start1"]

        if y_forward is False:
            df = flip_strand(df)
            df = reverse_coordinates(df, "start2+", "end2+")
            df["start2+"], df["end2+"] = df["end2+"], df["start2+"]

        outpath = output_dir / filepath.name
        df.to_csv(outpath, sep="\t", index=False)


def process_self_files(
    input_dir: Path,
    output_dir: Path,
    fasta_input_dir: Path,
    fasta_output_dir: Path,
    locus: str,
    orientation: dict[int, bool],
) -> None:
    self_files = sorted(
        p for p in input_dir.iterdir()
        if p.is_file() and p.name.startswith("self")
    )

    for filepath in self_files:
        df = pd.read_csv(filepath, sep="\t")

        x_idx, sample_name = parse_self_index_and_sample(filepath.name)
        x_forward = orientation[x_idx]

        input_fasta = fasta_input_dir / f"{sample_name}-{locus}.fasta"
        output_fasta = fasta_output_dir / f"{sample_name}-{locus}.fasta"

        if not input_fasta.exists():
            raise FileNotFoundError(
                f"Expected FASTA file not found: {input_fasta}\n"
                f"Expected naming pattern: <sample>-<locus>.fasta, e.g. s01-h1-IGH.fasta"
            )

        if x_forward is False:
            df = reverse_coordinates(df, "start1", "end1")
            df["start1"], df["end1"] = df["end1"], df["start1"]

            df = reverse_coordinates(df, "start2+", "end2+")
            df["start2+"], df["end2+"] = df["end2+"], df["start2+"]

            reverse_complement_fasta(input_fasta, output_fasta)
        else:
            shutil.copy(input_fasta, output_fasta)

        outpath = output_dir / filepath.name
        df.to_csv(outpath, sep="\t", index=False)


def validate_inputs(pairwise_input_dir: Path, fasta_input_dir: Path) -> None:
    if not pairwise_input_dir.exists():
        raise FileNotFoundError(f"Pairwise input directory does not exist: {pairwise_input_dir}")
    if not fasta_input_dir.exists():
        raise FileNotFoundError(f"FASTA input directory does not exist: {fasta_input_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Orient FASTA and PWP outputs to forward strand."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset name, e.g. dataset-01",
    )
    parser.add_argument(
        "--locus",
        required=True,
        help="Locus name, e.g. IGH",
    )
    parser.add_argument(
        "--base-dir",
        default="data/processed",
        help="Base directory for processed data (default: data/processed)",
    )

    args = parser.parse_args()

    dataset = args.dataset
    locus = args.locus
    base_dir = Path(args.base_dir)

    pairwise_input_dir = base_dir / dataset / locus / "pwp" / "pairwise_alignments"
    pairwise_output_dir = base_dir / dataset / locus / "pwp" / "pairwise_alignments_forward"

    fasta_input_dir = base_dir / dataset / "fasta" / "mixed"
    fasta_output_dir = base_dir / dataset / "fasta" / "mixed_forward"

    validate_inputs(pairwise_input_dir, fasta_input_dir)

    pairwise_output_dir.mkdir(parents=True, exist_ok=True)
    fasta_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing dataset={dataset}, locus={locus}")
    print(f"Pairwise input:  {pairwise_input_dir}")
    print(f"Pairwise output: {pairwise_output_dir}")
    print(f"FASTA input:     {fasta_input_dir}")
    print(f"FASTA output:    {fasta_output_dir}")

    orientation = get_orientation_flags(pairwise_input_dir)

    process_pair_files(
        input_dir=pairwise_input_dir,
        output_dir=pairwise_output_dir,
        orientation=orientation,
    )

    process_self_files(
        input_dir=pairwise_input_dir,
        output_dir=pairwise_output_dir,
        fasta_input_dir=fasta_input_dir,
        fasta_output_dir=fasta_output_dir,
        locus=locus,
        orientation=orientation,
    )


if __name__ == "__main__":
    main()