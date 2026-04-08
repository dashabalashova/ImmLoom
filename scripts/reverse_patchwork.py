#!/usr/bin/env python3

import os
import shutil
import argparse
import pandas as pd
from Bio import SeqIO


def reverse_complement_fasta(input_file: str, output_file: str) -> None:
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        record.seq = record.seq.reverse_complement()
        records.append(record)
    SeqIO.write(records, output_file, "fasta")


def reverse_coordinates(df: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    max_end = df[end_col].max()
    df[start_col + "_new"] = max_end + 1 - df[start_col]
    df[end_col + "_new"] = max_end + 1 - df[end_col]
    df[start_col] = df[start_col + "_new"]
    df[end_col] = df[end_col + "_new"]
    df = df.drop(columns=[start_col + "_new", end_col + "_new"])
    return df


def flip_strand(df: pd.DataFrame, strand_col: str = "strand2") -> pd.DataFrame:
    df[strand_col] = df[strand_col].map({"+": "-", "-": "+"})
    return df


def get_orientation_flags(pairwise_dir: str) -> dict:
    q1 = sorted(os.listdir(pairwise_dir))
    q1 = [q for q in q1 if q.startswith("pair_0")]

    d1 = {0: True}

    for n, filename in enumerate(q1):
        filepath = os.path.join(pairwise_dir, filename)
        df = pd.read_csv(filepath, sep="\t")

        plus_sum = df.loc[df["strand2"] == "+", "length1"].sum()
        minus_sum = df.loc[df["strand2"] == "-", "length1"].sum()

        d1[n + 1] = bool(plus_sum > minus_sum)

    print(f"Orientation flags: {d1}  # True = forward (+), False = reverse (-)")
    return d1


def process_pair_files(input_dir: str, output_dir: str, d1: dict) -> None:
    q2 = sorted(os.listdir(input_dir))
    q2 = [q for q in q2 if q.startswith("pair")]

    for q in q2:
        filepath = os.path.join(input_dir, q)
        df = pd.read_csv(filepath, sep="\t")

        q3 = q.split("_")[1:4]
        x_idx = int(q3[0].split("-")[0])
        y_idx = int(q3[2].split("-")[0])

        x_true = d1[x_idx]
        y_true = d1[y_idx]

        if x_true is False:
            df = flip_strand(df)
            df = reverse_coordinates(df, "start1", "end1")
            df["start1"], df["end1"] = df["end1"], df["start1"]

        if y_true is False:
            df = flip_strand(df)
            df = reverse_coordinates(df, "start2+", "end2+")
            df["start2+"], df["end2+"] = df["end2+"], df["start2+"]

        outpath = os.path.join(output_dir, q)
        df.to_csv(outpath, sep="\t", index=False)


def process_self_files(
    input_dir: str,
    output_dir: str,
    fasta_input_dir: str,
    fasta_output_dir: str,
    locus: str,
    d1: dict,
) -> None:
    q2 = sorted(os.listdir(input_dir))
    q2 = [q for q in q2 if q.startswith("self")]

    for q in q2:
        filepath = os.path.join(input_dir, q)
        df = pd.read_csv(filepath, sep="\t")

        q3 = q.split("_")[1]
        qf = q.split("-")[1][:-4]
        x_idx = int(q3.split("-")[0])
        x_true = d1[x_idx]

        input_fasta = os.path.join(fasta_input_dir, f"{qf}_{locus}.fasta")
        output_fasta = os.path.join(fasta_output_dir, f"{qf}_{locus}.fasta")

        if x_true is False:
            df = reverse_coordinates(df, "start1", "end1")
            df["start1"], df["end1"] = df["end1"], df["start1"]
            df = reverse_coordinates(df, "start2+", "end2+")
            df["start2+"], df["end2+"] = df["end2+"], df["start2+"]

            reverse_complement_fasta(input_fasta, output_fasta)
        else:
            shutil.copy(input_fasta, output_fasta)

        outpath = os.path.join(output_dir, q)
        df.to_csv(outpath, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Reverse patchwork alignments and FASTA sequences based on strand orientation."
    )

    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset name (e.g., dataset_01_subset)",
    )

    parser.add_argument(
        "--locus",
        required=True,
        help="Locus name (e.g., IGH)",
    )

    parser.add_argument(
        "--base_dir",
        default="data/processed",
        help="Base directory for processed data (default: data/processed)",
    )

    args = parser.parse_args()

    dataset = args.dataset
    locus = args.locus
    base_dir = args.base_dir

    pairwise_input_dir = os.path.join(base_dir, dataset, "patchwork_output", locus, "pairwise_alignments")
    pairwise_output_dir = os.path.join(base_dir, dataset, "patchwork_output_reversed", locus, "pairwise_alignments")
    fasta_input_dir = os.path.join(base_dir, dataset, "fasta")
    fasta_output_dir = os.path.join(base_dir, dataset, "fasta_reversed")

    os.makedirs(pairwise_output_dir, exist_ok=True)
    os.makedirs(fasta_output_dir, exist_ok=True)

    print(f"Processing dataset: {dataset}, locus: {locus}")

    d1 = get_orientation_flags(pairwise_input_dir)

    process_pair_files(pairwise_input_dir, pairwise_output_dir, d1)
    process_self_files(
        pairwise_input_dir,
        pairwise_output_dir,
        fasta_input_dir,
        fasta_output_dir,
        locus,
        d1,
    )

    print("Done.")


if __name__ == "__main__":
    main()