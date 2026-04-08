#!/usr/bin/env python3

import os
import re
import shutil
import argparse
import pandas as pd
from pathlib import Path


def normalize_sample_name(folder_name: str) -> str:
    return re.sub(r"_igdetective$", "", folder_name)


def split_sample_name(sample_original: str):
    m = re.match(r"^(.*)_hap(\d+)$", sample_original)
    if not m:
        raise ValueError(f"Unexpected sample name format: {sample_original}")
    base_sample = m.group(1)
    hap = m.group(2)
    return base_sample, hap


def build_sample_mapping(sample_original_list):
    base_samples = sorted({split_sample_name(s)[0] for s in sample_original_list})
    base_to_std = {
        base_sample: f"s{i:02d}"
        for i, base_sample in enumerate(base_samples, start=1)
    }

    sample_mapping = {}
    for sample_original in sample_original_list:
        base_sample, hap = split_sample_name(sample_original)
        sample_mapping[sample_original] = f"{base_to_std[base_sample]}_h{hap}"

    return sample_mapping


def find_matching_fasta(igloci_fasta_dir: Path, locus: str, max_contig: str, max_numv: int) -> Path:
    fname = f"{locus}_{max_contig}_{int(max_numv)}Vs.fasta"
    path = igloci_fasta_dir / fname

    if not path.exists():
        raise FileNotFoundError(f"Fasta file not found: {path}")

    return path


def main():
    parser = argparse.ArgumentParser(description="Build locus-contig table from IgDetective outputs")
    parser.add_argument("--dataset", required=True, help="Dataset name (e.g. dataset_01_subset)")
    parser.add_argument("--input-dir", default=None, help="Optional override for input dir")
    parser.add_argument("--output-dir", default=None, help="Optional override for output dir")

    args = parser.parse_args()

    dataset = args.dataset

    input_dir = Path(args.input_dir) if args.input_dir else Path(f"data/raw/{dataset}/igdetective_outputs")
    output_dir = Path(args.output_dir) if args.output_dir else Path(f"data/processed/{dataset}")

    output_tsv = output_dir / "locus_contig_table.tsv"
    output_fasta_dir = output_dir / "fasta"
    output_config_dir = output_dir / "config"

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    folders = sorted(
        x for x in os.listdir(input_dir)
        if x.endswith("_igdetective") and (input_dir / x).is_dir()
    )

    sample_original_list = [normalize_sample_name(folder) for folder in folders]
    sample_mapping = build_sample_mapping(sample_original_list)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_fasta_dir.mkdir(parents=True, exist_ok=True)
    output_config_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    for folder in folders:
        sample_original = normalize_sample_name(folder)
        sample_std = sample_mapping[sample_original]

        summary_csv = input_dir / folder / "refined_ig_loci" / "summary.csv"
        igloci_fasta_dir = input_dir / folder / "refined_ig_loci" / "igloci_fasta"

        if not summary_csv.exists():
            print(f"WARNING: summary.csv not found, skipping: {summary_csv}")
            continue

        df = pd.read_csv(summary_csv)

        required_cols = {"Locus", "Contig", "NumV"}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"{summary_csv} is missing columns: {missing}")

        idx = df.groupby("Locus")["NumV"].idxmax()
        df_best = (
            df.loc[idx, ["Locus", "Contig", "NumV"]]
            .rename(columns={"Contig": "MaxContig", "NumV": "MaxNumV"})
            .sort_values(["Locus", "MaxContig"])
            .reset_index(drop=True)
        )

        for _, row in df_best.iterrows():
            locus = row["Locus"]
            max_contig = row["MaxContig"]
            max_numv = int(row["MaxNumV"])

            dst_fasta = (output_fasta_dir / f"{sample_std}_{locus}.fasta").resolve()

            src_fasta = find_matching_fasta(igloci_fasta_dir, locus, max_contig, max_numv)
            shutil.copy2(src_fasta, dst_fasta)

            print(f"Copied: {src_fasta} -> {dst_fasta}")

            all_rows.append({
                "sample_original": sample_original,
                "sample_std": sample_std,
                "Locus": locus,
                "MaxContig": max_contig,
                "MaxNumV": max_numv,
                "Fasta": str(dst_fasta),
            })

    result_df = pd.DataFrame(all_rows)
    result_df = result_df.sort_values(["sample_std", "Locus"]).reset_index(drop=True)
    result_df.to_csv(output_tsv, sep="\t", index=False)

    print(f"\nSaved TSV: {output_tsv}")
    print(result_df)

    # configs
    for locus in sorted(result_df["Locus"].unique()):
        config_df = (
            result_df[result_df["Locus"] == locus][["sample_std", "Fasta"]]
            .rename(columns={"sample_std": "SampleID"})
            .copy()
        )
        config_df["Label"] = config_df["SampleID"]
        config_df = config_df[["SampleID", "Label", "Fasta"]] \
            .sort_values("SampleID") \
            .reset_index(drop=True)

        config_path = output_config_dir / f"config_{locus}.csv"
        config_df.to_csv(config_path, index=False)

        print(f"Saved config: {config_path}")


if __name__ == "__main__":
    main()