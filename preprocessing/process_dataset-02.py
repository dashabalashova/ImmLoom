#!/usr/bin/env python3

import re
import shutil
from pathlib import Path

import pandas as pd


RAW_DATASET_DIR = Path("data/raw/dataset-02/igdetective_outputs")
OUT_DIR = Path("data/processed/dataset-02")
OUT_FASTA_DIR = OUT_DIR / "fasta" / "mixed"

METADATA_DIR = Path("data/metadata")
MAPPING_PATH = METADATA_DIR / "dataset-02_mapping.tsv"


def normalize_sample_name(folder_name: str) -> str:
    return re.sub(r"_igdetective$", "", folder_name)


def split_sample_name(sample_original: str):
    """
    Examples:
      AG05253_1_hap1 -> ("AG05253_1", "1")
      PR00251_hap2   -> ("PR00251", "2")
    """
    m = re.match(r"^(.*)_hap(\d+)$", sample_original)
    if not m:
        raise ValueError(f"Unexpected sample name format: {sample_original}")
    base_sample = m.group(1)
    hap = m.group(2)
    return base_sample, hap


def build_sample_mapping(sample_original_list):
    """
    Maps:
      AG05253_1_hap1   -> s01-h1
      AG05253_1_hap2   -> s01-h2
      AG06939_13_hap1  -> s02-h1
      AG06939_13_hap2  -> s02-h2
      ...
    """
    base_samples = sorted({split_sample_name(s)[0] for s in sample_original_list})
    base_to_std = {
        base_sample: f"s{i:02d}"
        for i, base_sample in enumerate(base_samples, start=1)
    }

    sample_mapping = {}
    for sample_original in sample_original_list:
        base_sample, hap = split_sample_name(sample_original)
        sample_mapping[sample_original] = f"{base_to_std[base_sample]}-h{hap}"

    return sample_mapping


def find_matching_fasta(igloci_fasta_dir: Path, locus: str, max_contig: str, max_numv: int) -> Path:
    """
    Expected file pattern:
      {Locus}_{Contig}_{NumV}Vs.fasta

    Example:
      IGH_h1tg000081l_89Vs.fasta
    """
    fname = f"{locus}_{max_contig}_{int(max_numv)}Vs.fasta"
    path = igloci_fasta_dir / fname

    if not path.exists():
        raise FileNotFoundError(f"Fasta file not found: {path}")

    return path


def save_mapping_tsv(sample_mapping: dict):
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    mapping_rows = []
    for sample_original, sample_std in sample_mapping.items():
        base_sample, hap = split_sample_name(sample_original)
        mapping_rows.append(
            {
                "sample_original": sample_original,
                "sample_std": sample_std,
                "base_sample": base_sample,
                "hap": hap,
            }
        )

    mapping_df = (
        pd.DataFrame(mapping_rows)
        .sort_values(["sample_std"])
        .reset_index(drop=True)
    )

    mapping_df.to_csv(MAPPING_PATH, sep="\t", index=False)
    print(f"Saved mapping: {MAPPING_PATH}")
    print(mapping_df)


def main():
    if not RAW_DATASET_DIR.exists():
        raise FileNotFoundError(f"Input directory not found: {RAW_DATASET_DIR}")

    OUT_FASTA_DIR.mkdir(parents=True, exist_ok=True)

    folders = sorted(
        x for x in RAW_DATASET_DIR.iterdir()
        if x.is_dir() and x.name.endswith("_igdetective")
    )

    if not folders:
        raise ValueError(f"No *_igdetective folders found in {RAW_DATASET_DIR}")

    sample_original_list = [normalize_sample_name(folder.name) for folder in folders]
    sample_mapping = build_sample_mapping(sample_original_list)

    save_mapping_tsv(sample_mapping)

    copied = []

    for folder in folders:
        sample_original = normalize_sample_name(folder.name)
        sample_std = sample_mapping[sample_original]

        summary_csv = folder / "refined_ig_loci" / "summary.csv"
        igloci_fasta_dir = folder / "refined_ig_loci" / "igloci_fasta"

        if not summary_csv.exists():
            print(f"WARNING: summary.csv not found, skipping: {summary_csv}")
            continue

        if not igloci_fasta_dir.exists():
            print(f"WARNING: igloci_fasta dir not found, skipping: {igloci_fasta_dir}")
            continue

        df = pd.read_csv(summary_csv)

        required_cols = {"Locus", "Contig", "NumV"}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"{summary_csv} is missing columns: {missing}")

        # Choose the contig with maximal NumV for each locus
        idx = df.groupby("Locus")["NumV"].idxmax()
        df_best = (
            df.loc[idx, ["Locus", "Contig", "NumV"]]
            .rename(columns={"Contig": "MaxContig", "NumV": "MaxNumV"})
            .sort_values(["Locus", "MaxContig"])
            .reset_index(drop=True)
        )

        for _, row in df_best.iterrows():
            locus = str(row["Locus"])
            max_contig = str(row["MaxContig"])
            max_numv = int(row["MaxNumV"])

            src_fasta = find_matching_fasta(igloci_fasta_dir, locus, max_contig, max_numv)
            dst_fasta = OUT_FASTA_DIR / f"{sample_std}-{locus}.fasta"

            shutil.copy2(src_fasta, dst_fasta)
            copied.append(
                {
                    "sample_original": sample_original,
                    "sample_std": sample_std,
                    "locus": locus,
                    "src_fasta": str(src_fasta),
                    "dst_fasta": str(dst_fasta),
                }
            )

            print(f"Copied: {src_fasta} -> {dst_fasta}")

    print("\nDone.")
    print(f"Saved {len(copied)} fasta files to: {OUT_FASTA_DIR}")

    if copied:
        copied_df = pd.DataFrame(copied).sort_values(["sample_std", "locus"]).reset_index(drop=True)
        print("\nCopied files summary:")
        print(copied_df[["sample_std", "locus", "dst_fasta"]])


if __name__ == "__main__":
    main()