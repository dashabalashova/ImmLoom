#!/usr/bin/env python3
"""
Create per-IG CSV configs from processed FASTA files.

For each fasta file in src_dir, this script determines the IG token as the
last underscore-separated token in the filename (before .fasta), e.g.:
  s01_h1_IGK.fasta -> IG token = IGK

Then it writes files:
  <out_dir>/config_IGK.csv
with header: SampleID,Label,Fasta
SampleID = filename without extension (e.g. s01_h1_IGK)
Label = same as SampleID
Fasta = absolute path to the .fasta file
"""
from pathlib import Path
import argparse
import csv
import sys

def parse_args():
    p = argparse.ArgumentParser(description="Create per-IG config CSVs from processed fasta files.")
    p.add_argument("--src", default="data/processed/dataset_01/fasta", help="Directory with processed .fasta files")
    p.add_argument("--out", default="data/processed/dataset_01/configs", help="Output directory for config CSVs")
    p.add_argument("--force", action="store_true", help="Overwrite existing CSV files")
    return p.parse_args()

def main():
    args = parse_args()
    src_dir = Path(args.src)
    out_dir = Path(args.out)
    force = args.force

    if not src_dir.exists() or not src_dir.is_dir():
        print(f"ERROR: source directory does not exist: {src_dir}", file=sys.stderr)
        sys.exit(2)

    fasta_files = sorted([p for p in src_dir.iterdir() if p.is_file() and p.suffix.lower() == ".fasta"], key=lambda p: p.name)
    if not fasta_files:
        print(f"No .fasta files found in {src_dir}", file=sys.stderr)
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    # Group by IG token (last underscore-separated token)
    groups = {}
    for f in fasta_files:
        stem = f.stem  # filename without suffix
        parts = stem.split("_")
        if len(parts) < 2:
            ig = "UNKNOWN"
        else:
            ig = parts[-1].upper()
        groups.setdefault(ig, []).append(f)

    # For each IG group write CSV
    for ig, files in sorted(groups.items(), key=lambda kv: kv[0]):
        csv_name = f"config_{ig}.csv"
        csv_path = out_dir / csv_name
        if csv_path.exists() and not force:
            print(f"Skipping existing {csv_path} (use --force to overwrite)")
            continue

        # Prepare rows: SampleID, Label, Fasta (absolute path)
        rows = []
        for f in sorted(files, key=lambda p: p.name):
            sample_id = f.stem
            label = sample_id
            fasta_abs = str(f.resolve())
            rows.append((sample_id, label, fasta_abs))

        # Write CSV (comma-separated)
        with csv_path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["SampleID", "Label", "Fasta"])
            for r in rows:
                writer.writerow(r)

        print(f"Wrote {csv_path} ({len(rows)} entries)")

    print("Done.")

if __name__ == "__main__":
    main()
