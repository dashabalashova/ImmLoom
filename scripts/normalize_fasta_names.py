#!/usr/bin/env python3
"""
Normalize FASTA filenames and copy to processed folder.
Rules:
- sample prefix (before first underscore) -> assigned sequential IDs s01, s02, ... (sorted by prefix name)
- haplotypes: 'mat' or 'hap1' -> h1 ; 'pat' or 'hap2' -> h2
- IG token = last underscore-separated token before extension (e.g. IGK, IGH, TRA, TRB, TRG, IGL)
- new filename: {sample_id}_{hap}_{IG}.fasta  (e.g. s01_h1_IGK.fasta)
- sample_name_was: original filename without hap token and without the IG token
- sample_name_became: sample_id (e.g. s01)
Outputs:
- copied FASTA files in dst directory
- TSV mapping with columns: original_fullname, new_fullname, sample_name_was, sample_name_became, src_path, dst_path
"""
from pathlib import Path
import shutil
import csv
import argparse
import sys

HAP_MAP = {"mat": "h1", "hap1": "h1", "pat": "h2", "hap2": "h2"}

def parse_args():
    p = argparse.ArgumentParser(description="Normalize fasta filenames and copy to processed directory.")
    p.add_argument("--src", default="data/raw/dataset_01/fasta", help="Source directory with .fasta files")
    p.add_argument("--dst", default="data/processed/dataset_01/fasta", help="Destination directory for processed files")
    p.add_argument("--tsv", default="data/processed/dataset_01/rename_map.tsv", help="Output TSV mapping file")
    # p.add_argument("--force", action="store_true", help="Overwrite existing processed files")
    return p.parse_args()

def collect_fasta_files(src_dir: Path):
    if not src_dir.exists() or not src_dir.is_dir():
        raise FileNotFoundError(f"Source directory not found: {src_dir}")
    files = sorted([p for p in src_dir.iterdir() if p.is_file() and p.suffix.lower() == ".fasta"], key=lambda p: p.name)
    return files

def build_sample_map(files):
    prefixes = [p.name.split("_", 1)[0] for p in files]
    unique = sorted(sorted(set(prefixes), key=str.lower))
    return {pref: f"s{(i+1):02d}" for i, pref in enumerate(unique)}

def find_hap_and_ig(tokens):
    # tokens: list of parts of filename (without .fasta)
    hap_token = None
    hap_code = None
    hap_index = None
    for i, t in enumerate(tokens):
        if t in HAP_MAP:
            hap_token = t
            hap_code = HAP_MAP[t]
            hap_index = i
            break
    # fallback: substrings "hap1"/"hap2"
    if hap_token is None:
        for i, t in enumerate(tokens):
            if "hap1" in t:
                hap_token = t; hap_code = "h1"; hap_index = i; break
            if "hap2" in t:
                hap_token = t; hap_code = "h2"; hap_index = i; break
    if hap_code is None:
        hap_code = "h?"  # unknown
    ig_token = tokens[-1] if len(tokens) >= 2 else ""
    return hap_token, hap_code, hap_index, ig_token

def sample_name_without_hap_and_ig(tokens, hap_index):
    parts = []
    for i, t in enumerate(tokens):
        if i == hap_index:
            continue
        if i == len(tokens) - 1:
            continue
        parts.append(t)
    return "_".join(parts)

def main():
    args = parse_args()
    src = Path(args.src)
    dst = Path(args.dst)
    tsv_out = Path(args.tsv)
    force = True

    try:
        fasta_files = collect_fasta_files(src)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        sys.exit(2)

    if not fasta_files:
        print(f"No .fasta files found in {src}", file=sys.stderr)
        sys.exit(1)

    dst.mkdir(parents=True, exist_ok=True)
    tsv_out.parent.mkdir(parents=True, exist_ok=True)

    sample_map = build_sample_map(fasta_files)

    rows = []
    collisions = {}
    for p in fasta_files:
        name = p.name
        stem = p.stem
        tokens = stem.split("_")
        prefix = tokens[0]
        sample_id = sample_map[prefix]
        hap_token, hap_code, hap_index, ig_token = find_hap_and_ig(tokens)
        sample_was = sample_name_without_hap_and_ig(tokens, hap_index)
        sample_became = sample_id
        new_name = f"{sample_id}_{hap_code}_{ig_token}.fasta"
        target = dst / new_name
        # if target.exists() and not force:
        #     collisions.setdefault(new_name, 0)
        #     collisions[new_name] += 1
        #     new_name = f"{sample_id}_{hap_code}_{ig_token}_dup{collisions[new_name]}.fasta"
        #     target = dst / new_name
        shutil.copy2(p, target)
        rows.append({
            "original_fullname": name,
            "new_fullname": new_name,
            "sample_name_was": sample_was,
            "sample_name_became": sample_became,
            # "src_path": str(p.resolve()),
            # "dst_path": str(target.resolve())
        })

    # Write TSV
    with tsv_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["original_fullname","new_fullname","sample_name_was","sample_name_became"], delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Processed {len(rows)} files.")
    print(f"Processed fasta copies saved to: {dst}")
    print(f"TSV mapping saved to: {tsv_out}")
    print("\nSample prefix -> sample_id mapping:")
    for k, v in sorted(sample_map.items(), key=lambda kv: kv[1]):
        print(f"  {k} -> {v}")

if __name__ == "__main__":
    main()
