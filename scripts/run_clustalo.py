#!/usr/bin/env python3
import subprocess
from pathlib import Path

IN_DIR = Path("outputs/dataset_01/IGH/components/sequences")
OUT_DIR = Path("outputs/dataset_01/IGH/components/alignments")
OUT_DIR.mkdir(parents=True, exist_ok=True)

for fasta in sorted(IN_DIR.glob("*.fasta"))[18:]:
    base = fasta.stem  # comp_segments_comp-XX
    outname = base.replace("comp_segments_comp-", "comp-") + ".aln"
    outpath = OUT_DIR / outname
    print(f"Running: {fasta} -> {outpath}")
    cmd = ["clustalo", "-i", str(fasta), "-o", str(outpath), "--outfmt=clu", "--force"]
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f"ERROR: clustalo failed on {fasta} (rc={r.returncode})")
