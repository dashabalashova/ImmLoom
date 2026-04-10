from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATASET = "dataset-02"
SRC = PROJECT_ROOT / "data" / "raw" / "dataset-01" / "igdetective_outputs"
DST = PROJECT_ROOT / "data" / "raw" / DATASET / "igdetective_outputs"

SAMPLES = [
    "AG05253_1_hap1",
    "AG05253_1_hap2",
    "AG06939_13_hap1",
    "AG06939_13_hap2",
    "PR00251_hap1",
    "PR00251_hap2",
    "PR00366_hap1",
    "PR00366_hap2",
]


def ensure_dataset_01() -> None:
    if SRC.exists():
        print("dataset-01 already exists.")
        return

    print("dataset-01 not found. Creating it...")

    script_path = PROJECT_ROOT / "preprocessing" / "create_dataset-01.py"

    if not script_path.exists():
        raise FileNotFoundError(
            f"Cannot find script to create dataset-01: {script_path}"
        )

    subprocess.run(
        [sys.executable, str(script_path)],
        check=True,
    )

    if not SRC.exists():
        raise RuntimeError("dataset-01 creation failed.")


def copy_sample(sample: str) -> None:
    src_dir = SRC / f"{sample}_igdetective"
    dst_dir = DST / f"{sample}_igdetective"

    if not src_dir.exists():
        raise FileNotFoundError(f"Source sample directory not found: {src_dir}")

    if dst_dir.exists():
        print(f"Skipping existing sample: {sample}")
        return

    print(f"Copying {sample}...")
    shutil.copytree(src_dir, dst_dir)


def main() -> None:
    ensure_dataset_01()

    DST.mkdir(parents=True, exist_ok=True)

    for sample in SAMPLES:
        copy_sample(sample)

    print(f"{DATASET} is ready.")


if __name__ == "__main__":
    main()