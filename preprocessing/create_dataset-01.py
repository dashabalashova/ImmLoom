from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = PROJECT_ROOT / "data" / "raw" / "dataset-01"
ZIP_PATH = DATASET_DIR / "primate_igdetective_outputs.zip"
DOWNLOAD_URL = "https://drive.google.com/uc?id=1hPsaegHbII8sGf6toy8Hr0NBM80OtDZw"


def run_command(cmd: list[str]) -> None:
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_gdown() -> None:
    if shutil.which("gdown") is not None:
        print("gdown is already available.")
        return

    print("gdown not found. Installing with uv pip...")
    run_command([sys.executable, "-m", "uv", "pip", "install", "gdown"])


def download_zip() -> None:
    DATASET_DIR.mkdir(parents=True, exist_ok=True)

    if ZIP_PATH.exists():
        print(f"Zip already exists, skipping download: {ZIP_PATH}")
        return

    run_command(
        [
            "gdown",
            "-O",
            str(ZIP_PATH),
            DOWNLOAD_URL,
        ]
    )


def unpack_zip() -> None:
    extract_dir = DATASET_DIR / "igdetective_outputs"
    if extract_dir.exists():
        print(f"Output directory already exists, skipping unzip: {extract_dir}")
        return

    print(f"Unpacking {ZIP_PATH} to {DATASET_DIR}")
    shutil.unpack_archive(str(ZIP_PATH), str(DATASET_DIR))


def main() -> None:
    ensure_gdown()
    download_zip()
    unpack_zip()
    print("dataset-01 is ready.")


if __name__ == "__main__":
    main()