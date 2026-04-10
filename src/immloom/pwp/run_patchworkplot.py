from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run patchworkplot for a given dataset and locus."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="Dataset name, e.g. dataset-02",
    )
    parser.add_argument(
        "--locus",
        required=True,
        help="Locus name, e.g. IGH",
    )
    parser.add_argument(
        "--config",
        default=None,
        help=(
            "Optional explicit path to patchworkplot config TSV. "
            "Default: data/processed/{dataset}/config/pwp_{locus}.tsv"
        ),
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Optional explicit output directory. "
            "Default: data/processed/{dataset}/{locus}/pwp"
        ),
    )
    parser.add_argument(
        "--cmap",
        default="PuBuGn",
        help="Colormap for patchworkplot.",
    )
    parser.add_argument(
        "--reverse-cmap",
        default="false",
        choices=["true", "false"],
        help="Whether to reverse the colormap.",
    )
    parser.add_argument(
        "--lwidth",
        type=int,
        default=2,
        help="Line width for patchworkplot.",
    )
    parser.add_argument(
        "--lower",
        action="store_true",
        help="Pass --lower to patchworkplot.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the command without executing it.",
    )
    return parser.parse_args()


def build_paths(dataset: str, locus: str, config: str | None, outdir: str | None) -> tuple[Path, Path]:
    if config is None:
        config_path = Path(f"data/processed/{dataset}/config/pwp_{locus}.tsv")
    else:
        config_path = Path(config)

    if outdir is None:
        outdir_path = Path(f"data/processed/{dataset}/{locus}/pwp")
    else:
        outdir_path = Path(outdir)

    return config_path, outdir_path


def validate_inputs(config_path: Path) -> None:
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    if not config_path.is_file():
        raise ValueError(f"Config path is not a file: {config_path}")

    if shutil.which("patchworkplot") is None:
        raise RuntimeError(
            "Executable 'patchworkplot' was not found in PATH. "
            "Please install it or activate the correct environment."
        )


def build_command(
    config_path: Path,
    outdir_path: Path,
    cmap: str,
    reverse_cmap: str,
    lwidth: int,
    lower: bool,
) -> list[str]:
    cmd = [
        "patchworkplot",
        "-i",
        str(config_path),
        "-o",
        str(outdir_path),
        "--cmap",
        cmap,
        "--reverse-cmap",
        reverse_cmap,
        "--lwidth",
        str(lwidth),
    ]

    if lower:
        cmd.append("--lower")

    return cmd


def main() -> None:
    args = parse_args()

    config_path, outdir_path = build_paths(
        dataset=args.dataset,
        locus=args.locus,
        config=args.config,
        outdir=args.outdir,
    )

    outdir_path.mkdir(parents=True, exist_ok=True)
    validate_inputs(config_path)

    cmd = build_command(
        config_path=config_path,
        outdir_path=outdir_path,
        cmap=args.cmap,
        reverse_cmap=args.reverse_cmap,
        lwidth=args.lwidth,
        lower=args.lower,
    )

    print("Running command:")
    print(" ".join(cmd))

    if args.dry_run:
        return

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()