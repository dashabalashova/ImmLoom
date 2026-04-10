from pathlib import Path
import pandas as pd


def parse_fasta_filename(path: Path):
    """
    Expect filename like:
        s01-h1-IGH.fasta

    Returns:
        sample_id: s01-h1
        locus: IGH
    """
    name = path.stem  # s01-h1-IGH
    parts = name.split("-")

    if len(parts) < 3:
        raise ValueError(f"Unexpected fasta filename format: {path}")

    sample_id = "-".join(parts[:2])  # s01-h1
    locus = parts[2]                # IGH

    return sample_id, locus


def make_relative(path: Path, root: Path):
    """Convert absolute path to relative path from project root"""
    try:
        return path.relative_to(root)
    except ValueError:
        return path


def generate_pwp_configs(
    dataset: str,
    project_root: Path = Path("."),
):
    """
    Generate pwp config TSV files for each locus.

    Output:
        data/processed/<dataset>/config/pwp_<LOCUS>.tsv

    Columns:
        SampleID, Label, Fasta
    """

    processed_dir = project_root / "data" / "processed" / dataset
    fasta_dir = processed_dir / "fasta" / "mixed"
    config_dir = processed_dir / "config"

    if not fasta_dir.exists():
        raise FileNotFoundError(f"Fasta directory not found: {fasta_dir}")

    config_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(fasta_dir.glob("*.fasta"))

    if not fasta_files:
        raise ValueError(f"No fasta files found in: {fasta_dir}")

    rows_by_locus = {}

    for fasta_path in fasta_files:
        sample_id, locus = parse_fasta_filename(fasta_path)

        fasta_rel = make_relative(fasta_path, project_root)

        row = {
            "SampleID": sample_id,
            "Label": sample_id,
            "Fasta": str(fasta_rel),
        }

        rows_by_locus.setdefault(locus, []).append(row)

    # write one file per locus
    for locus, rows in rows_by_locus.items():
        df = (
            pd.DataFrame(rows)
            .sort_values("SampleID")
            .reset_index(drop=True)
        )

        out_path = config_dir / f"pwp_{locus}.tsv"
        df.to_csv(out_path, sep="\t", index=False)

        print(f"Saved: {out_path}")
        print(df)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate pwp config files")
    parser.add_argument("--dataset", required=True, help="dataset name (e.g. dataset-02)")
    parser.add_argument("--project-root", default=".", help="project root path")

    args = parser.parse_args()

    generate_pwp_configs(
        dataset=args.dataset,
        project_root=Path(args.project_root),
    )


if __name__ == "__main__":
    main()