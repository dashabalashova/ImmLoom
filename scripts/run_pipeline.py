#!/usr/bin/env python3
import argparse
from pathlib import Path
from immloom.pipeline import run_pipeline

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-path", required=True)
    p.add_argument("--input-file", required=True)
    p.add_argument("--out-dir", default="./immloom_out")
    p.add_argument("--pi-min", default=80.0)
    p.add_argument("--segm-length-min", default=5000)
    p.add_argument("--dist-max", default=5000)
    p.add_argument("--block-length-min", default=5000)
    p.add_argument("--inv-plot", action="store_true", 
                   help="If set, produce inversion plot (default: false)")
    args = p.parse_args()
    run_pipeline(Path(args.input_path), args.input_file, Path(args.out_dir),
                 pi_min=args.pi_min, 
                 segm_length_min=args.segm_length_min, 
                 dist_max=args.dist_max, 
                 block_length_min=args.block_length_min, 
                 inv_plot=args.inv_plot)

if __name__ == "__main__":
    main()
