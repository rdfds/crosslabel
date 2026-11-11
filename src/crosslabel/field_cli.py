"""Command line for the manuscript field estimator."""
import argparse
import json
from pathlib import Path

from .field import field_check
from .io import read_picasso, read_table


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--csv", type=Path, help="Table with x_nm,y_nm,code columns")
    input_group.add_argument("--picasso", nargs="+", type=Path, help="One HDF5 file per code")
    parser.add_argument("--pixel-nm", type=float, help="Picasso camera pixel size in nm")
    parser.add_argument("--radius-nm", type=float, required=True, choices=(10., 30.))
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--reference-scores", type=Path, help="JSON array of independent density scores")
    parser.add_argument("--output", type=Path, help="Write a new JSON report")
    args = parser.parse_args(argv)
    if args.csv is not None and args.pixel_nm is not None:
        parser.error("--pixel-nm is only for Picasso; CSV coordinates must already be in nm")
    if args.picasso is not None and (args.pixel_nm is None or not 0 < args.pixel_nm < float("inf")):
        parser.error("Picasso input requires a finite positive --pixel-nm")
    try:
        if args.csv is not None:
            xy, codes, names = read_table(args.csv, x="x_nm", y="y_nm", label="code")
        else:
            xy, codes = read_picasso(args.picasso, args.pixel_nm)
            names = [p.name for p in args.picasso]
        references = None if args.reference_scores is None else json.loads(args.reference_scores.read_text())
        report = field_check(xy, codes, code_count=len(names), radius_nm=args.radius_nm,
                             seed=args.seed, reference_scores=references)
        report["code_names"] = names
        payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output is None:
            print(payload, end="")
        else:
            with args.output.open("x") as f:
                f.write(payload)
    except (ValueError, KeyError, OSError, ImportError, OverflowError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
