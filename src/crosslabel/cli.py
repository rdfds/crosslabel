"""Command line: crosslabel report ...

Examples
  crosslabel report --pixel-nm 130 R1.hdf5 R2.hdf5 R3.hdf5 R4.hdf5
  crosslabel report --csv molecules.csv --label-column round --efficiency 0.47 --json out.json
"""
import argparse
import json
import sys

from .io import read_picasso, read_table
from .report import analyze, format_report


def main(argv=None):
    ap = argparse.ArgumentParser(prog='crosslabel', description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='command', required=True)
    r = sub.add_parser('report', help='analyse one molecule map')
    r.add_argument('files', nargs='*', help='Picasso HDF5 files, one per label (round)')
    r.add_argument('--pixel-nm', type=float, help='camera pixel size for Picasso files (nm)')
    r.add_argument('--csv', help='CSV of molecules instead of Picasso files')
    r.add_argument('--x', default='x'); r.add_argument('--y', default='y')
    r.add_argument('--label-column', default='label')
    r.add_argument('--scale', type=float, default=1.0, help='multiply CSV positions by this to get nm')
    r.add_argument('--efficiency', type=float, help='labelling efficiency p, if known')
    r.add_argument('--bootstrap', type=int, default=400)
    r.add_argument('--seed', type=int, default=0)
    r.add_argument('--json', help='write the full result to this file')
    a = ap.parse_args(argv)

    if a.csv:
        xy, labels, names = read_table(a.csv, a.x, a.y, a.label_column, a.scale)
    elif a.files:
        if a.pixel_nm is None:
            ap.error('--pixel-nm is required for Picasso files')
        xy, labels = read_picasso(a.files, a.pixel_nm)
        names = list(a.files)
    else:
        ap.error('give Picasso files or --csv')
    res = analyze(xy, labels, efficiency=a.efficiency, seed=a.seed, B=a.bootstrap)
    res['label_names'] = names
    print(format_report(res))
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(res, fh, indent=1)
    return 0


if __name__ == '__main__':   # pragma: no cover
    sys.exit(main())
