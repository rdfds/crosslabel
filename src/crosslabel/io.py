"""Reading molecule maps.

Picasso HDF5 files hold a structured `locs` table with `x` and `y` in camera
pixels. One file per label (imaging round) is the usual RESI layout; the label
of each molecule is the index of its file.
"""
import csv

import numpy as np


def read_picasso(paths, pixel_nm):
    """Concatenate Picasso HDF5 molecule (or localization) tables, one per label.
    Returns positions in nm (N x 2) and labels (file index)."""
    try:
        import h5py
    except ImportError as e:   # pragma: no cover
        raise ImportError('reading Picasso files needs h5py: pip install "crosslabel[hdf5]"') from e
    xy, lab = [], []
    for i, path in enumerate(paths):
        with h5py.File(path, 'r') as f:
            locs = f['locs'][...]
        xy.append(np.column_stack([locs['x'], locs['y']]).astype(float) * pixel_nm)
        lab.append(np.full(len(locs), i))
    return np.concatenate(xy), np.concatenate(lab)


def read_table(path, x='x', y='y', label='label', scale=1.0):
    """Molecule positions and labels from a CSV with named columns. Labels may be
    any strings or integers; they are mapped to 0..K-1 in sorted order.
    `scale` converts positions to nm."""
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f'{path} has no rows')
    xy = np.array([[float(r[x]), float(r[y])] for r in rows]) * scale
    raw = [r[label] for r in rows]
    names = sorted(set(raw), key=lambda s: (len(s), s))
    index = {n: i for i, n in enumerate(names)}
    return xy, np.array([index[v] for v in raw]), names
