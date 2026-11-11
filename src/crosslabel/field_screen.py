"""Coarse record density diagnostics for a rectangular coordinate extent."""
import numpy as np


def density_score(xy):
    """Mean of the squared count CVs on equal-area 4x4 and 8x8 grids."""
    points = np.asarray(xy, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not len(points) or not np.isfinite(points).all():
        raise ValueError("Expected finite two-dimensional coordinates")
    origin = points.min(axis=0)
    width = points.max(axis=0) - origin
    if not np.isfinite(width).all() or np.any(width <= 0):
        raise ValueError("The rectangle must have finite positive width and height")
    unit = (points - origin) / width
    scores = {}
    for divisions in (4, 8):
        cells = np.minimum((unit * divisions).astype(int), divisions - 1)
        counts = np.bincount(cells[:, 0] * divisions + cells[:, 1], minlength=divisions**2)
        scores[str(divisions)] = float(counts.var() / counts.mean()**2)
    return {"score": float(np.mean(list(scores.values()))), "cv_squared_by_grid": scores}


def reference_tail(score, controls):
    """Inclusive upper-tail rank against separately supplied reference fields."""
    values = np.asarray(controls, dtype=float)
    if values.ndim != 1 or len(values) < 19 or not np.isfinite(values).all():
        raise ValueError("At least 19 finite reference scores are required")
    if not np.isfinite(score):
        raise ValueError("The observed score must be finite")
    return (1 + int(np.sum(values >= score))) / (1 + len(values))
