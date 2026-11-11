"""Fitted marked radial pair summaries for the fixed component libraries."""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_triangular
from scipy.spatial import cKDTree

EDGES = np.array([0, 5, 10, 15, 20, 25, 30, 40, 60, 80, 100], dtype=float)


def summary(xy, codes, code_count):
    """Counts and the complete unordered code-pair by distance-bin histogram."""
    xy = np.asarray(xy, dtype=float)
    codes = np.asarray(codes)
    if (
        xy.shape != (len(codes), 2)
        or len(codes) == 0
        or codes.ndim != 1
        or codes.dtype.kind not in "iu"
        or not np.isfinite(xy).all()
        or not np.all((codes >= 0) & (codes < code_count))
    ):
        raise ValueError("Expected a nonempty, finite, coded 2D catalogue")
    left, right = np.triu_indices(code_count)
    pair_index = np.empty((code_count, code_count), dtype=int)
    pair_index[left, right] = np.arange(len(left))
    pair_index[right, left] = np.arange(len(left))
    pairs = cKDTree(xy).query_pairs(EDGES[-1], output_type="ndarray")
    i, j = pairs.T
    d2 = np.sum((xy[i] - xy[j]) ** 2, axis=1)
    n_bins = len(EDGES) - 1
    # Include 100 nm in the final bin. A duplicated position is in the first.
    bins = np.minimum(np.searchsorted(EDGES**2, d2, side="right") - 1, n_bins - 1)
    hist = np.bincount(
        n_bins * pair_index[codes[i], codes[j]] + bins,
        minlength=len(left) * n_bins,
    )
    return np.concatenate((np.bincount(codes, minlength=code_count), hist))


@dataclass(frozen=True)
class Library:
    origin: np.ndarray
    direction: np.ndarray
    factor: np.ndarray

    @classmethod
    def train(cls, monomer, dimer):
        monomer = np.asarray(monomer, dtype=float)
        dimer = np.asarray(dimer, dtype=float)
        if monomer.ndim != 2 or dimer.shape != monomer.shape or len(monomer) < 2:
            raise ValueError("Expected matching independent endpoint banks")
        m0, m1 = monomer.mean(0), dimer.mean(0)
        centered = np.vstack((monomer - m0, dimer - m1))
        covariance = centered.T @ centered / (len(centered) - 2)
        target = np.diag(np.maximum(np.diag(covariance), 1.0))
        factor = np.linalg.cholesky(0.5 * covariance + 0.5 * target)
        direction = solve_triangular(factor, m1 - m0, lower=True)
        if not np.isfinite(direction).all() or np.dot(direction, direction) <= 0:
            raise ValueError("Endpoints do not identify the fitted direction")
        return cls(m0, direction, factor)

    def fit(self, values):
        values = np.atleast_2d(np.asarray(values, dtype=float))
        if not np.isfinite(values).all():
            raise ValueError("Expected finite marked-pair summaries")
        whitened = solve_triangular(self.factor, (values - self.origin).T, lower=True).T
        q = np.clip(whitened @ self.direction / np.dot(self.direction, self.direction), 0, 1)
        residual = whitened - q[:, None] * self.direction
        return q, np.sum(residual**2, axis=1)


def rank_grid(score, grid_scores):
    grid_scores = np.asarray(grid_scores, dtype=float)
    if grid_scores.ndim != 2 or not grid_scores.size or not np.isfinite(grid_scores).all():
        raise ValueError("Expected a nonempty, finite score bank at every grid point")
    if not np.isfinite(score):
        raise ValueError("Expected a finite observed discrepancy")
    return (1 + (grid_scores >= score).sum(axis=1)) / (1 + grid_scores.shape[1])
