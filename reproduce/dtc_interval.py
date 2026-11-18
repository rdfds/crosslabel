"""Full-field connected triple correlation with the Crosslabel interval."""
from unittest.mock import patch

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import t

from crosslabel import _field_kernel as kernel


def value_gradient(s, W, area, disc):
    """Differentiate the saved per-anchor connected triple-correlation score."""
    s = np.asarray(s, float)
    K = len(W)
    n = s[..., :K]
    q = s[..., K:-1].reshape(s.shape[:-1] + (K, K))
    triples = s[..., -1]
    count = n.sum(-1)
    denominator = np.einsum("abc,...a,...b,...c->...", W, n, n, n)
    if (n < 2).any() or (denominator <= 0).any():
        raise ValueError("Distinct-code and count support required")
    dP = (np.einsum("abc,...b,...c->...a", W, n, n)
          + np.einsum("abc,...a,...c->...b", W, n, n)
          + np.einsum("abc,...a,...b->...c", W, n, n))
    correction = np.einsum("abc,...a,...bc->...", W, n, q)
    factor = count**2 / denominator
    part = factor * (triples - .5 * correction)
    point = part + (count * disc / area)**2
    dn = part[..., None] * (2/count[..., None] - dP/denominator[..., None])
    dn -= .5 * factor[..., None] * np.einsum("abc,...bc->...a", W, q)
    dn += (2*count*(disc/area)**2)[..., None]
    dq = -.5 * factor[..., None, None] * np.einsum("abc,...a->...bc", W, n)
    gradient = np.concatenate((dn, dq.reshape(s.shape[:-1]+(K*K,)), factor[..., None]), -1)
    if not np.isfinite(point).all() or not np.isfinite(gradient).all():
        raise ValueError("Finite point and derivative required")
    return point, gradient


def sufficient(xy, codes, radius=30.):
    """Partition the exact translation-corrected score by anchor cell."""
    xy, codes, _ = kernel.validate_catalogue(xy, codes, 4, radius, "nm")
    K = 4
    span = np.ptp(xy, axis=0)
    area = float(np.prod(span))
    disc = np.pi * radius**2
    ij = np.floor((xy-xy.min(0))/kernel.CELL).astype(int)
    shape = tuple(ij.max(0)+1)
    cells = ij[:, 0]*shape[1]+ij[:, 1]
    total_cells = int(np.prod(shape))
    coarse_positions = np.unique(ij//20, axis=0)
    lookup = {tuple(p): i for i, p in enumerate(coarse_positions)}
    grid = np.indices(shape).reshape(2, -1).T
    mapping = np.array([lookup.get(tuple(p), -1) for p in grid//20])
    x = np.zeros((total_cells, K+K*K+1))
    np.add.at(x, (cells, codes), 1)
    q = np.zeros((total_cells, K, K))

    pairs = cKDTree(xy).query_pairs(2*radius, output_type="ndarray")
    left, right = pairs.T
    keep = codes[left] != codes[right]
    left, right = left[keep], right[keep]
    delta = np.abs(xy[left]-xy[right])
    distances = np.linalg.norm(delta, axis=1)
    exposure = np.prod(1-delta/span, axis=1)
    h = np.clip(distances/(2*radius), 0, 1)
    overlap = 2*radius**2*(np.arccos(h)-h*np.sqrt(1-h*h))
    short = distances <= radius
    weights = (2*disc*short+overlap)/(area*exposure)
    np.add.at(q, (cells[left], codes[left], codes[right]), weights)
    np.add.at(q, (cells[right], codes[right], codes[left]), weights)
    x[:, K:-1] = q.reshape(total_cells, K*K)

    first = np.r_[left[short], right[short]]
    second = np.r_[right[short], left[short]]
    neighbors = second[np.argsort(first, kind="stable")]
    offsets = np.r_[0, np.cumsum(np.bincount(first, minlength=len(xy)))]
    for anchor in np.flatnonzero(np.diff(offsets) >= 2):
        near = neighbors[offsets[anchor]:offsets[anchor+1]]
        i, j = np.triu_indices(len(near), k=1)
        distinct = codes[near[i]] != codes[near[j]]
        u = xy[near[i[distinct]]]-xy[anchor]
        v = xy[near[j[distinct]]]-xy[anchor]
        extent = np.maximum(np.maximum(u, v), 0)-np.minimum(np.minimum(u, v), 0)
        x[cells[anchor], -1] += np.sum(1/np.prod(1-extent/span, axis=1))
    mask = ~np.eye(K, dtype=bool)
    W = (mask[:, :, None] & mask[:, None, :] & mask[None, :, :]).astype(float)
    return x, mapping, shape, len(coarse_positions), W, area, disc


def spatial(x, W, shape, coarse_blocks, objective):
    point, gradient = objective(x.sum(0), W)
    scores = x @ gradient
    # For Crosslabel Euler's identity makes this mean zero already. The
    # fixed-window triple statistic is not homogeneous in its additive totals.
    scores -= scores.mean()
    variance, pairs = kernel.finite_lag_variance(scores, shape)
    if not np.isfinite(variance) or variance <= 0:
        return dict(status="unsupported", reason="nonpositive finite-lag variance")
    se = float(np.sqrt(variance))
    critical = float(t.ppf(.975, coarse_blocks-1))
    return dict(status="estimated", value=float(point), se=se, variance=variance,
                df=coarse_blocks-1, included_ordered_cell_pairs=pairs,
                spatial_t_ci=[float(point-critical*se), float(point+critical*se)])


def estimate(xy, codes, seed):
    x, mapping, shape, blocks, W, area, disc = sufficient(xy, codes)
    objective = lambda sums, weights: value_gradient(sums, weights, area, disc)
    spatial_result = spatial(x, W, shape, blocks, objective)
    coarse = kernel.coarse_sum(x, mapping, blocks)
    with patch.object(kernel, "value_gradient", objective):
        result, arrays = kernel.nested_estimate(coarse, W, spatial_result, seed)
    return result, arrays, dict(fine=x, coarse=coarse, W=W, area=area, disc=disc,
                               shape=shape, coarse_blocks=blocks)
