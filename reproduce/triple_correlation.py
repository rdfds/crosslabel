"""Coordinate triple correlation integrated over two anchor-neighbor discs."""

import numpy as np
from scipy.spatial import cKDTree


def disc_overlap(distance, radius):
    """Area shared by two equal discs at separations from 0 to 2R."""
    half = np.clip(distance / (2 * radius), 0, 1)
    return 2 * radius**2 * (np.arccos(half) - half * np.sqrt(1 - half**2))


def score(xy, codes, radius=30.0):
    """Return distinct-code triple and lower-order moments on one rectangle."""
    xy = np.asarray(xy, dtype=np.float64)
    codes = np.asarray(codes)
    if (xy.ndim != 2 or xy.shape[1] != 2 or not len(xy)
            or codes.shape != (len(xy),) or codes.dtype.kind not in "iu"
            or not np.isfinite(xy).all() or not np.isfinite(radius) or radius <= 0
            or (codes < 0).any() or (codes >= 4).any()):
        raise ValueError("Expected finite xy positions and four integer round codes")
    counts = np.bincount(codes, minlength=4)
    width = np.ptp(xy, axis=0)
    if (counts < 2).any() or (width <= 2 * radius).any():
        raise ValueError("Insufficient code or rectangular translation support")
    area = float(np.prod(width))
    disc_area = np.pi * radius**2
    denominator = sum(
        int(counts[a]) * int(counts[b]) * int(counts[c])
        for a in range(4) for b in range(4) for c in range(4)
        if len({a, b, c}) == 3
    )

    pairs = cKDTree(xy).query_pairs(2 * radius, output_type="ndarray")
    left, right = pairs.T
    keep = codes[left] != codes[right]
    left, right = left[keep], right[keep]
    delta = np.abs(xy[left] - xy[right])
    distance = np.linalg.norm(delta, axis=1)
    pair_overlap_fraction = np.prod(1 - delta / width, axis=1)
    other_codes = len(xy) - counts[codes[left]] - counts[codes[right]]
    pair_weight = 2 * other_codes / pair_overlap_fraction
    short = distance <= radius
    anchor_pair = float(area * pair_weight[short].sum() / (denominator * disc_area))
    neighbor_pair = float(
        area * np.dot(pair_weight, disc_overlap(distance, radius))
        / (denominator * disc_area**2)
    )

    # Build sorted adjacency from the same distinct-code short pair list.
    first = np.concatenate((left[short], right[short]))
    second = np.concatenate((right[short], left[short]))
    order = np.argsort(first, kind="stable")
    neighbors = second[order]
    offsets = np.r_[0, np.cumsum(np.bincount(first, minlength=len(xy)))]
    weighted_triples = 0.0
    for anchor in np.flatnonzero(np.diff(offsets) >= 2):
        near = neighbors[offsets[anchor]:offsets[anchor + 1]]
        i, j = np.triu_indices(len(near), k=1)
        distinct = codes[near[i]] != codes[near[j]]
        u, v = xy[near[i[distinct]]] - xy[anchor], xy[near[j[distinct]]] - xy[anchor]
        span = np.maximum(np.maximum(u, v), 0) - np.minimum(np.minimum(u, v), 0)
        weighted_triples += float(np.sum(1 / np.prod(1 - span / width, axis=1)))

    raw = 2 * (area / disc_area)**2 * weighted_triples / denominator
    connected = raw - 2 * anchor_pair - neighbor_pair + 2
    return {
        "dtc_raw": raw,
        "dtc_connected": connected,
        "anchor_pair": anchor_pair,
        "neighbor_pair": neighbor_pair,
        "dtc_connected_per_anchor": connected * (len(xy) * disc_area / area)**2 / 2,
    }
