"""Observation models and record-level targets for the interval benchmark."""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from crosslabel.simulate import acquire, complexes, smooth_field


CONDITIONS = (
    "sparse_dimer", "sparse_trimer", "dense_dimer", "dense_four",
    "raw_dimer", "raw_trimer",
)


def schedule(repeats):
    for index in range(1, repeats + 1):
        for base, condition in zip(range(301000, 307000, 1000), CONDITIONS):
            yield base + index, condition, condition.startswith("dense") and index <= 10


def single_record_target(xy, components, radius):
    neighbors = cKDTree(xy).query_ball_point(xy, radius)
    counts = np.array([
        np.sum(components[indices] == components[i]) - 1
        for i, indices in enumerate(neighbors)
    ])
    return float((counts * (counts - 1) / 2).mean())


def repeated_record_target(xy, molecules, components, radius):
    """Count within-component stars whose three records have distinct sources."""
    order = np.argsort(components, kind="stable")
    sorted_components = components[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_components)) + 1, len(components)]
    total = 0.
    for low, high in zip(starts[:-1], starts[1:]):
        indices = order[low:high]
        points, sources = xy[indices], molecules[indices]
        distinct = np.unique(sources)
        if len(distinct) < 3:
            continue
        near = ((points[:, None] - points[None, :]) ** 2).sum(-1) <= radius ** 2
        counts = np.column_stack([
            near[:, sources == source].sum(1) * (sources != source)
            for source in distinct
        ])
        total += (counts.sum(1) ** 2 - (counts ** 2).sum(1)).sum()
    return float(total / (2 * len(xy)))


def sparse_field(rng, fractions):
    """Generate one record per retained molecule in the 12 nm ring model."""
    mixture = np.asarray(fractions, float) / np.arange(1, len(fractions) + 1)
    mixture /= mixture.sum()
    points, components = complexes(rng, 10000., 70. / .47, mixture, 12., shape="ring")
    data = acquire(rng, points, components, 4, efficiency=.47, link_sd=5.,
                   locs_per_molecule=1, loc_sd=0.)
    xy = data["site"] + rng.normal(0, .5, data["site"].shape)
    return xy, data["site_label"], single_record_target(xy, data["site_complex"], 30.)


def field(seed, condition):
    rng = np.random.default_rng(seed)
    if condition.startswith("dense"):
        intensity = smooth_field(rng, 10000., n_bumps=8, scale=1500., floor=.4)
        dimer = condition == "dense_dimer"
        fractions = [.5, .5] if dimer else [.6, 0, 0, .4]
        points, components = complexes(
            rng, 10000., 1300. / .5 / .7, fractions, 3. if dimer else 2.,
            intensity=intensity, shape="chain" if dimer else "blob",
        )
        data = acquire(
            rng, points, components, 6, efficiency=.5,
            pi=[.18, .20, .23, .05, .17, .17], locs_per_molecule=1, loc_sd=0.,
        )
        xy = data["site"] + rng.normal(0, 1., data["site"].shape)
        return xy, data["site_label"], single_record_target(xy, data["site_complex"], 10.), 10.

    fractions = np.array([.5, .2, .3] if condition.endswith("trimer") else [.5, .5])
    if condition.startswith("sparse"):
        return (*sparse_field(rng, fractions), 30.)
    mixture = fractions / np.arange(1, len(fractions) + 1)
    mixture /= mixture.sum()
    points, components = complexes(rng, 10000., 70. / .47, mixture, 12., shape="ring")
    data = acquire(
        rng, points, components, 4, efficiency=.47, link_sd=5.,
        locs_per_molecule=20., loc_sd=4.,
    )
    xy = data["xy"]
    return xy, data["label"], repeated_record_target(xy, data["mol"], data["complex"], 30.), 30.


def merge_same_code(xy, codes, distance=10.):
    merged, labels = [], []
    for code in np.unique(codes):
        points = xy[codes == code]
        pairs = cKDTree(points).query_pairs(distance, output_type="ndarray")
        graph = coo_matrix(
            (np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
            shape=(len(points), len(points)),
        )
        _, groups = connected_components(graph, directed=False)
        counts = np.bincount(groups)
        merged.append(np.column_stack([
            np.bincount(groups, column) / counts for column in points.T
        ]))
        labels.append(np.full(len(counts), code))
    return np.concatenate(merged), np.concatenate(labels)
