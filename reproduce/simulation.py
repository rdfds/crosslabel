"""Simulate the marked fields used in the held-out comparisons."""

import numpy as np

import spinna_design as old


BASE = dict(p=.47, density=55., K=4, codes=(.15, .26, .48, .11))


NEW = {
    "matched_moments": dict(BASE, shape="iid", spacing=7., mask=False, null=2/3, alt=(7/12, 2/12, 3/12)),
    "membrane_variation": dict(BASE, shape="line", spacing=13.5, mask=True, mask_pattern="two_dimensional",
                               null=.47, alt=(.7, .2, .1)),
}


def generate(c, fractions, rng):
    positions, marks = [], []
    own = 0
    if c["mask"]:
        pmf = old.mask_pmf()
        if c.get("mask_pattern") == "two_dimensional":
            x = 0.5 * (np.arange(20) + .5)
            cos = np.cos(2 * np.pi * x / 5)
            pmf = 1 + .7 * np.outer(cos, cos)
        pmf = pmf.ravel() / pmf.sum()
    for m, fraction in enumerate(fractions, 1):
        n = rng.poisson(100 * c["density"] * fraction / (m * c["p"]))
        if n == 0:
            continue
        if c["mask"]:
            bins = rng.choice(400, size=n, p=pmf)
            centers = 500 * (np.column_stack((bins % 20, bins // 20)) + rng.random((n, 2)))
        else:
            centers = old.SIDE * rng.random((n, 2))
        if c["shape"] == "iid":
            offsets = rng.normal(0, c["spacing"], (n, m, 2))
        else:
            angle = 2 * np.pi * rng.random(n)
            cs, sn = np.cos(angle)[:, None], np.sin(angle)[:, None]
            x, y = old.shape(m, c["spacing"], c["shape"]).T
            offsets = np.stack((cs*x - sn*y, sn*x + cs*y), axis=2)
        xy = centers[:, None, :] + offsets + rng.normal(0, old.SIGMA, (n, m, 2))
        keep = rng.random((n, m)) < c["p"]
        codes = rng.choice(c["K"], size=(n, m), p=c.get("codes"))
        own += old.component_contribution(xy, keep)
        positions.append(xy[keep])
        marks.append(codes[keep])
    xy, codes = np.concatenate(positions), np.concatenate(marks)
    return xy, codes, own / len(xy)


def restricted(q):
    return 1 - q, q, 0.
