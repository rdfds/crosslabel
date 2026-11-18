"""Compound-Poisson simulator used in the retention study."""
import numpy as np

SIDE = 10000.0
SIGMA = 5.024937810560445
RADIUS = 30.0

def shape(m, spacing, kind):
    if m == 1:
        return np.zeros((1, 2))
    if m == 2 or kind == 'line':
        return np.column_stack((spacing * (np.arange(m) - (m - 1) / 2), np.zeros(m)))
    a = 2 * np.pi * np.arange(m) / m
    r = spacing / (2 * np.sin(np.pi / m))
    return r * np.column_stack((np.cos(a), np.sin(a)))


def component_contribution(x, keep, radius=RADIUS):
    """Count observed own-component anchored pairs, independent of marks."""
    distance2 = np.sum((x[:, :, None, :] - x[:, None, :, :]) ** 2, axis=-1)
    n = ((distance2 <= radius ** 2) & keep[:, None, :]).sum(-1) - keep
    return int(np.sum(n * (n - 1) // 2 * keep))


def generate(c, fractions, rng):
    positions, marks = ([], [])
    own = 0
    for m, fraction in enumerate(fractions, 1):
        n = rng.poisson(100 * c['density'] * fraction / (m * c['p']))
        if n == 0:
            continue
        centers = SIDE * rng.random((n, 2))
        if c['shape'] == 'iid':
            offsets = rng.normal(0, c['spacing'], (n, m, 2))
        else:
            angle = 2 * np.pi * rng.random(n)
            cs, sn = (np.cos(angle)[:, None], np.sin(angle)[:, None])
            x, y = shape(m, c['spacing'], c['shape']).T
            offsets = np.stack((cs * x - sn * y, sn * x + cs * y), axis=2)
        xy = centers[:, None, :] + offsets + rng.normal(0, SIGMA, (n, m, 2))
        keep = rng.random((n, m)) < c['p']
        codes = rng.choice(c['K'], size=(n, m), p=c.get('codes'))
        own += component_contribution(xy, keep)
        positions.append(xy[keep])
        marks.append(codes[keep])
    xy, codes = (np.concatenate(positions), np.concatenate(marks))
    return (xy, codes, own / len(xy))
