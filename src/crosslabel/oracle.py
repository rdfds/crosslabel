"""Oracle quantities (molecule identities known), used only to verify the exact
label identities on simulated data."""
import numpy as np
from scipy.spatial import cKDTree

from .labels import label_moments

def oracle_pairs(xy, mol, edges):
    """Ordered same-molecule (S) and between-molecule (W) pair histograms, plus
    per-molecule-pair counts within each bin for the exact label variance."""
    edges = np.asarray(edges, float)
    tree = cKDTree(xy)
    i, j = tree.query_pairs(edges[-1], output_type='ndarray').T
    d = np.hypot(*(xy[i] - xy[j]).T)
    b = np.digitize(d, edges) - 1
    ok = (b >= 0) & (b < len(edges) - 1)
    i, j, b = i[ok], j[ok], b[ok]
    same = mol[i] == mol[j]
    nb = len(edges) - 1
    S = 2 * np.bincount(b[same], minlength=nb).astype(float)
    W = 2 * np.bincount(b[~same], minlength=nb).astype(float)
    lo, hi = np.minimum(mol[i], mol[j])[~same], np.maximum(mol[i], mol[j])[~same]
    return S, W, (lo, hi, b[~same])


def exact_cross_variance(pairs, pi, nbins, nmol):
    """Exact label-draw variance of T_cross in each bin, given the configuration.

    T_cross = sum_{s<t} v_st (1 - 1[c_s = c_t]) with v_st the ordered pair count.
    Indicators of disjoint molecule pairs are independent; pairs sharing one
    molecule have covariance s3 - s2**2 (zero for uniform labels).
    """
    m = label_moments(pi)
    var1, cov1 = m['s2'] - m['s2'] ** 2, m['s3'] - m['s2'] ** 2
    lo, hi, b = pairs
    out = np.zeros(nbins)
    for k in range(nbins):
        sel = b == k
        if not sel.any():
            continue
        key = lo[sel] * nmol + hi[sel]
        u, v = np.unique(key, return_counts=True)
        v = 2.0 * v  # ordered pairs
        s, t = u // nmol, u % nmol
        total = var1 * (v ** 2).sum()
        if cov1 != 0:
            # sum over molecules of (sum_t v_st)^2 - sum_t v_st^2, i.e. pairs sharing s
            per = np.bincount(np.r_[s, t], weights=np.r_[v, v], minlength=nmol)
            per2 = np.bincount(np.r_[s, t], weights=np.r_[v, v] ** 2, minlength=nmol)
            total += cov1 * (per ** 2 - per2).sum()
        out[k] = total
    return out


def oracle_star(xy, mol, radius):
    """Star-triplet sum over three distinct molecules (labels ignored)."""
    tree = cKDTree(xy)
    nb = tree.query_ball_point(xy, radius)
    total = 0.0
    for i, js in enumerate(nb):
        js = np.asarray(js)
        js = js[mol[js] != mol[i]]
        if len(js) < 2:
            continue
        _, c = np.unique(mol[js], return_counts=True)
        total += c.sum() ** 2 - (c ** 2).sum()
    return total


def star_weights(xy, mol, radius):
    """Star triplets aggregated by molecule triple: for every ordered localization
    triple (centre i, neighbours j != k within `radius`, three distinct molecules)
    add one to the sorted molecule triple. Returns {(a, b, c): w}; sum(w) equals
    `oracle_star`."""
    tree = cKDTree(xy)
    w = {}
    for i, js in enumerate(tree.query_ball_point(xy, radius)):
        js = np.asarray(js, int)
        js = js[mol[js] != mol[i]]
        if len(js) < 2:
            continue
        u, c = np.unique(mol[js], return_counts=True)
        for x in range(len(u)):
            for y in range(x + 1, len(u)):
                key = tuple(sorted((int(mol[i]), int(u[x]), int(u[y]))))
                w[key] = w.get(key, 0.0) + 2.0 * c[x] * c[y]
    return w


def triple_label_covariances(pi):
    """Covariances of three-distinct-label indicators D_T = 1[labels of the three
    molecules in T are distinct] under independent labels ~ pi, by the number of
    molecules two triples share (3, 2, 1; disjoint triples are independent).

      shared 3:  D3 (1 - D3),                     D3 = 1 - 3 s2 + 2 s3
      shared 2:  sum_{x != y} pi_x pi_y (1 - pi_x - pi_y)^2 - D3^2
      shared 1:  sum_x pi_x g_x^2 - D3^2,         g_x = (1 - pi_x)^2 - (s2 - pi_x^2)

    For equiprobable labels the shared-1 covariance vanishes, as the shared-
    molecule covariance of pair indicators does."""
    pi = np.asarray(pi, float) / np.sum(pi)
    m = label_moments(pi)
    D3, s2 = m['distinct3'], m['s2']
    E2 = float((np.outer(pi, pi) * (1 - pi[:, None] - pi[None, :]) ** 2).sum() - (pi ** 2 * (1 - 2 * pi) ** 2).sum())
    g = (1 - pi) ** 2 - (s2 - pi ** 2)
    E1 = float((pi * g ** 2).sum())
    return D3 * (1 - D3), E2 - D3 ** 2, E1 - D3 ** 2


def exact_star_variance(weights, pi):
    """Exact label-draw variance of the three-distinct-label star sum
    T_star = sum_T w_T D_T given the configuration (positions, counts and
    molecule identities fixed; labels redrawn)."""
    c3, c2, c1 = triple_label_covariances(pi)
    if not weights:
        return 0.0
    pair, pair_sq, single = {}, {}, {}
    A3 = 0.0
    for (a, b, c), v in weights.items():
        A3 += v * v
        for p in ((a, b), (a, c), (b, c)):
            pair[p] = pair.get(p, 0.0) + v
            pair_sq[p] = pair_sq.get(p, 0.0) + v * v
        for s in (a, b, c):
            single[s] = single.get(s, 0.0) + v
    P2 = sum(pair[p] ** 2 - pair_sq[p] for p in pair)          # ordered triple pairs sharing exactly 2
    P1 = sum(v ** 2 for v in single.values()) - 3 * A3 - 2 * P2  # ... sharing exactly 1
    return float(c3 * A3 + c2 * P2 + c1 * P1)
