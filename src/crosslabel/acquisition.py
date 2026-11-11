"""Acquisition-level raw-triplet moments under known independent retention.

This module separates detection+label randomization from the historical
label-only oracle conditional on an already detected configuration. It does
not supply end-to-end confidence intervals for background-corrected estimates.
The variance calculation specializes classical inclusion-probability sampling.
"""
import numpy as np
from .labels import label_moments
from .oracle import triple_label_covariances


def _overlaps(weights):
    pair, pair_sq, single = {}, {}, {}
    A3 = total = 0.0
    seen = set()
    for triple, weight in weights.items():
        if len(triple) != 3 or len(set(triple)) != 3:
            raise ValueError('each key must identify three distinct molecular sources')
        triple = tuple(sorted(triple))
        if triple in seen:
            raise ValueError('duplicate unordered source triple under permutation')
        seen.add(triple)
        weight = float(weight)
        if not np.isfinite(weight) or weight < 0:
            raise ValueError('raw-count weights must be finite and nonnegative')
        A3 += weight ** 2
        total += weight
        for p in ((triple[0], triple[1]), (triple[0], triple[2]), (triple[1], triple[2])):
            pair[p] = pair.get(p, 0.0) + weight
            pair_sq[p] = pair_sq.get(p, 0.0) + weight ** 2
        for source in triple:
            single[source] = single.get(source, 0.0) + weight
    P2 = sum(pair[p] ** 2 - pair_sq[p] for p in pair)
    P1 = sum(v ** 2 for v in single.values()) - 3 * A3 - 2 * P2
    return total, A3, P2, P1


def _moments(overlaps, pi, detection):
    pi = np.asarray(pi, float)
    if pi.ndim != 1 or not np.isfinite(pi).all() or (pi < 0).any() or pi.sum() <= 0:
        raise ValueError('pi must contain finite nonnegative label mass')
    if not np.isfinite(detection) or not 0 < detection <= 1:
        raise ValueError('known independent detection probability must be in (0, 1]')
    D = label_moments(pi)['distinct3']
    if D <= 0:
        raise ValueError('at least three nonzero labels are required')
    _, c2_label, c1_label = triple_label_covariances(pi)
    q = detection ** 3 * D
    coeff = np.array([q - q ** 2,
                      detection ** 4 * (c2_label + D ** 2) - q ** 2,
                      detection ** 5 * (c1_label + D ** 2) - q ** 2])
    total, A3, P2, P1 = overlaps
    variance = float(coeff @ [A3, P2, P1])
    return dict(inclusion_probability=q, raw_mean=q * total, raw_variance=variance,
                target_total=total, unbiased_total_variance=variance / q ** 2,
                relative_variance=variance / (q * total) ** 2 if total else None,
                overlap_sums=dict(A3=A3, P2=P2, P1=P1),
                normalized_coefficients=(coeff / q ** 2).tolist(),
                exact_scope='known fixed source weights, independent retention and permanent random labels',
                background_corrected_ci=False)


def triplet_sampling_moments(weights, pi, detection=1.0):
    """Exact mean/variance of a raw source-triplet sum and its inverse-inclusion estimate.

    weights maps unordered triples of known source IDs to fixed nonnegative
    weights. Unobserved source geometry/weights cannot be recovered by calling
    this function. detection is known independent per-source retention, including
    any explicitly assumed availability and acquisition detection factors.
    Its uncertainty, correlated losses, source-dependent kernels, background
    estimation and biological sampling variation are not included.
    """
    return _moments(_overlaps(weights), pi, detection)


def fixed_budget_precision(K_values, p_det, weights, availability=1.0, K_ref=4):
    """Compare expected triple yield with exact raw-total precision at fixed time.

    p_det(phi) gives an assumed per-source detection probability at relative
    round length phi=K_ref/K. All scenarios use the same known source weights.
    Predictions are conditional, not experimentally validated recommendations.
    """
    if not np.isfinite(availability) or not 0 < availability <= 1:
        raise ValueError('availability must be in (0, 1]')
    if int(K_ref) != K_ref or K_ref < 3:
        raise ValueError('K_ref must be an integer >= 3')
    overlaps = _overlaps(weights)
    if overlaps[0] <= 0:
        raise ValueError('a positive target total is required to compare precision')
    designs = {}
    for K in K_values:
        if int(K) != K or K < 2:
            raise ValueError('K values must be integers >= 2')
        K = int(K)
        if K < 3:
            designs[K] = dict(status='unsupported_fewer_than_three_labels')
            continue
        p = availability * float(p_det(K_ref / K))
        if p == 0:
            designs[K] = dict(status='unsupported_zero_detection')
            continue
        designs[K] = dict(status='model_prediction', detection=p,
                          **_moments(overlaps, np.ones(K), p))
    supported = [K for K, d in designs.items() if d['status'] == 'model_prediction']
    if not supported:
        raise ValueError('no supported design with at least three labels')
    ky = max(supported, key=lambda K: designs[K]['inclusion_probability'])
    kv = min(supported, key=lambda K: designs[K]['relative_variance'])
    return dict(designs=designs, yield_optimum=ky, variance_optimum=kv,
                predicted_relative_variance_reduction=(1 - designs[kv]['relative_variance'] /
                                                       designs[ky]['relative_variance']),
                empirical_validation=False,
                assumptions='equal round cost; known independent retention; fixed source geometry; no fitted background')


def geometry_free_precision(K_values, p_det, availability=1.0, K_ref=4):
    """Conservative raw-total variance screening without supplied source geometry.

    A design is weakly dominated if another has no larger coefficient for any
    of the three nonnegative overlap sums. Multiple retained designs mean the
    geometry-free rule has not selected an optimum. This sufficient partial
    order assumes known independent retention and uniform permanent labels;
    it is not a ratio-MSE or end-to-end spatial-inference rule.
    """
    base = fixed_budget_precision(K_values, p_det, {(0, 1, 2): 1.}, availability, K_ref)
    designs = {}
    for K, d in base['designs'].items():
        if d['status'] != 'model_prediction':
            designs[K] = dict(d)
            continue
        p = d['detection']
        D = (1 - 1/K) * (1 - 2/K)
        designs[K] = dict(status='known_model_coefficient_comparison', detection=p,
                          coefficients=[1/(p**3*D)-1, K/((K-1)*p*p)-1, (1-p)/p])
    supported = [K for K, d in designs.items() if 'coefficients' in d]
    for K in supported:
        b = np.asarray(designs[K]['coefficients'])
        dominated = []
        for other in supported:
            a = np.asarray(designs[other]['coefficients'])
            tolerance = 1e-12 * np.maximum(1., np.maximum(np.abs(a), np.abs(b)))
            if np.all(a <= b) and np.any(a < b - tolerance):
                dominated.append(other)
        designs[K]['weakly_dominated_by'] = dominated
    return dict(designs=designs,
                retained_designs=[K for K in supported if not designs[K]['weakly_dominated_by']],
                scope='sufficient geometry-free order for raw-total variance; not a ratio optimum',
                empirical_validation=False, assumptions=base['assumptions'].replace('fixed source geometry; ', ''))


def _ratio_moments(overlaps, n_sources, pi, detection):
    from scipy.stats import binom
    if int(n_sources) != n_sources or n_sources < 3:
        raise ValueError('full known source count must be an integer >= 3')
    # Shared validation, and the same label-only overlap probabilities.
    _moments(overlaps, pi, detection)
    n = int(n_sources)
    W, A3, P2, P1 = overlaps
    D = label_moments(pi)['distinct3']
    _, c2, c1 = triple_label_covariances(pi)
    E2, E1 = c2 + D ** 2, c1 + D ** 2
    m = np.arange(n + 1, dtype=float)
    prob = binom.pmf(m, n, detection)
    if not np.isclose(prob.sum(), 1., atol=1e-10):
        raise ValueError('binomial probability-mass calculation failed')
    inv = np.divide(1., m, out=np.zeros_like(m), where=m > 0)

    def f(r, power):
        if r > n:
            return 0.0
        inclusion = np.ones_like(m)
        for j in range(r):
            inclusion *= np.maximum(m - j, 0) / (n - j)
        return float(prob @ (inclusion * inv ** power))

    P0 = W ** 2 - A3 - P2 - P1
    mean = W * f(3, 1) / (2 * detection ** 2)
    second = (D * A3 * f(3, 2) + E2 * P2 * f(4, 2) + E1 * P1 * f(5, 2)
              + D ** 2 * P0 * f(6, 2)) / (4 * detection ** 4 * D ** 2)
    variance = second - mean ** 2
    if variance < -1e-9 * max(1., second):
        raise ValueError('negative ratio variance beyond numerical tolerance')
    variance = max(0., variance)
    target = W / (2 * n)
    bias = mean - target
    mse = variance + bias ** 2
    return dict(mean=mean, target=target, variance=variance, bias=bias, mse=mse,
                relative_mse=mse / target ** 2 if target else None,
                cov_raw_sum_and_retained_count=3 * detection ** 3 * D * (1 - detection) * W,
                exact_scope='known raw source-triple ratio with random retained count; not a spatial-estimator CI')


def triplet_ratio_moments(weights, pi, detection, n_sources):
    """Exact mean, variance and MSE of S/(2 p^2 D3 M), with value zero at M=0.

    n_sources is the full known source count, including sources absent from every
    weighted triple. The ratio is not generally unbiased; its reported bias is
    part of MSE. This is a raw geometric/source-group target. Estimated local
    backgrounds, unknown source counts, estimated retention and observation-map
    selection remain outside this exact calculation.
    """
    sources = set(source for triple in weights for source in triple)
    if n_sources < len(sources):
        raise ValueError('n_sources is smaller than the source IDs supplied')
    return _ratio_moments(_overlaps(weights), n_sources, pi, detection)


def fixed_budget_ratio_precision(K_values, p_det, weights, n_sources, availability=1.0, K_ref=4):
    """Known-model round design using exact raw-ratio MSE, including its denominator.

    The total-variance and ratio-MSE optima are separate questions. No fitted
    background or empirically validated acquisition improvement is implied.
    """
    base = fixed_budget_precision(K_values, p_det, weights, availability, K_ref)
    overlaps = _overlaps(weights)
    if n_sources < len(set(s for t in weights for s in t)):
        raise ValueError('n_sources is smaller than the supplied source set')
    designs = {}
    for K, d in base['designs'].items():
        if d['status'] != 'model_prediction':
            designs[K] = dict(d)
            continue
        designs[K] = dict(status='model_prediction', detection=d['detection'],
                          inclusion_probability=d['inclusion_probability'],
                          **_ratio_moments(overlaps, n_sources, np.ones(K), d['detection']))
    supported = [K for K, d in designs.items() if d['status'] == 'model_prediction']
    ky = base['yield_optimum']
    km = min(supported, key=lambda K: designs[K]['mse'])
    return dict(designs=designs, yield_optimum=ky, ratio_mse_optimum=km,
                predicted_relative_mse_reduction=1 - designs[km]['mse'] / designs[ky]['mse'],
                empirical_validation=False,
                assumptions=base['assumptions'], target='raw source-normalized ratio; includes finite ratio bias')
