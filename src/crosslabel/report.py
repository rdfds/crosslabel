"""One-call analysis of a stochastically labelled molecule map.

`analyze(xy, labels)` runs every cross-label statistic on one point per
molecule (for example RESI or per-round clustered molecule positions) and
returns estimates, approximate 95% block-bootstrap intervals, descriptive
diagnostics and the assumptions needed for biological interpretation.
"""
import numpy as np

from scipy.spatial import cKDTree

from .correlation import (ANGLE_EDGES, EDGES, excess_neighbours, excess_with_ci, null_shifts, plateau_density,
                          point_pair_g, triplet_angles, triplet_m2)
from .labels import label_moments
from .stoichiometry import effective_size

LINEAR_FROM_DEG = 150.0
DUPLICATE_RADIUS = 10.0     # nm
DUPLICATE_WARN = 0.05       # same-label excess ordered pairs per input record
COINCIDENT_NM = 2.0
COINCIDENT_WARN = 0.01      # excess fraction of cross-label nearest neighbours within COINCIDENT_NM
MIN_BLOCKS = 16             # warning threshold, NOT a coverage guarantee above it


def duplicate_pairs(res, rmax=DUPLICATE_RADIUS):
    """Descriptive same-label excess ordered pairs per input record below rmax.

    Under exchangeable detection this can reveal splitting or merging, but
    label-dependent structure/detection can also change it. It is not a measured
    duplicate-molecule prevalence. The finite-count n_a^2 normalization and
    estimated cross-label background are retained for historical reproduction.
    """
    N, n = res['N'].sum(0), res['n'].sum(0)
    K = len(n)
    off = ~np.eye(K, dtype=bool)
    cross_rate = N[off].sum(0) / max(np.outer(n, n)[off].sum(), 1)
    same = np.array([N[a, a] for a in range(K)]).sum(0)
    inner = EDGES[1:] <= rmax
    return float((same[inner] - (n ** 2).sum() * cross_rate[inner]).sum() / max(n.sum(), 1))


def coincidences(xy, labels, rng, close=COINCIDENT_NM, n_shift=4, n_labels=None):
    """Close cross-label-neighbour fraction, NOT a double-recording rate bound.

    The denominator counts (input record, other occupied label) queries, not
    molecules. Each query contributes at most once. Duplicates displaced beyond
    close are missed; genuine neighbours and registration also affect the score.
    A quiet diagnostic does not establish label permanence.
    """
    if not np.isfinite(close) or close <= 0 or int(n_shift) != n_shift or n_shift < 1:
        raise ValueError('close must be positive and n_shift a positive integer')
    K = int(labels.max()) + 1 if n_labels is None else int(n_labels)
    sets = [xy[labels == a] for a in range(K)]

    def frac(shift):
        close_n = total = 0
        for a in range(K):
            for b in range(K):
                if a == b or not len(sets[a]) or not len(sets[b]):
                    continue
                pts = sets[b] if shift is None else sets[b] + shift[b] - shift[a]
                d, _ = cKDTree(pts).query(sets[a])
                close_n += int((d < close).sum()); total += len(d)
        return close_n, total

    hits, queries = frac(None)
    obs = hits / max(queries, 1)
    null = float(np.mean([h / max(n, 1) for h, n in
                         (frac(s) for s in null_shifts(rng, n_shift, n_labels=K))]))
    return dict(observed=obs, null=null, excess=obs - null, close_nm=float(close),
                close_queries=hits, queries=queries,
                denominator='input record x other occupied label',
                prevalence_bound=False, permanence_established=False)


def label_pair_excess(res, rho):
    """Labelled partners within 30 nm from each label pair separately. Under
    random labelling all pairs estimate the same quantity; a spread well beyond
    the intervals indicates label-dependent detection or registration."""
    K = res['pair_g'].shape[0]
    return {f'{a}-{b}': float(excess_neighbours(res['pair_g'][a, b], rho))
            for a in range(K) for b in range(a + 1, K) if res['n'][:, a].sum() and res['n'][:, b].sum()}


def analyze(xy, labels, efficiency=None, seed=0, B=400, angles=True, n_labels=None):
    """Cross-label analysis of molecule positions `xy` (N x 2, nm) with integer
    `labels` (0..K-1, the docking sequence or imaging round of each molecule).

    efficiency  labelling-and-detection efficiency p, if known (for example from
                a co-imaged standard via `assembly_stoichiometry`); enables the
                effective assembly size.
    Intervals are approximate block-bootstrap intervals, not exact label-redraw
    confidence sets. higher_order is a legacy alias for m2_positive_interval;
    it does not establish assembly identity or verify the model assumptions.
    """
    xy = np.asarray(xy, float)
    labels = np.asarray(labels)
    if xy.ndim != 2 or xy.shape[1] != 2 or labels.ndim != 1 or len(xy) != len(labels):
        raise ValueError('xy must be (N, 2) and labels (N,)')
    if not len(xy) or not np.isfinite(xy).all():
        raise ValueError('xy must be nonempty and finite')
    if not np.issubdtype(labels.dtype, np.number) or not np.isfinite(labels).all() or (labels < 0).any() or (labels != np.floor(labels)).any():
        raise ValueError('labels must be nonnegative integers')
    labels = labels.astype(int)
    if n_labels is not None and (int(n_labels) != n_labels or n_labels <= labels.max()):
        raise ValueError('n_labels must be an integer greater than every label')
    if len(np.unique(labels)) < 2:
        raise ValueError('at least two occupied labels are required')
    if int(B) != B or B < 2:
        raise ValueError('B must be an integer >= 2')
    if efficiency is not None and (not np.isfinite(efficiency) or not 0 < efficiency <= 1):
        raise ValueError('efficiency must be in (0, 1]')
    B = int(B)
    rng = np.random.default_rng(seed)
    K = int(labels.max()) + 1 if n_labels is None else int(n_labels)
    counts = np.bincount(labels, minlength=K)
    pi = counts / counts.sum()
    warnings = ['Biological interpretation is conditional on independent permanent labels, '
                'exchangeable detection, adequate background removal and appropriate spatial support. '
                'These diagnostics do not verify those assumptions.',
                'Intervals are approximate; sufficient block count alone does not guarantee coverage.']

    res = point_pair_g(xy, labels, rng, n_labels=K)
    rho = plateau_density(xy)
    n_blocks = int((res['n'].sum(1) > 0).sum())
    if n_blocks < MIN_BLOCKS:
        warnings.append(f'only {n_blocks} occupied 2-um blocks: intervals for m2 are too narrow in fields this small '
                        '(simulated coverage 0.76-0.84 with 4-9 blocks). More blocks do not guarantee coverage; '
                        'independent cells require a hierarchical analysis.')
    partners, partners_ci = excess_with_ci(res, rho, rng, B=B)
    out = dict(molecules=int(len(xy)), labels=int((counts > 0).sum()), label_fractions=pi.tolist(), blocks=n_blocks,
               density_per_um2=rho * 1e6, partners=partners, partners_ci=partners_ci,
               label_pair_partners=label_pair_excess(res, rho),
               duplicates_per_molecule=duplicate_pairs(res),
               coincident_cross_label=coincidences(xy, labels, rng, n_labels=K))
    out['records'] = len(xy)
    out['same_label_excess_pairs_per_record'] = out['duplicates_per_molecule']
    out['inference'] = dict(interval_method='approximate spatial block percentile bootstrap',
                            nominal_level=0.95, coverage_guaranteed=False,
                            assumptions_verified=False, physical_motion_validated=False,
                            weighting='input records; molecule weighting requires one record per molecule',
                            legacy_keys={'molecules': 'records',
                                         'duplicates_per_molecule': 'same_label_excess_pairs_per_record',
                                         'higher_order': 'm2_positive_interval'})
    if out['coincident_cross_label']['excess'] > COINCIDENT_WARN:
        warnings.append(f"{100 * out['coincident_cross_label']['excess']:.1f}% excess cross-label neighbours within "
                        f"{COINCIDENT_NM:g} nm: some molecules may be recorded under two labels (or partners are "
                        "closer than 2 nm). This is a fraction of record-label queries, not a prevalence bound.")
    if out['duplicates_per_molecule'] > DUPLICATE_WARN:
        warnings.append(f"{out['duplicates_per_molecule']:.2f} same-label excess ordered pairs per record at "
                        f"< {DUPLICATE_RADIUS:g} nm: compatible with duplicate detections, but not specific to them. "
                        "Same-label copies cannot form a three-distinct-molecule tuple; other tuples and weighting can change.")
    elif out['duplicates_per_molecule'] < -DUPLICATE_WARN:
        warnings.append(f"same-label pairs at < {DUPLICATE_RADIUS:g} nm are depleted "
                        f"({-out['duplicates_per_molecule']:.2f} ordered pairs per record): compatible with merging "
                        "or label-dependent detection. Merging can alter geometry and label selection; "
                        "cross-label statistics are not immune to it.")

    if (counts > 0).sum() >= 3:
        m2, m2_ci, p3 = triplet_m2(xy, labels, res, rng, B=B)
        out.update(m2=m2, m2_ci=m2_ci, distinct3_probability=p3, m2_per_partner=m2 / partners if partners > 0 else None)
        higher = m2_ci[0] > 0
        out['higher_order'] = bool(higher)
        out['m2_positive_interval'] = bool(higher)
        if efficiency is not None and higher and partners > 0:
            out['effective_size'] = effective_size(m2, partners, efficiency)
            out['effective_size_scope'] = ('model-conditional partner-weighted size point estimate; requires common '
                                           'independent efficiency and complete within-radius assemblies; '
                                           'not a minimum, unique complex size or confidence bound')
        if angles and higher:
            est, boot = triplet_angles(xy, labels, res, rng, B=B)
            lin = np.rad2deg((ANGLE_EDGES[:-1] + ANGLE_EDGES[1:]) / 2) > LINEAR_FROM_DEG
            totals = boot.sum(1)
            if not np.isfinite(est).all() or est.sum() <= 0 or not np.isfinite(boot).all() or (totals <= 0).any():
                warnings.append('angle normalization unsupported: nonpositive or nonfinite total in estimate/bootstrap; profile not reported.')
            else:
                norm = est / est.sum()
                nb = boot / totals[:, None]
                out.update(angle_edges_deg=np.rad2deg(ANGLE_EDGES).tolist(), angle_profile=norm.tolist(),
                           angle_profile_se=nb.std(0).tolist(), linearity=float(norm[lin].sum()),
                           linearity_ci=np.percentile(nb[:, lin].sum(1), [2.5, 97.5]).tolist(),
                           angle_scope='normalized signed connected profile, not an identified assembly geometry')
        elif angles:
            warnings.append('no significant connected triplets: angle profile not reported.')
    else:
        warnings.append('fewer than three labels in use: connected triplets need three distinct labels.')
    if efficiency is not None and partners > efficiency:
        warnings.append('partner point estimate exceeds the supplied efficiency: incompatible with the ideal '
                        'monomer/dimer expectation if the map, background and common-efficiency assumptions hold. '
                        'This point comparison is not a calibrated hypothesis test.')
    p = out['label_pair_partners']
    if len(p) > 1:
        out['label_pair_spread'] = float(np.max(list(p.values())) - np.min(list(p.values())))
    out['moments'] = {k: float(v) for k, v in label_moments(pi).items()}
    out['warnings'] = warnings
    return out


def format_report(r):
    """Human-readable summary of `analyze` output."""
    f = lambda v, ci: f'{v:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]'
    lines = [f"records {r['molecules']}  labels {r['labels']}  density {r['density_per_um2']:.1f} per um^2",
             f"labelled partners within 30 nm   {f(r['partners'], r['partners_ci'])}"]
    if 'm2' in r:
        lines.append(f"connected triplets m2            {f(r['m2'], r['m2_ci'])}"
                     + ('   positive approximate interval (model-conditional)' if r['higher_order'] else '   interval does not exclude zero'))
    if 'effective_size' in r:
        lines.append(f"partner-weighted effective size  {r['effective_size']:.1f} (conditional point estimate)")
    if 'linearity' in r:
        lines.append(f"linearity (angles >= 150 deg)    {f(r['linearity'], r['linearity_ci'])}")
    lines.append(f"same-label excess pairs / record {r['duplicates_per_molecule']:.3f}")
    c = r['coincident_cross_label']
    lines.append(f"close record-label queries < {c.get('close_nm', COINCIDENT_NM):g} nm "
                 f"{100 * c['observed']:.2f}% (null {100 * c['null']:.2f}%; not a prevalence bound)")
    if 'label_pair_spread' in r:
        lines.append(f"label-pair partners spread       {r['label_pair_spread']:.3f}")
    lines += [f'warning: {w}' for w in r['warnings']]
    return '\n'.join(lines)
