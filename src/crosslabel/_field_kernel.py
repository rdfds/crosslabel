"""Geometry, sufficient statistics, and interval kernels for the field check."""
import hashlib
import math

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import t as student_t

from . import correlation as C
from .labels import star_counts

CELL = 100.
DEPENDENCE = 600.
OUTER = 999
INNER = 199
BATCH = 16
MAX_RECORDS = 5_000_000
MAX_CELLS = 1_000_000
MAX_COARSE = 2048



# Disc overlap


def disc_overlap(separation, r, s):
    z = np.asarray(separation, float)
    out = np.zeros_like(z)
    enclosed = z <= abs(r-s)
    out[enclosed] = np.pi * min(r,s)**2
    partial = (z > abs(r-s)) & (z < r+s)
    d = z[partial]
    if len(d):
        x = np.clip((d*d+r*r-s*s)/(2*d*r), -1., 1.)
        y = np.clip((d*d+s*s-r*r)/(2*d*s), -1., 1.)
        q = np.maximum((-d+r+s)*(d+r-s)*(d-r+s)*(d+r+s), 0.)
        out[partial] = r*r*np.arccos(x)+s*s*np.arccos(y)-.5*np.sqrt(q)
    return out


# Coarse-block standard errors


def block_se(x, gradients, multiplicities=None):
    x = np.asarray(x, float)
    nb = len(x)
    if nb < 2:
        raise ValueError('at least two spatial blocks required')
    scores = np.asarray(gradients) @ x.T
    if multiplicities is None:
        centered = scores - scores.mean(-1, keepdims=True)
        variance = nb / (nb - 1) * (centered ** 2).sum(-1)
    else:
        m = np.asarray(multiplicities, float)
        if m.shape != scores.shape or (m < 0).any() or not np.all(m.sum(-1) == nb):
            raise ValueError('multiplicities must be nonnegative and sum to block count')
        centered = scores - (m * scores).sum(-1, keepdims=True) / nb
        variance = nb / (nb - 1) * (m * centered ** 2).sum(-1)
    se = np.sqrt(variance)
    if not np.isfinite(se).all() or (se <= 0).any():
        raise ValueError('unusable block standard error; no draws discarded')
    return se


# Pair compensation and derivatives


def compensated_weights(radius, annulus=(100., 250.)):
    mid = (C.EDGES[:-1] + C.EDGES[1:]) / 2
    area = np.pi * np.diff(C.EDGES ** 2)
    lo, hi = annulus
    if not 0 < radius <= lo < hi:
        raise ValueError('positive radius and disjoint ordered annulus required')
    ratio = radius ** 2 / (hi ** 2 - lo ** 2)
    dd = disc_overlap(mid, radius, radius)
    aa = disc_overlap(mid, hi, hi) - 2 * disc_overlap(mid, hi, lo) + disc_overlap(mid, lo, lo)
    da = disc_overlap(mid, radius, hi) - disc_overlap(mid, radius, lo)
    kernel = dd + ratio ** 2 * aa - 2 * ratio * da
    base = (mid >= C.BASELINE[0]) & (mid < C.BASELINE[1])
    # Assumption: connected pair excess beyond the supplied100nm support is
    # zero. The baseline removes the constant density term, not unknown tails.
    weights = kernel.copy()
    weights[base] -= np.sum(kernel * area) / area[base].sum()
    return weights


def value_gradient(s, W):
    s = np.asarray(s, float); W = np.asarray(W, float); K = len(W)
    if W.shape != (K, K, K) or s.shape[-1] != 2 * K + K * K + 1 or not np.isfinite(s).all():
        raise ValueError('incompatible or nonfinite sufficient vector')
    n = s[..., :K]
    Q = s[..., K:K + K * K].reshape(s.shape[:-1] + (K, K))
    A = s[..., K + K * K]
    I = s[..., K + K * K + 1:]
    count = n.sum(-1)
    P = np.einsum('abc,...a,...b,...c->...', W, n, n, n)
    if (n < 2).any() or (P <= 0).any():
        raise ValueError('at least two records per label and distinct-triple support required')
    factor = count ** 2 / (2 * P)
    correction = np.einsum('abc,...a,...bc->...', W, I, Q) / count
    point = factor * (A - correction)
    dP = (np.einsum('abc,...b,...c->...a', W, n, n) +
          np.einsum('abc,...a,...c->...b', W, n, n) +
          np.einsum('abc,...a,...b->...c', W, n, n))
    dn = point[..., None] * (2 / count[..., None] - dP / P[..., None])
    dn += (factor * correction / count)[..., None]
    dq = -(factor / count)[..., None, None] * np.einsum('abc,...a->...bc', W, I)
    di = -(factor / count)[..., None] * np.einsum('abc,...bc->...a', W, Q)
    gradient = np.concatenate([dn, dq.reshape(s.shape[:-1] + (K * K,)), factor[..., None], di], axis=-1)
    if not np.isfinite(point).all() or not np.isfinite(gradient).all():
        raise ValueError('nonfinite density-repaired point or gradient')
    return point, gradient


# Window geometry


def quadrant_area(x, y, radius):
    """Area of a radius-r disc in [0,x] x [0,y], x/y nonnegative."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(x).all() or not np.isfinite(y).all() or (x < 0).any() or (y < 0).any():
        raise ValueError('finite nonnegative quadrant extents and positive radius required')
    x = np.minimum(x, radius); y = np.minimum(y, radius)
    inside = x*x + y*y <= radius*radius
    cutoff = np.sqrt(np.maximum(radius*radius-y*y, 0.))
    primitive_x = .5 * (x*np.sqrt(np.maximum(radius*radius-x*x, 0.)) +
                        radius*radius*np.arcsin(np.clip(x/radius, 0., 1.)))
    primitive_c = .5 * (cutoff*y + radius*radius*np.arcsin(np.clip(cutoff/radius, 0., 1.)))
    return np.where(inside, x*y, y*cutoff + primitive_x - primitive_c)


def circle_rectangle_area(xy, radius, window):
    """Exact analytic disc intersection with window=(xmin,xmax,ymin,ymax)."""
    xy = np.asarray(xy, float)
    window = np.asarray(window, float)
    if xy.shape[-1] != 2 or window.shape != (4,) or not np.isfinite(xy).all() or not np.isfinite(window).all():
        raise ValueError('finite centres and four window bounds required')
    xmin, xmax, ymin, ymax = window
    if xmin >= xmax or ymin >= ymax: raise ValueError('window must have positive width and height')
    def signed(x, y):
        return np.sign(x)*np.sign(y)*quadrant_area(np.abs(x), np.abs(y), radius)
    x0, x1 = xmin-xy[...,0], xmax-xy[...,0]
    y0, y1 = ymin-xy[...,1], ymax-xy[...,1]
    result = signed(x1,y1)-signed(x0,y1)-signed(x1,y0)+signed(x0,y0)
    # No clipping conceals a geometry error; allow floating roundoff in checks.
    tolerance = 1e-10*max(np.pi*radius*radius, 1.)
    if (result < -tolerance).any() or (result > np.pi*radius*radius+tolerance).any():
        raise ValueError('invalid circle/window area')
    return result


def radial_pair_exposure(edges, width, height):
    """Exact annular-average overlap fraction, provided rmax<=min(width,height)."""
    edges = np.asarray(edges, float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or (np.diff(edges) <= 0).any() or edges[0] < 0:
        raise ValueError('ordered finite nonnegative radial edges required')
    if not np.isfinite([width,height]).all() or min(width,height) <= 0 or edges[-1] > min(width,height):
        raise ValueError('positive window dimensions and supported radii required')
    a, b = edges[:-1], edges[1:]
    mean_r = (2./3.)*(b**3-a**3)/(b*b-a*a)
    mean_r2 = .5*(b**4-a**4)/(b*b-a*a)
    return 1 - (2/np.pi)*(1/width+1/height)*mean_r + mean_r2/(np.pi*width*height)


# Annulus covariance


def covariance_amplitude(disc_area,annulus_area,full_disc,full_annulus):
    disc_area=np.asarray(disc_area,float);annulus_area=np.asarray(annulus_area,float)
    if not np.isfinite(disc_area).all() or not np.isfinite(annulus_area).all() or (disc_area<=0).any() or (annulus_area<=0).any():
        raise ValueError('positive finite effective disc and annulus areas required')
    return (disc_area+disc_area**2/annulus_area)/(full_disc+full_disc**2/full_annulus)


# Spatial sufficient statistics


def fine_sufficient(xy,labels,radius,K,annulus=(100.,250.)):
    xy=np.asarray(xy,float);labels=np.asarray(labels)
    if xy.ndim!=2 or xy.shape[1]!=2 or len(xy)!=len(labels) or not len(xy) or not np.isfinite(xy).all():
        raise ValueError('matching nonempty finite coordinates/labels required')
    if not np.issubdtype(labels.dtype,np.integer) or (labels<0).any() or (labels>=K).any():
        raise ValueError('integer labels within declared support required')
    window=np.array([xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max()])
    width,height=window[1]-window[0],window[3]-window[2]
    if min(width,height)<=2*annulus[1]:raise ValueError('assumed footprint too small')
    ij=np.floor((xy-xy.min(0))/CELL).astype(int)
    shape=tuple(int(v) for v in ij.max(0)+1);nb=int(np.prod(shape));cell=ij[:,0]*shape[1]+ij[:,1]
    coarse_positions,coarse=np.unique(ij//int(C.BLOCK/CELL),axis=0,return_inverse=True)
    n=np.zeros((nb,K));np.add.at(n,(cell,labels),1)
    d=np.pi*radius**2;a=np.pi*(annulus[1]**2-annulus[0]**2)
    Di=circle_rectangle_area(xy,radius,window)
    Ai=circle_rectangle_area(xy,annulus[1],window)-circle_rectangle_area(xy,annulus[0],window)
    if (Di<=0).any() or (Ai<=0).any():raise ValueError('positive observed areas required')
    factor=covariance_amplitude(Di,Ai,d,a)
    ann=star_counts(xy,labels,annulus[1],K)-star_counts(xy,labels,annulus[0],K)
    disc=star_counts(xy,labels,radius,K)
    lam=ann/Ai[:,None];lam_old=ann/a
    e=disc-lam*Di[:,None];e_old=disc-lam_old*d
    I=np.zeros((nb,K));Iold=np.zeros_like(I)
    np.add.at(I,(cell,labels),lam.sum(1)*factor)
    np.add.at(Iold,(cell,labels),lam_old.sum(1))
    weights=compensated_weights(radius,annulus)
    exposure=radial_pair_exposure(C.EDGES,width,height)
    if (exposure<=0).any():raise ValueError('positive pair exposure required')
    Q=np.zeros((nb,K,K));Qold=np.zeros_like(Q)
    tree=cKDTree(xy);start=0;pair_count=0
    # Bound temporary neighbour materialization, never subsample a record.
    while start<len(xy):
        end=min(start+512,len(xy))
        while True:
            counts=tree.query_ball_point(xy[start:end],C.EDGES[-1],return_length=True)
            if int(counts.sum())<=262144:break
            if end-start==1:raise ValueError('one-anchor neighbour list exceeds fixed memory support')
            end=start+max(1,(end-start)//2)
        lists=tree.query_ball_point(xy[start:end],C.EDGES[-1])
        ii=np.repeat(np.arange(start,end),[len(v) for v in lists])
        jj=np.concatenate([np.asarray(v,dtype=int) for v in lists])
        keep=ii!=jj;ii=ii[keep];jj=jj[keep]
        distance=np.linalg.norm(xy[ii]-xy[jj],axis=1)
        shell=np.maximum(0,np.searchsorted(C.EDGES,distance,side='left')-1)
        keep=shell<len(weights);ii=ii[keep];jj=jj[keep];shell=shell[keep]
        np.add.at(Q,(cell[ii],labels[ii],labels[jj]),weights[shell]/exposure[shell])
        np.add.at(Qold,(cell[ii],labels[ii],labels[jj]),weights[shell])
        pair_count+=len(ii);start=end
    sep=np.abs(np.subtract.outer(np.arange(K),np.arange(K)));results={}
    for name,mask in [('all',sep>=1),('non_consecutive',sep>=2)]:
        W=(mask[:,:,None]&mask[:,None,:]&mask[None,:,:]).astype(float)
        A=np.bincount(cell,weights=np.einsum('ib,ibc,ic->i',e,W[labels],e),minlength=nb)
        Aold=np.bincount(cell,weights=np.einsum('ib,ibc,ic->i',e_old,W[labels],e_old),minlength=nb)
        Z=np.zeros((len(coarse_positions),K,K))
        for b in range(len(coarse_positions)):
            select=coarse==b
            Z[b]=np.einsum('ib,ic,ibc->bc',lam_old[select],lam_old[select],W[labels[select]])
        results[name]=dict(x=np.column_stack([n,Q.reshape(nb,K*K),A[:,None],I]),
                           x_s2=np.column_stack([n,Qold.reshape(nb,K*K),Aold[:,None],Iold]),W=W,Z=Z)
    # Coarse-map contribution aggregation is exact; drop only empty *cells*,
    # never records, when matching the original occupied-block definition.
    grid=np.indices(shape).reshape(2,-1).T
    fine_coarse=grid//int(C.BLOCK/CELL)
    lookup={tuple(p):i for i,p in enumerate(coarse_positions)}
    mapping=np.array([lookup.get(tuple(p),-1) for p in fine_coarse])
    geometry=dict(window=window.tolist(),shape=list(shape),cell_nm=CELL,dependence_nm=DEPENDENCE,
                  records_retained=len(xy),ordered_pairs_enumerated=pair_count,
                  coarse_blocks=len(coarse_positions),true_footprint_verified=False,
                  finite_separation_covariance_approximation=True)
    return results,mapping,geometry


def coarse_sum(x,mapping,nb):
    out=np.zeros((nb,x.shape[1]));valid=mapping>=0
    if np.any(x[~valid]!=0):raise ValueError('nonzero contribution outside occupied coarse support')
    np.add.at(out,mapping[valid],x[valid]);return out


def finite_lag_variance(scores,shape,cell=CELL,dependence=DEPENDENCE):
    scores=np.asarray(scores,float)
    if not np.isfinite(scores).all() or tuple(shape)==() or scores.size!=int(np.prod(shape)):
        raise ValueError('finite scores on the complete declared grid required')
    y=scores.reshape(shape);value=0.;included=0
    limit=int(np.ceil(dependence/cell))+1
    for dx in range(-limit,limit+1):
        for dy in range(-limit,limit+1):
            gap=cell*np.hypot(max(abs(dx)-1,0),max(abs(dy)-1,0))
            if gap>dependence:continue
            xlo=max(0,-dx);xhi=min(shape[0],shape[0]-dx)
            ylo=max(0,-dy);yhi=min(shape[1],shape[1]-dy)
            if xlo>=xhi or ylo>=yhi:continue
            value+=np.sum(y[xlo:xhi,ylo:yhi]*y[xlo+dx:xhi+dx,ylo+dy:yhi+dy])
            included+=(xhi-xlo)*(yhi-ylo)
    return float(value),int(included)


def spatial_estimate(x,W,shape,coarse_blocks):
    n=x[:,:len(W)].sum(0)
    if np.einsum('abc,a,b,c->',W,n,n,n)<=0:
        return dict(status='unsupported_no_distinct_allowed_triples')
    try:
        if coarse_blocks<2:raise ValueError('at least two original spatial blocks required')
        point,gradient=value_gradient(x.sum(0),W)
        variance,pairs=finite_lag_variance(x@gradient,shape)
        if not np.isfinite(variance) or variance<=0:raise ValueError('nonpositive finite-lag variance; no clipping')
        se=float(np.sqrt(variance));critical=float(student_t.ppf(.975,coarse_blocks-1))
        return dict(status='estimated',value=float(point),se=se,
                    spatial_t_ci=[float(point-critical*se),float(point+critical*se)],
                    normal_ci=[float(point-1.959963984540054*se),float(point+1.959963984540054*se)],
                    df=coarse_blocks-1,nominal=.95,exact_coverage=False,
                    variance=variance,included_ordered_cell_pairs=pairs,
                    cell_nm=CELL,dependence_nm=DEPENDENCE,method='finite-lag score covariance; coarse-df t interval')
    except (ValueError,FloatingPointError) as error:
        return dict(status='unsupported',reason=str(error),nominal=.95,exact_coverage=False)


# Bootstrap-t inversion


def invert(spatial,bootstrap):
    if spatial.get('status')!='estimated' or bootstrap.get('status')!='estimated':
        return dict(status='unsupported',reason='both spatial and bootstrap support required',exact_coverage=False)
    try:
        point=float(spatial['value']);other=float(bootstrap['value']);se=float(spatial['se'])
        lo,hi=map(float,bootstrap['pivot_quantiles'])
        if not all(math.isfinite(v) for v in [point,other,se,lo,hi]) or se<=0 or lo>hi:
            raise ValueError('finite point, positive spatial SE and ordered quantiles required')
        if not math.isclose(point,other,rel_tol=1e-10,abs_tol=1e-8):
            raise ValueError('bootstrap and spatial point targets differ')
        interval=[point-hi*se,point-lo*se]
        if not all(map(math.isfinite,interval)):raise ValueError('nonfinite interval')
        return dict(status='estimated',value=point,se=se,spatial_bootstrap_t_ci=interval,
                    pivot_quantiles=[lo,hi],nominal=.95,exact_coverage=False,
                    method='coarse-block bootstrap-t shape with finite-lag spatial covariance scale',
                    original_spatial_se_unchanged=True,truth_fitted_scale=False)
    except (KeyError,TypeError,ValueError,OverflowError) as error:
        return dict(status='unsupported',reason=str(error),exact_coverage=False)


# Nested bootstrap calibration


def calibrated_quantiles(outer_pivots,prepivot_values):
    outer=np.asarray(outer_pivots,float);u=np.asarray(prepivot_values,float)
    if outer.ndim!=1 or outer.shape!=u.shape or len(outer)<2 or not np.isfinite(outer).all() or not np.isfinite(u).all() or (u<=0).any() or (u>=1).any():
        raise ValueError('finite matching outer pivots and internal CDF values required')
    levels=np.quantile(u,[.025,.975]);q=np.quantile(outer,levels)
    return levels,q


def nested_estimate(x,W,spatial,seed):
    if spatial.get('status')!='estimated':return dict(status='unsupported',reason='spatial support required'),{}
    try:
        x=np.asarray(x,float);W=np.asarray(W,float);nb=len(x)
        if nb<2 or not np.isfinite(x).all():raise ValueError('finite block contributions required')
        point,gradient=value_gradient(x.sum(0),W)
        # Match the exact original999 outer block draws, without new random tuning.
        picks=np.random.default_rng(200000+seed).integers(0,nb,(OUTER,nb))
        outer=np.array([np.bincount(p,minlength=nb) for p in picks])
        values,gradients=value_gradient(outer@x,W)
        outer_se=block_se(x,gradients,outer);pivots=(values-point)/outer_se
        rng=np.random.default_rng(500000+seed);u=np.empty(OUTER);digest=hashlib.sha256()
        for lo in range(0,OUTER,BATCH):
            hi=min(OUTER,lo+BATCH);b=hi-lo
            inner=rng.multinomial(nb,outer[lo:hi]/nb,size=(INNER,b)).transpose(1,0,2)
            if not np.all(inner.sum(-1)==nb):raise ValueError('inner record-block multiplicity mismatch')
            digest.update(np.ascontiguousarray(inner,dtype='<i8').tobytes())
            flat=inner.reshape(b*INNER,nb)
            iv,ig=value_gradient(flat@x,W);ise=block_se(x,ig,flat)
            it=((iv.reshape(b,INNER)-values[lo:hi,None])/ise.reshape(b,INNER))
            u[lo:hi]=(.5+(it<=pivots[lo:hi,None]).sum(1))/(INNER+1)
        levels,q=calibrated_quantiles(pivots,u)
        row=invert(spatial,dict(status='estimated',value=float(point),pivot_quantiles=q.tolist()))
        if row['status']!='estimated':return row,{}
        row.update(method='nested coarse-bootstrap calibration with finite-lag spatial covariance scale',
                   outer_draws=OUTER,inner_draws_per_outer=INNER,calibrated_cdf_levels=levels.tolist(),
                   inner_multiplicities_sha256=digest.hexdigest(),simulation_truth_used=False)
        return row,dict(outer_pivots=pivots,prepivot_values=u)
    except (ValueError,FloatingPointError,OverflowError) as error:
        return dict(status='unsupported',reason=str(error),exact_coverage=False),{}


# Input support


def validate_catalogue(xy,labels,K,radius,units):
    xy=np.asarray(xy,float);labels=np.asarray(labels)
    if units!='nm':raise ValueError('Explicit nm coordinates required; no implicit conversion')
    if K not in (4,6) or radius not in (10.,30.):raise ValueError('Only original engineering K/radius values supported')
    if xy.ndim!=2 or xy.shape[1]!=2 or len(xy)!=len(labels) or labels.ndim!=1 or len(xy)==0:
        raise ValueError('Matching nonempty coordinate/label rows required')
    if not np.isfinite(xy).all() or labels.dtype.kind not in 'iu' or (labels<0).any() or (labels>=K).any():
        raise ValueError('Finite coordinates and explicit integer code support required')
    if len(xy)>MAX_RECORDS:raise ValueError('Record cap exceeded; no subsampling fallback')
    counts=np.bincount(labels,minlength=K)
    if (counts<2).any():raise ValueError('Original two-record minimum per declared code required')
    ij=np.floor((xy-xy.min(0))/100.).astype(int);shape=tuple(int(x)+1 for x in ij.max(0))
    coarse=len(np.unique(ij//20,axis=0))
    if shape[0]*shape[1]>MAX_CELLS or coarse>MAX_COARSE:
        raise ValueError('Declared grid/bootstrap memory support exceeded; no crop fallback')
    if min(np.ptp(xy,axis=0))<=500.:raise ValueError('Original footprint-size support failed')
    return xy,labels,dict(records=len(xy),K=K,counts=counts.tolist(),fine_cells=shape[0]*shape[1],
        coarse_blocks=coarse,units='nm',no_records_discarded=True)
