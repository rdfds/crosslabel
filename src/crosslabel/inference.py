"""Joint nuisance propagation for connected stars; see qualification records.

The historical point estimator is retained. This is an approximate conditional
block bootstrap, not a correction of label-dependent merging or an exact
confidence set. Known shifted-null/plateau inputs remain conditional. No
automatic molecular/assembly interpretation is supplied.
"""
import numpy as np
from .correlation import EDGES, baseline_excess
from .labels import lens_integral, star_counts


def triplet_sufficient(xy, label, res, radius=30., annulus=(100.,250.)):
    """Full-centre block statistics, including the separately named label masks.

    Input records are not automatically molecules. Localizations produce a
    record-weighted statistic; a processed map does not establish independent
    source selection. The non-consecutive mask is not automatically qualified.
    """
    xy=np.asarray(xy,float);label=np.asarray(label);block=np.asarray(res['block'])
    if xy.ndim!=2 or xy.shape[1]!=2 or label.shape!=(len(xy),) or block.shape!=(len(xy),):
        raise ValueError('coordinates, labels and centre blocks must describe the same records')
    if not len(xy) or not np.isfinite(xy).all() or not np.isfinite(label).all() or (label<0).any() or (label!=np.floor(label)).any():
        raise ValueError('finite coordinates and nonnegative integer labels required')
    if not np.isfinite([radius,*annulus]).all() or radius<=0 or annulus[0]<radius or annulus[1]<=annulus[0]:
        raise ValueError('positive radius and disjoint ordered background annulus required')
    K, nb = res['n'].shape[1], res['n'].shape[0]
    if (label>=K).any() or not np.isfinite(block).all() or (block<0).any() or (block>=nb).any() or (block!=np.floor(block)).any():
        raise ValueError('labels or block IDs outside the supplied count support')
    label=label.astype(int);centre_block=block.astype(int)
    counts=np.zeros((nb,K));np.add.at(counts,(centre_block,label),1)
    if not np.array_equal(counts,res['n']):raise ValueError('coordinate block/label counts do not match pair inputs')
    sep=np.abs(np.subtract.outer(np.arange(K),np.arange(K)))
    local=star_counts(xy,label,annulus[1],K)-star_counts(xy,label,annulus[0],K)
    lam=local/(np.pi*(annulus[1]**2-annulus[0]**2))
    e=star_counts(xy,label,radius,K)-lam*np.pi*radius**2
    ans={}
    for name,M in {'all':sep>=1,'non_consecutive':sep>=2}.items():
        W=(M[:,:,None]&M[:,None,:]&M[None,:,:]).astype(float)
        per=np.einsum('ib,ibc,ic->i',e,W[label],e)
        A=np.bincount(centre_block,weights=per,minlength=nb)
        Z=np.zeros((nb,K,K))
        for block in range(nb):
            ix=centre_block==block
            Z[block]=np.einsum('ib,ic,ibc->bc',lam[ix],lam[ix],W[label[ix]])
        ans[name]=dict(A=A,Z=Z,W=W)
    return ans


def _refit(N,n,A,Z,W,nd,pr,radius):
    counts=n.sum(-1);pi=n/counts[...,None]
    D=np.einsum('...a,...b,...c,abc->...',pi,pi,pi,W)
    den=n[..., :,None]*n[...,None,:]
    g=N/np.maximum(den[...,None],1)/(np.pi*np.diff(EDGES**2)*nd*pr)
    L=lens_integral(EDGES,baseline_excess(g),radius)
    return (A-(Z*L).sum((-2,-1)))/(2*D*counts)


def joint_from_sufficient(res,sufficient,rng,B=999,radius=30.):
    """Refit pair/label nuisances jointly in every block-bootstrap draw.

    No invalid bootstrap draws are discarded. All masks and unsupported
    outcomes are returned. This propagates N,n,A,Z covariance while conditioning
    on local neighbourhoods and shifted-null/plateau normalizations. It does not
    remove finite-background or source-selection bias.
    """
    if int(B)!=B or B<2:raise ValueError('B must be an integer >=2')
    if not np.isfinite(radius) or radius<=0:raise ValueError('radius must be positive')
    N,n=np.asarray(res['N']),np.asarray(res['n'])
    if n.ndim!=2 or N.shape!=(len(n),n.shape[1],n.shape[1],len(EDGES)-1):
        raise ValueError('incompatible block count arrays')
    if not len(n) or not np.isfinite(n).all() or not np.isfinite(N).all() or (n<0).any() or (N<0).any() or n.sum()<=0:
        raise ValueError('finite nonnegative counts and positive record support required')
    nd,pr=res['null_density'],res['plateau_ratio']
    if not np.isfinite([nd,pr]).all() or nd<=0 or pr<=0:raise ValueError('positive finite normalization required')
    nb=len(n);pick=rng.integers(0,nb,(int(B),nb))
    mult=np.array([np.bincount(p,minlength=nb) for p in pick])
    nr=mult@n;Nr=np.tensordot(mult,N,axes=(1,0));ng=n.sum(0);Ng=N.sum(0)
    pi=ng/ng.sum();ans={}
    for mask,s in sufficient.items():
        A,Z,W=s['A'],s['Z'],s['W']
        D=np.einsum('a,b,c,abc',pi,pi,pi,W)
        if D<=0:
            ans[mask]=dict(status='unsupported_no_distinct_allowed_triples');continue
        Ar=mult@A;Zr=np.tensordot(mult,Z,axes=(1,0))
        point=_refit(Ng,ng,A.sum(),Z.sum(0),W,nd,pr,radius)
        draws=_refit(Nr,nr,Ar,Zr,W,nd,pr,radius)
        if not np.isfinite(point) or not np.isfinite(draws).all():
            raise ValueError('nonfinite joint point/bootstrap; no draws discarded')
        ans[mask]=dict(status='estimated',value=float(point),ci=np.percentile(draws,[2.5,97.5]).tolist(),
                       se=float(draws.std(ddof=1)),nominal=.95,exact_coverage=False)
    return ans
