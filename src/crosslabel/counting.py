"""Closed-radius ordered pairs, including distinct coincident records.

The historical correlation implementation remains available unchanged. This
adapter fixes only the zero-distance counting boundary; it does not certify
any spatial estimator, background model or confidence interval.
"""
import numpy as np
from . import correlation as C


def block_pair_counts_closed(xy, label, edges, block=None, n_labels=None):
    """Count i!=j by anchor block/labels in closed outer-radius shells.

    The first shell is [0, edges[1]]; subsequent shells are
    (edges[k], edges[k+1]]. Every record remains eligible. Coincident records
    contribute, but a record never pairs with itself. No shifted-null option.
    """
    xy=np.asarray(xy,float);label=np.asarray(label);edges=np.asarray(edges,float)
    if xy.ndim!=2 or xy.shape[1]!=2 or not len(xy) or not np.isfinite(xy).all():
        raise ValueError('nonempty finite Nx2 coordinates required')
    if label.shape!=(len(xy),) or not np.issubdtype(label.dtype,np.integer) or (label<0).any():
        raise ValueError('matching nonnegative integer labels required')
    K=int(label.max())+1 if n_labels is None else n_labels
    if not isinstance(K,(int,np.integer)) or K<1 or (label>=K).any():
        raise ValueError('positive label support must contain every record label')
    if edges.ndim!=1 or len(edges)<2 or edges[0]!=0 or not np.isfinite(edges).all() or (np.diff(edges)<=0).any():
        raise ValueError('strictly increasing finite radius edges starting at zero required')
    block=C.blocks(xy) if block is None else np.asarray(block)
    if block.shape!=(len(xy),) or not np.issubdtype(block.dtype,np.integer) or (block<0).any():
        raise ValueError('matching nonnegative integer anchor blocks required')
    if block.max()>np.iinfo(np.intp).max:
        raise ValueError('anchor block index exceeds platform support')
    block=block.astype(np.intp,copy=False)
    N,n=C.block_pair_counts(xy,label,edges=edges,block=block,n_labels=K)
    # Historical cumulative differencing drops all distances exactly equal
    # to edges[0]=0. Recover those ordered pairs by exact coordinate groups.
    _,group=np.unique(xy,axis=0,return_inverse=True)
    counts=np.zeros((int(group.max())+1,K),dtype=np.int64)
    np.add.at(counts,(group,label),1)
    for a in range(K):
        select=label==a
        for b in range(K):
            weight=counts[group[select],b]-(a==b)
            N[:,a,b,0]+=np.bincount(block[select],weights=weight,minlength=len(n))
    return N,n
