"""Explicit observation-stage reporting for the approximate joint workflow."""
import numpy as np
from .correlation import point_pair_g
from .inference import joint_from_sufficient,triplet_sufficient

STAGES=('localizations','native_positions','processed_map','unknown')


def analyze_joint(xy,labels,input_stage='unknown',seed=0,B=999,radius=30.,n_labels=None):
    """Full-centre joint connected-star summary, with declared observation stage.

    The registered all-domain confirmation failed. Every observation stage
    receives research diagnostics without a calibrated molecular-null decision.
    The historical analyze/CLI remain available for reproduction.
    """
    if input_stage not in STAGES:raise ValueError('unknown observation stage')
    xy=np.asarray(xy,float);labels=np.asarray(labels)
    if xy.ndim!=2 or xy.shape[1]!=2 or labels.shape!=(len(xy),) or not len(xy) or not np.isfinite(xy).all():
        raise ValueError('finite nonempty (N,2) coordinates and matching labels required')
    if not np.issubdtype(labels.dtype,np.number) or not np.isfinite(labels).all() or (labels<0).any() or (labels!=np.floor(labels)).any():
        raise ValueError('nonnegative integer labels required')
    labels=labels.astype(int)
    if n_labels is not None and (int(n_labels)!=n_labels or n_labels<=labels.max()):
        raise ValueError('n_labels must exceed every label')
    K=int(labels.max())+1 if n_labels is None else int(n_labels)
    if int(B)!=B or B<2:raise ValueError('B must be an integer >=2')
    if not np.isfinite(radius) or not 0<radius<=100:raise ValueError('radius must be in (0,100] nm')
    native=input_stage in ['localizations','native_positions']
    out=dict(records=len(xy),label_counts=np.bincount(labels,minlength=K).tolist(),
             input_stage_declared=input_stage,assumptions_verified=False,
             weighting='input records; repeated localizations are not molecule counts',
             molecular_target_qualification=('conditional on native-observation assumptions' if native
                                              else 'unsupported after unknown/label-dependent map selection'),
             method='joint block percentile bootstrap; historical point estimator',
             nominal=.95,coverage_guaranteed=False,radius_nm=float(radius),
             qualification_status='failed registered all-domain confirmation; research diagnostic only',
             molecular_null_decision=None,positive_interval=None,
             propagated=['pair counts and normalization by label counts','pair baseline',
                         'empirical label probabilities','triplet numerator and record denominator'],
             conditional=['local neighbourhoods','shifted-null and plateau normalizations'],
             required_assumptions=['one permanent code per interpreted source; a probe is not automatically one protein',
                                   'applicable independent-source/cluster observation model',
                                   'adequate background and spatial support','no unmodelled source-selection coupling'])
    if len(np.unique(labels))<3:
        return dict(out,status='unsupported',reason='three occupied labels required')
    try:
        with np.errstate(divide='ignore',invalid='ignore'):
            res=point_pair_g(xy,labels,np.random.default_rng(100000+seed),n_labels=K,n_shift=12,fraction=1.)
            out['occupied_blocks']=len(res['n'])
            suff=triplet_sufficient(xy,labels,res,radius=radius)
            fits=joint_from_sufficient(res,suff,np.random.default_rng(200000+seed),B=B,radius=radius)
        fit=fits['all']
        out.update(occupied_blocks=len(res['n']),masks=fits,status=fit['status'])
        if fit['status']=='estimated':
            out.update(value=fit['value'],ci=fit['ci'],numerical_interval_positive=bool(fit['ci'][0]>0))
    except (ValueError,FloatingPointError) as error:
        out.update(status='unsupported',reason=str(error),positive_interval=None)
    return out


def format_joint_report(result):
    lines=[f"records {result['records']}; declared stage {result['input_stage_declared']}",
           f"molecular target: {result['molecular_target_qualification']}"]
    lines.append(result['qualification_status'])
    lines.append(f"occupied labels {sum(c>0 for c in result['label_counts'])}; occupied 2-um blocks "
                 f"{result.get('occupied_blocks','not computed')}; block count alone does not certify coverage")
    if result['status']=='estimated':
        lo,hi=result['ci']
        lines.append(f"joint connected star {result['value']:.6g} [{lo:.6g}, {hi:.6g}]; approximate nominal95% interval")
        if result['positive_interval'] is None:
            lines.append('Molecular-null decision withheld; record-level interval is descriptive.')
        else:
            lines.append('Interval '+('excludes zero' if result['positive_interval'] else 'includes zero')+
                         ' under the specified observation model.')
    else:lines.append('Unsupported: '+result.get('reason',result['status']))
    lines.append('Propagated: '+', '.join(result['propagated']))
    lines.append('Conditional inputs: '+', '.join(result['conditional']))
    lines.append('Assumptions remain unverified; no literal assembly size or geometry is identified.')
    return '\n'.join(lines)
