"""Official SPINNA and an explicitly separate full-refit adequacy wrapper."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from picasso import spinna
from spinna_library_design import rank_p,bootstrap_seed



def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')


def mixer():
    structures=[spinna.Structure('monomer').define_coordinates('CD20',[0.],[0.]),
                spinna.Structure('dimer').define_coordinates('CD20',[-6.,6.],[0.,0.])]
    return spinna.StructureMixer(structures,label_unc={'CD20':float(np.sqrt(25.25))},le={'CD20':.47},
        width=10000.,height=10000.,random_rot_mode='2D',nn_counts={'CD20-CD20':2})


def fit(xy,seed,out,native_bootstrap):
    np.random.seed(seed);m=mixer();space=spinna.generate_N_structures(m.structures,{'CD20':int(len(xy)/.47)},21)
    f=spinna.SPINNA(m,{'CD20':xy},N_sim=10)
    if f.dists_gt[0].shape!=(len(xy),2):raise ValueError('Declared two-neighbour score changed')
    original=f.NN_scorer;calls=[]
    def observe(counts,callback=None):
        start=time.time();n,s=original(counts,callback=callback);i=int(np.argmin(s))
        calls.append(dict(counts=n.tolist(),scores=s.tolist(),best_counts=n[i].tolist(),
            best_props=m.convert_counts_to_props(n[i]).tolist(),best_score=float(s[i]),seconds=time.time()-start))
        write(out/'CALLS.json',calls);return n,s
    f.NN_scorer=observe;start=time.time()
    props,score,scores=f.fit(space,fitting_mode='brute-force',asynch=False,bootstrap=native_bootstrap,return_scores=True)
    if len(calls)!=(21 if native_bootstrap else 1):raise ValueError('Official refit count differs')
    if native_bootstrap:
        best_props=np.asarray(props[0]);std=np.asarray(props[1]);best_score=float(score[0]);score_std=float(score[1])
        bp=np.array([c['best_props'] for c in calls[1:]],dtype=best_props.dtype);bs=np.array([c['best_score'] for c in calls[1:]])
        if not np.array_equal(bp.std(0),std) or bs.std()!=score_std:raise ValueError('Native bootstrap uncertainty replay differs')
        unc=dict(proportion_std=std.tolist(),KS2_std=score_std,bootstrap_props=bp.tolist(),bootstrap_scores=bs.tolist())
    else:best_props=np.asarray(props);best_score=float(score);unc=None
    result=dict(proportions=best_props.tolist(),best_score=best_score,best_counts=calls[0]['best_counts'],
        uncertainty=unc,seed=seed,seconds=time.time()-start,candidate_evaluations=sum(len(c['scores']) for c in calls),
        official_scorer_calls=len(calls),native_bootstrap_refits=20 if native_bootstrap else 0,
        N_sim=10,NND_orders=2,asynch=False,full_grid_refit=True)
    write(out/'RESULT.json',result);return result


def main(input_path,out,seed):
    if out.exists():raise FileExistsError('Fresh comparator output required')
    if sha(Path(spinna.__file__))!='f097830230cc13fa18993a17abf4335d90d11768ffb63c5cdf1169a4e2190c79':raise ValueError('Official source changed')
    with np.load(input_path,allow_pickle=False) as d:xy=d['xy']
    if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all():raise ValueError('Finite shared input coordinates required')
    out.mkdir();start=time.time();state=dict(state='running',pid=os.getpid(),started_epoch=start,
        completed_parametric_refits=0,native_fit_complete=False,declared_parametric_refits=19)
    write(out/'STATE.json',state)
    try:
        native=out/'native';native.mkdir();original=fit(xy,1000000+seed,native,True)
        state['native_fit_complete']=True;write(out/'STATE.json',state);refits=[]
        for i in range(19):
            p=out/f'null-refit-{i:02d}';p.mkdir()
            np.random.seed(bootstrap_seed(seed,i,True));m=mixer()
            boot=m.run_simulation(np.asarray(original['best_counts'],dtype=int))['CD20']
            np.savez_compressed(p/'input.npz',xy=boot)
            value=fit(boot,bootstrap_seed(seed,i,False),p,False)
            refits.append(dict(index=i,simulation_seed=bootstrap_seed(seed,i,True),fit_seed=bootstrap_seed(seed,i,False),
                records=len(boot),input_sha256=sha(p/'input.npz'),best_score=value['best_score'],
                proportions=value['proportions'],candidate_evaluations=value['candidate_evaluations'],seconds=value['seconds']))
            write(out/'PARAMETRIC_REFITS.json',refits);state['completed_parametric_refits']+=1;write(out/'STATE.json',state)
        pvalue=rank_p(original['best_score'],[r['best_score'] for r in refits])
        result=dict(native=original,parametric_refits=refits,model_bootstrap_p=pvalue,
            decision='reject_restricted_library' if pvalue<=.05 else 'not_rejected',
            library='monomer/12nm dimer',alpha=.05,bootstrap_p_is_approximate=True,
            native_uncertainty_and_custom_adequacy_distinct=True,coordinates_sha256=sha(input_path),
            candidate_evaluations=original['candidate_evaluations']+sum(r['candidate_evaluations'] for r in refits),
            official_fitting_calls=20,native_bootstrap_refits=20,extra_full_null_refits=19,seconds=time.time()-start)
        write(out/'RESULT.json',result);state['state']='complete'
    except BaseException as e:state.update(state='failed_or_partial',error=repr(e));raise
    finally:state.update(ended_epoch=time.time(),seconds=time.time()-start);write(out/'STATE.json',state)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,required=True);a=p.parse_args();main(a.input,a.output,a.seed)
