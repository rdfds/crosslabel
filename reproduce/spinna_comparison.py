"""Official M/D and M/D/T fit plus explicit refitted-null model comparison."""
import argparse,hashlib,json,os,time
from pathlib import Path
from contextlib import contextmanager
import numpy as np
from picasso import spinna
from spinna_design import FAMILIES,shape,mask_pmf,rank_p,SIDE,SIGMA

def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clock():return time.time(),time.monotonic(),time.process_time()
def elapsed(a):return dict(zip(['wall','monotonic','cpu'],(b-c for b,c in zip(clock(),a))))

@contextmanager
def streams(seed):
    """Supply seeds to otherwise unseeded native mask draws; restore on exit."""
    np.random.seed(seed);original=np.random.default_rng;ss=np.random.SeedSequence([seed,887]);count=[0]
    def bound(given=None):
        if given is not None:return original(given)
        count[0]+=1;return original(ss.spawn(1)[0])
    np.random.default_rng=bound
    try:yield count
    finally:np.random.default_rng=original

def mixer(c,expanded):
    structures=[]
    for m in range(1,4 if expanded else 3):
        t=shape(m,c['spacing'],c['shape']);structures.append(spinna.Structure(str(m)).define_coordinates('CD20',t[:,0],t[:,1]))
    extra={}
    if c['mask']:
        extra=dict(mask_dict={'mask':{'ALL':mask_pmf()},'info':{'ALL':{
            'Camera pixelsize (nm)':1.,'Binsize (nm)':[500.,500.],
            'x_min':0.,'y_min':0.,'x_max':SIDE,'y_max':SIDE}}})
    return spinna.StructureMixer(structures,label_unc={'CD20':SIGMA},le={'CD20':c['p']},width=SIDE,height=SIDE,
        random_rot_mode='2D',nn_counts={'CD20-CD20':2},**extra)

def fit(xy,c,expanded,native,seed,out):
    out.mkdir();began=clock();m=mixer(c,expanded);count=int(len(xy)/c['p'])
    grid=m.convert_N_structures_to_array(spinna.generate_N_structures(m.structures,{'CD20':count},21))
    # Picasso starts with up to 21 candidates per free dimension and
    # removes infeasible rows before returning integer structure counts.
    # In particular, the expanded grid size depends on the observed count.
    columns = 3 if expanded else 2
    if (grid.ndim != 2 or grid.shape[1] != columns
            or not 2 <= len(grid) <= 21 ** (columns - 1)
            or grid.dtype.kind not in 'iu' or np.any(grid < 0)
            or np.any(~np.any(grid > 0, axis=1))):
        raise ValueError(f'Invalid official candidate grid: {grid.shape}')
    f=spinna.SPINNA(m,{'CD20':xy},N_sim=10)
    if f.dists_gt[0].shape!=(len(xy),2):raise ValueError('Fixed two neighbor orders required')
    calls=[];original=f.NN_scorer
    def observe(counts,callback=None):
        t=clock();n,s=original(counts,callback=callback);i=int(np.argmin(s))
        calls.append(dict(counts=n.tolist(),scores=s.tolist(),best_counts=n[i].tolist(),best_score=float(s[i]),
            best_props=m.convert_counts_to_props(n[i]).tolist(),spans=elapsed(t)));write(out/'CALLS.json',calls);return n,s
    f.NN_scorer=observe
    with streams(seed) as seeded:
        props,score=f.fit(grid,fitting_mode='coarse-to-fine' if expanded else 'brute-force',asynch=False,bootstrap=native)
    initial=2 if expanded else 1
    if len(calls)!=initial+(20 if native else 0):raise ValueError('Unexpected fit/refit calls')
    if native:
        p=np.asarray(props[0]);s=float(score[0]);sd=np.asarray(props[1]);ssd=float(score[1])
        bp=np.asarray([r['best_props'] for r in calls[initial:]],dtype=p.dtype)
        bs=np.asarray([r['best_score'] for r in calls[initial:]])
        if not np.array_equal(bp.std(0),sd) or bs.std()!=ssd:raise ValueError('Native uncertainty capture differs')
        uncertainty=dict(prop_std=sd.tolist(),score_std=ssd,dtype=str(p.dtype))
    else:p=np.asarray(props);s=float(score);uncertainty=None
    if s!=calls[initial-1]['best_score']:raise ValueError('Returned native fit differs')
    row=dict(proportions=p.tolist(),score=s,best_counts=calls[initial-1]['best_counts'],uncertainty=uncertainty,
        seed=seed,expanded=expanded,native_bootstrap=native,full_grid_size=len(grid),scorer_calls=len(calls),
        score_entries=sum(len(r['scores']) for r in calls),seeded_unseeded_rng_calls=seeded[0],spans=elapsed(began))
    write(out/'RESULT.json',row);return row

def run(path,family,seed,out):
    if sha(Path(spinna.__file__))!='f097830230cc13fa18993a17abf4335d90d11768ffb63c5cdf1169a4e2190c79':raise ValueError('Official source differs')
    c=FAMILIES[family];out.mkdir();begin=clock();state=dict(state='running',pid=os.getpid(),completed_refits=0,new_independent_fields=0)
    write(out/'STATE.json',state)
    try:
        with np.load(path,allow_pickle=False) as d:xy=d['xy']
        if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all():raise ValueError('Finite input XY required')
        # Disjoint, declared RNG seeds. Neither class nor realized truth is read.
        base=100000000+100*seed
        a=fit(xy,c,False,True,base,out/'MD');b=fit(xy,c,True,True,base+1,out/'MDT')
        samples=[]
        for i in range(19):
            p=out/f'null-{i:02d}';p.mkdir();s=base+2+3*i
            with streams(s) as seeded:boot=mixer(c,False).run_simulation(np.asarray(a['best_counts'],int))['CD20']
            np.savez_compressed(p/'input.npz',xy=boot)
            v=fit(boot,c,False,False,s+1,p/'MD');w=fit(boot,c,True,False,s+2,p/'MDT')
            samples.append(dict(index=i,sim_seed=s,coordinates_sha256=sha(p/'input.npz'),records=len(boot),
                seeded_sim_mask_rng_calls=seeded[0],MD=v,MDT=w,improvement=v['score']-w['score']))
            state['completed_refits']+=1;write(out/'STATE.json',state);write(out/'NULL_REFITS.json',samples)
        dp=a['score']-b['score'];p_abs=rank_p(a['score'],[r['MD']['score'] for r in samples]);p_delta=rank_p(dp,[r['improvement'] for r in samples])
        fits=[a,b]+[v for r in samples for v in [r['MD'],r['MDT']]]
        result=dict(MD=a,MDT=b,improvement=dp,approximate_omnibus_p=p_abs,approximate_improvement_p=p_delta,
            reject_omnibus=p_abs<=.05,reject_improvement=p_delta<=.05,input_sha256=sha(path),family=family,
            native_bootstrap_refits=40,extra_full_procedure_refits=38,top_level_fits=len(fits),
            observed_scorer_calls=sum(r['scorer_calls'] for r in fits),candidate_entries=sum(r['score_entries'] for r in fits),
            known_generator_nuisances_supplied=True,spans=elapsed(begin),official_source_unchanged=True,
            mask_rng_seed_hook=c['mask'],both_p_values_are_custom_approximate_wrappers=True)
        write(out/'RESULT.json',result);state['state']='complete'
    except BaseException as e:state.update(state='failed_or_partial',error=repr(e));raise
    finally:state['spans']=elapsed(begin);write(out/'STATE.json',state)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--family',choices=list(FAMILIES),required=True);p.add_argument('--seed',type=int,required=True);a=p.parse_args();run(a.input,a.family,a.seed,a.output)
