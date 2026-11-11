"""Fixed prospective comparison schedule and discrete adequacy-rule guards."""
COHORTS=[
    ('adequate_monomer',701000,40,[1.,0.,0.],True),
    ('adequate_monomer_dimer',702000,40,[.5,.5,0.],True),
    ('missing_trimer_5',703000,20,[.75,.2,.05],False),
    ('missing_trimer_10',704000,20,[.7,.2,.1],False),
    ('missing_trimer_20',705000,20,[.6,.2,.2],False),
    ('missing_trimer_30',706000,20,[.5,.2,.3],False)]


def schedule():
    return [dict(key=str(base+i),seed=base+i,index=i,cohort=name,fractions=f,adequate_library=adequate)
            for i in range(1,41) for name,base,reps,f,adequate in COHORTS if i<=reps]


def rank_p(observed,bootstrap):
    import math
    if len(bootstrap)!=19 or not all(math.isfinite(v) for v in [observed,*bootstrap]):
        raise ValueError('Exactly19 complete finite refitted-null scores required')
    return (1+sum(v>=observed for v in bootstrap))/20


def bootstrap_seed(field_seed,index,simulation=False):
    if not 0<=index<19:raise ValueError('Declared19 refits only')
    return (300000000 if simulation else 200000000)+field_seed*100+index
