"""Independent compound-Poisson generator; fixed new schedule, no outcome tuning."""
import numpy as np

FAMILIES={
 'reference':dict(shape='triangle',spacing=12.,p=.47,density=70.,K=4,mask=False),
 'linear':dict(shape='line',spacing=12.,p=.47,density=70.,K=4,mask=False),
 'wide':dict(shape='triangle',spacing=24.,p=.47,density=70.,K=4,mask=False),
 'low_retention':dict(shape='triangle',spacing=12.,p=.25,density=70.,K=4,mask=False),
 'low_density_six_codes':dict(shape='triangle',spacing=12.,p=.47,density=35.,K=6,mask=False),
 'known_intensity':dict(shape='triangle',spacing=12.,p=.47,density=70.,K=4,mask=True),
}
SIDE=10000.;SIGMA=float(np.sqrt(25.25));RADIUS=30.;REPS=100;PAIRED=12
ADEQUATE=(.5,.5,0.);ALTERNATIVE=(.7,.2,.1)

def schedule():
    return [dict(key=str(801000+1000*j+500*alt+i),seed=801000+1000*j+500*alt+i,
        family=f,index=i,adequate=not alt,paired=i<=PAIRED)
        for i in range(1,REPS+1) for j,f in enumerate(FAMILIES) for alt in [0,1]]

def mask_pmf():
    i=np.arange(20)+.5
    a=1.+.5*np.cos(2*np.pi*i/20)
    return np.repeat(a[None,:],20,axis=0)

def shape(m,spacing,kind):
    if m==1:return np.zeros((1,2))
    if m==2 or kind=='line':return np.column_stack((spacing*(np.arange(m)-(m-1)/2),np.zeros(m)))
    a=2*np.pi*np.arange(m)/m;r=spacing/(2*np.sin(np.pi/m))
    return r*np.column_stack((np.cos(a),np.sin(a)))

def component_contribution(x,keep,radius=RADIUS):
    """Count observed own-component anchored pairs, independent of marks."""
    distance2=np.sum((x[:,:,None,:]-x[:,None,:,:])**2,axis=-1)
    n=((distance2<=radius**2)&keep[:,None,:]).sum(-1)-keep
    return int(np.sum((n*(n-1)//2)*keep))

def generate(plan):
    c=FAMILIES[plan['family']];rng=np.random.default_rng(plan['seed']);f=ADEQUATE if plan['adequate'] else ALTERNATIVE
    xy=[];labels=[];n_centers=[];own=0
    # Independent Poisson counts for each parent type; f is a SUBUNIT fraction.
    for m,fraction in enumerate(f,1):
        n=int(rng.poisson(c['density']*(SIDE/1000)**2*fraction/(m*c['p'])));n_centers.append(n)
        if not n:continue
        if c['mask']:
            w=mask_pmf().ravel();b=rng.choice(400,size=n,p=w/w.sum())
            centres=500.*(np.column_stack((b%20,b//20))+rng.random((n,2)))
        else:centres=SIDE*rng.random((n,2))
        ang=2*np.pi*rng.random(n);cos=np.cos(ang)[:,None];sin=np.sin(ang)[:,None];t=shape(m,c['spacing'],c['shape'])
        x=centres[:,None,:]+np.stack((cos*t[:,0]-sin*t[:,1],sin*t[:,0]+cos*t[:,1]),axis=2)
        keep=rng.random((n,m))<c['p'];x+=rng.normal(0.,SIGMA,x.shape)
        label=rng.integers(c['K'],size=(n,m),dtype=np.int32)
        own+=component_contribution(x,keep);xy.append(x[keep]);labels.append(label[keep])
    x=np.concatenate(xy);l=np.concatenate(labels)
    if not len(x):raise ValueError('Empty fixed field')
    return x,l,dict(m2=own/len(x),own_anchored_pairs=own,parent_counts=n_centers,records=len(x))

def rank_p(v,refits):
    if len(refits)!=19 or not np.isfinite([v,*refits]).all():raise ValueError('Exactly19 finite null refits required')
    return (1+sum(x>=v for x in refits))/20
