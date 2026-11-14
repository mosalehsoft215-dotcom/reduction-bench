import numpy as np
from space import construct
R=.00831446261815324
TR=450.

def subset(s,c,o):
    # Desorption decreases N_C-N_O; oxidation preserves it; hopping preserves
    # both counts. States above its initial value have exactly zero probability.
    keep=np.flatnonzero((s['obs'][:,0]<=c)&(s['obs'][:,1]<=o)&(s['obs'][:,0]-s['obs'][:,1]<=c-o))
    remap=np.full(len(s['reps']),-1);remap[keep]=np.arange(len(keep))
    edge=(remap[s['src']]>=0)&(remap[s['dst']]>=0)
    return dict(reps=s['reps'][keep],deg=s['deg'][keep],obs=s['obs'][keep],
                src=remap[s['src'][edge]],dst=remap[s['dst'][edge]],
                labels=s['labels'][edge],types=s['types'])


def rates(x,T,types):
    z=types[:,0];v=np.zeros(len(z),dtype=np.result_type(x,float));d=1/R*(1/TR-1/T)
    a=z==0;v[a]=10**x[1]*np.exp(x[0]*d-(x[2]*types[a,1]+x[3]*types[a,2])/(R*T))
    a=z==1;v[a]=10**x[5]*np.exp(x[4]*d)
    a=z>=2;v[a]=10**x[7]*np.exp(x[6]*d+(x[2]*types[a,1]+x[3]*types[a,2]+x[8]*types[a,3])/(2*R*T))
    a=z==3;v[a]*=.2*np.exp((x[9]-x[6])*d)
    return v


def rate_jac(x,T,types):
    v=rates(x,T,types);out=np.zeros((len(types),10));z=types[:,0];d=(1/TR-1/T)/R
    a=z==0;out[a,0]=v[a]*d;out[a,1]=v[a]*np.log(10)
    out[a,2]=-v[a]*types[a,1]/(R*T);out[a,3]=-v[a]*types[a,2]/(R*T)
    a=z==1;out[a,4]=v[a]*d;out[a,5]=v[a]*np.log(10)
    a=z==2;out[a,6]=v[a]*d;out[a,7]=v[a]*np.log(10)
    a=z==3;out[a,7]=v[a]*np.log(10);out[a,9]=v[a]*d
    a=z>=2;out[a,8]=v[a]*types[a,3]/(2*R*T)
    out[a,2]=v[a]*types[a,1]/(2*R*T);out[a,3]=v[a]*types[a,2]/(2*R*T)
    return out


def initial(s,prep):
    if isinstance(prep,list):
        # Find the orbit of a specified microscopic pattern by direct symmetry.
        from space import PERM,POW
        code=min(sum(prep[j]*int(POW[p[j]]) for j in range(12)) for p in PERM)
        p=np.zeros(len(s['reps']));p[np.flatnonzero(s['reps']==code)[0]]=1.;return p
    c,o=prep['counts'];p=s['deg'].astype(float)*((s['obs'][:,0]==c)&(s['obs'][:,1]==o));return p/p.sum()


