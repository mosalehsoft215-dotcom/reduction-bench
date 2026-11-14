"""Exact symmetry quotient of a small periodic catalytic lattice."""
import numpy as np
from numba import njit

ROWS,COLS=3,4
N=ROWS*COLS
POW=3**np.arange(N,dtype=np.int64)
BONDS=np.array(sorted({tuple(sorted((r*COLS+c,rr*COLS+cc)))
    for r in range(ROWS) for c in range(COLS)
    for rr,cc in [(r,(c+1)%COLS),((r+1)%ROWS,c)]}),dtype=np.int64)
PERM=np.array([[((sr*r+dr)%ROWS)*COLS+(sc*c+dc)%COLS
    for r in range(ROWS) for c in range(COLS)]
    for dr in range(ROWS) for dc in range(COLS) for sr in [1,-1] for sc in [1,-1]],dtype=np.int64)

@njit
def digits(code):
    a=np.zeros(N,np.int64)
    for j in range(N):a[j]=code%3;code//=3
    return a

@njit
def canonical_map():
    lookup=np.full(3**N,-1,np.int64);deg=np.zeros(3**N,np.int64)
    for code in range(3**N):
        if lookup[code]>=0:continue
        a=digits(code);c=0;o=0
        for v in a:
            if v==1:c+=1
            if v==2:o+=1
        if c>6 or o>4:continue
        for perm in PERM:
            image=0
            for j in range(N):image+=a[j]*POW[perm[j]]
            if lookup[image]<0:deg[code]+=1
            lookup[image]=code
    return lookup,deg

@njit
def pair_counts(a):
    out=np.zeros(3,np.int64)
    for i,j in BONDS:
        if a[i]==1 and a[j]==1:out[0]+=1
        if a[i]==2 and a[j]==2:out[2]+=1
        if (a[i]==1 and a[j]==2) or (a[i]==2 and a[j]==1):out[1]+=1
    return out

@njit
def raw_edges(reps,mapping):
    capacity=len(reps)*80
    source=np.zeros(capacity,np.int64);dest=np.zeros(capacity,np.int64)
    labels=np.zeros((capacity,4),np.int64);count=0
    obs=np.zeros((len(reps),5),np.int64)
    for idx in range(len(reps)):
        code=reps[idx];a=digits(code);bc=pair_counts(a)
        obs[idx,0]=np.sum(a==1);obs[idx,1]=np.sum(a==2);obs[idx,2:]=bc
        for i in range(N):
            if a[i]!=1:continue
            nc=0;no=0
            for b0,b1 in BONDS:
                j=-1
                if b0==i:j=b1
                if b1==i:j=b0
                if j>=0:
                    nc+=int(a[j]==1);no+=int(a[j]==2)
            source[count]=idx;dest[count]=mapping[code-POW[i]]
            labels[count,0]=0;labels[count,1]=nc;labels[count,2]=no;count+=1
        for i,j in BONDS:
            if (a[i]==1 and a[j]==2) or (a[i]==2 and a[j]==1):
                source[count]=idx;dest[count]=mapping[code-a[i]*POW[i]-a[j]*POW[j]]
                labels[count,0]=1;count+=1
            if a[i]==0 and a[j]>0:site=j;empty=i
            elif a[j]==0 and a[i]>0:site=i;empty=j
            else:continue
            species=a[site];a[empty]=species;a[site]=0
            after=pair_counts(a);a[site]=species;a[empty]=0
            target=code+species*(POW[empty]-POW[site])
            source[count]=idx;dest[count]=mapping[target]
            labels[count,0]=1+species;labels[count,1:]=after-bc;count+=1
    return source[:count],dest[:count],labels[:count],obs

def construct():
    lookup,deg=canonical_map();reps=np.flatnonzero(deg)
    idx=np.full(3**N,-1,np.int64);idx[reps]=np.arange(len(reps))
    mapping=np.full(3**N,-1,np.int64);valid=lookup>=0;mapping[valid]=idx[lookup[valid]]
    src,dst,labels,obs=raw_edges(reps,mapping)
    assert dst.min()>=0
    types,labels=np.unique(labels,axis=0,return_inverse=True)
    # Quotient self-jumps remain recorded: they contribute to hopping counts
    # even though their probability-generator contribution cancels exactly.
    return dict(reps=reps,mapping=mapping,deg=deg[reps],src=src,dst=dst,labels=labels,types=types,obs=obs)

