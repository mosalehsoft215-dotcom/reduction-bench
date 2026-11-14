"""Independent 24-element quotient and event-table evaluator.

No reference-solver imports. Translations and one row reflection preserve all public observables.
The reference uses both rectangular reflections and exact state reachability.
"""
import numpy as np
from numba import njit

KEYS=('E_des','log10_k_des_450','epsilon_CC','epsilon_CO','E_ox','log10_k_ox_450','E_hop','log10_k_hop_450','epsilon_OO','E_hop_O')
R=.00831446261815324
POWER=3**np.arange(12,dtype=np.int64)
NEIGH=np.array([[4*r+(c+1)%4,4*r+(c-1)%4,4*((r+1)%3)+c,4*((r-1)%3)+c] for r in range(3) for c in range(4)],np.int64)

@njit
def count_pairs(a):
    b=np.zeros(3,np.int64)
    for i in range(12):
        for j in NEIGH[i]:
            if j<=i:continue
            if a[i]==1 and a[j]==1:b[0]+=1
            elif a[i]==2 and a[j]==2:b[2]+=1
            elif a[i]+a[j]==3 and a[i]*a[j]==2:b[1]+=1
    return b

@njit
def edge_table(a,codes,lookup):
    n=len(codes);ii=np.empty(n*60,np.int64);jj=ii.copy();feat=np.zeros((n*60,4),np.int64);obs=np.zeros((n,5),np.int64);k=0
    for state in range(n):
        row=a[state].copy();before=count_pairs(row);obs[state,0]=np.sum(row==1);obs[state,1]=np.sum(row==2);obs[state,2:]=before
        for i in range(12):
            if row[i]==1:
                ii[k]=state;jj[k]=lookup[codes[state]-POWER[i]];feat[k,0]=0
                for j in NEIGH[i]:feat[k,1]+=int(row[j]==1);feat[k,2]+=int(row[j]==2)
                k+=1
                for j in NEIGH[i]:
                    if row[j]==2:
                        ii[k]=state;jj[k]=lookup[codes[state]-POWER[i]-2*POWER[j]];feat[k,0]=1;k+=1
            if row[i]!=0:
                for j in NEIGH[i]:
                    if row[j]!=0:continue
                    species=row[i];row[i]=0;row[j]=species
                    delta=count_pairs(row)-before;row[i]=species;row[j]=0
                    ii[k]=state;jj[k]=lookup[codes[state]+species*(POWER[j]-POWER[i])];feat[k,0]=species+1;feat[k,1:]=delta;k+=1
    return ii[:k],jj[:k],feat[:k],obs

def build(c=6,o=4,translations=True):
    codes=np.arange(3**12,dtype=np.int64);a=((codes[:,None]//POWER)%3).astype(np.int8)
    keep=(np.sum(a==1,axis=1)<=c)&(np.sum(a==2,axis=1)<=o);codes=codes[keep];a=a[keep]
    canonical=codes.copy()
    if translations:
        for dr in range(3):
            for dc in range(4):
                perm=[4*((r+dr)%3)+(cc+dc)%4 for r in range(3) for cc in range(4)]
                canonical=np.minimum(canonical,a@POWER[perm])
        # One additional row reflection gives a 24-element quotient. The
        # independently implemented reference uses both rectangular reflections.
        for dr in range(3):
            for dc in range(4):
                perm=[4*((-r+dr)%3)+(cc+dc)%4 for r in range(3) for cc in range(4)]
                canonical=np.minimum(canonical,a@POWER[perm])
    reps,inverse,deg=np.unique(canonical,return_inverse=True,return_counts=True)
    lookup=np.full(3**12,-1,np.int64);lookup[codes]=inverse
    ra=((reps[:,None]//POWER)%3).astype(np.int8);i,j,feat,obs=edge_table(ra,reps,lookup)
    assert j.min()>=0
    return Space(reps,deg,lookup,i,j,feat,obs)

class Space:
    def __init__(self,reps,deg,lookup,i,j,feat,obs):
        self.reps=reps;self.deg=deg;self.lookup=lookup;self.i=i;self.j=j;self.feat=feat;self.obs=obs;self.n=len(reps)
        self.hop=feat[:,0]>=2;self.des=feat[:,0]==0;self.ox=feat[:,0]==1
        self.types,self.labels=np.unique(feat,axis=0,return_inverse=True)
    def initial(self,prep):
        p=np.zeros(self.n)
        if isinstance(prep,list):p[self.lookup[np.array(prep)@POWER]]=1
        else:
            p=self.deg.astype(float)*np.all(self.obs[:,:2]==prep['counts'],axis=1);p/=p.sum()
        return p
    def rates(self,x,T):
        x=np.array([x[k] for k in KEYS]) if isinstance(x,dict) else np.asarray(x)
        f=self.types;q=np.empty(len(f));shift=(1/T-1/450)/R
        m=f[:,0]==0;q[m]=np.exp(np.log(10)*x[1]-x[0]*shift-(x[2]*f[m,1]+x[3]*f[m,2])/(R*T))
        m=f[:,0]==1;q[m]=np.exp(np.log(10)*x[5]-x[4]*shift)
        m=f[:,0]>=2;delta=-(x[2]*f[m,1]+x[3]*f[m,2]+x[8]*f[m,3]);q[m]=np.exp(np.log(10)*x[7]-x[6]*shift-delta/(2*R*T))
        oxygen=f[:,0]==3;q[oxygen]*=.2*np.exp(-(x[9]-x[6])*shift)
        return q[self.labels]
