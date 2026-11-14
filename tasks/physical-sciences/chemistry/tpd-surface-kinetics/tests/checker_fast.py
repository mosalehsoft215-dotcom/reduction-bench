"""Compiled event-flux accumulation for the independent checker graph."""
import numpy as np
from numba import njit
R=.00831446261815324

def coefficients(x,types,population):
    """Unique event rates are exp(intercept + inverse_temperature/T)."""
    z=types[:,0];energy=np.where(z==0,x[0],np.where(z==1,x[4],np.where(z==2,x[6],x[9])))
    logk=np.where(z==0,x[1],np.where(z==1,x[5],x[7]))
    intercept=np.log(10.)*logk+energy/(R*450.)+np.where(z==3,np.log(.2),0.)
    inverse=-energy/R
    inverse-=np.where(z==0,(x[2]*types[:,1]+x[3]*types[:,2])/R,0.)
    inverse+=np.where(z>=2,(x[2]*types[:,1]+x[3]*types[:,2]+x[8]*types[:,3])/(2*R),0.)
    if population:inverse-=np.where(z==0,x[16],np.where(z==1,x[17],x[18]))/R
    return intercept,inverse

@njit
def event_action(p,q,labels,src,dst,tags):
    dp=np.zeros(len(p));CO=0.;CO2=0.
    for e in range(len(src)):
        typ=labels[e];flow=q[typ]*p[src[e]]
        dp[src[e]]-=flow;dp[dst[e]]+=flow
        if tags[typ]==0:CO+=flow
        elif tags[typ]==1:CO2+=flow
    return dp,CO/12.,CO2/12.
