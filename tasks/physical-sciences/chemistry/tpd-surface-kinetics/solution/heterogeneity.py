"""Unresolved surface-population mixture with shared instrument response."""
import numpy as np
import experiment as core
from experiment import thermal_knots,temperature,whiten
KEYS=core.KEYS+('fraction_B','delta_E_des_B','delta_E_ox_B','delta_E_hop_B')
SCALE=np.r_[core.SCALE,.4,3.,3.,8.]
NP=len(KEYS)
R=.00831446261815324

def components(x):
    x=np.asarray(x);a=x[:15].copy();b=a.copy();mapping=np.zeros((15,NP));mapping[:,:15]=np.eye(15)
    for shift,E,log in [(16,0,1),(17,4,5),(18,6,7)]:
        b[E]+=x[shift];b[log]-=x[shift]/(np.log(10)*R*450)
        mapping[E,shift]=1.;mapping[log,shift]=-1/(np.log(10)*R*450)
    b[9]+=x[18];mapping[9,18]=1.
    return a,b,mapping

def predict(x,run,sens=False,rtol=3e-9,atol=2e-13,full=False):
    a,b,M=components(x);pa=core.predict(a,run,sens=sens,rtol=rtol,atol=atol,full=full);pb=core.predict(b,run,sens=sens,rtol=rtol,atol=atol,full=full);f=x[15];out={}
    for key in ('q','yields','coverages','pair_fractions','flux'):
        if key in pa:out[key]=(1-f)*pa[key]+f*pb[key]
    out['sample_temperature']=pa['sample_temperature']
    if full:out['mass_error']=max(pa['mass_error'],pb['mass_error'])
    for key,value in [('dq','q'),('dyields','yields')]:
        if sens:
            da=np.pad(pa[key],[(0,0),(0,0),(0,4)]);db=pb[key]@M
            out[key]=(1-f)*da+f*db;out[key][:,:,15]+=pb[value]-pa[value]
        else:out[key]=np.empty((*out[value].shape,0))
    return out

def design(x,run,pred):
    t=np.asarray(run['times']);v=t/t[-1]-.5;q=pred['q'];dq=pred['dq'];N=len(t)
    A=np.zeros((N,3,5));A[:,0,0]=1200*q[:,0]+x[14]*850*q[:,1];A[:,1,0]=850*q[:,1];A[:,2,0]=run['argon_reference_pA']
    A[:,0,1]=1.;A[:,0,2]=v;A[:,1,3]=1.;A[:,1,4]=v;dA=np.zeros((N,3,5,NP))
    if dq.shape[-1]:
        dA[:,0,0]=1200*dq[:,0]+x[14]*850*dq[:,1];dA[:,0,0,14]+=850*q[:,1];dA[:,1,0]=850*dq[:,1]
    return A,dA
