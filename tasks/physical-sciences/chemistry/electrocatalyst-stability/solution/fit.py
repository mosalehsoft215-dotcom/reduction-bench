"""Bounded joint inverse fit with latent equilibrium coverages."""
import json,time
import numpy as np
from scipy.optimize import least_squares
from model import F,FRT,rates,jacobian

SCALE=np.array([.05,1.,1.,1.,100.,1e-9,1.,100.,1e-4])
LOG=np.array([False,True,False,False,True,True,True,True,True])
def encode(p):
    v=np.asarray(p)/SCALE
    return np.where(LOG,np.log(np.maximum(v,1e-30)),v)
def decode(v):return np.where(LOG,np.exp(np.minimum(v,100)),v)*SCALE

def fit(case,records,mechanism,starts=5):
    conditions={s['id']:s for s in case['conditions']}
    cond=[conditions[r['condition_id']] for r in records]
    bulk=np.array([s['bulk_M'] for s in cond]);rpm=np.array([s['rpm'] for s in cond]);film=np.array([s['film_cm'] for s in cond])
    u=np.array([r['command_V'] for r in records]);curr=np.array([r['current_A_cm2'] for r in records]);csig=np.array([r['current_sigma_A_cm2'] for r in records])
    f=np.array([v['frequency_Hz'] for v in records[0]['spectrum']]);n=len(records)
    z=np.array([[v['real_ohm_cm2']+1j*v['imag_ohm_cm2'] for v in r['spectrum']] for r in records])
    sig=np.array([[v['sigma_ohm_cm2'] for v in r['spectrum']] for r in records])
    def residual(v,return_predictions=False):
        p=decode(v[:9]);E,k,a,g,kc,G,km,R,C=p;th=v[9:]
        transfer=km*np.sqrt(rpm/1000);h=1000*G/film
        c=transfer/(transfer+h*kc*th**mechanism)
        cat=kc*bulk*c*th**mechanism
        current=F*G*mechanism*cat
        phi=FRT*(u+R*current-E)
        eta=phi-g*(2*th-1);ef=np.exp(-a*eta);er=np.exp((1-a)*eta)
        ve=k*((1-th)*ef-th*er)
        vp=k*(-a*(1-th)*ef-(1-a)*th*er)
        vt=k*((-1+2*a*g*(1-th))*ef-(1-2*(1-a)*g*th)*er)
        ct=kc*bulk*c*mechanism*th**(mechanism-1);cc=kc*bulk*th**mechanism
        J=np.zeros((n,3,3));J[:,0,0]=-1/(R*C)+FRT*F*G*vp/C;J[:,0,1]=FRT*F*G*vt/C
        J[:,1,0]=vp;J[:,1,1]=vt-mechanism*ct;J[:,1,2]=-mechanism*cc
        J[:,2,1]=-h*kc*c*mechanism*th**(mechanism-1);J[:,2,2]=-transfer-h*kc*th**mechanism
        matrices=2j*np.pi*f[None,:,None,None]*np.eye(3)-J[:,None,:,:]
        b=np.zeros((n,len(f),3,1));b[:,:,0,0]=FRT/(R*C)
        res=np.linalg.solve(matrices,b)[:,:,0,0]
        zp=R/(1-res/FRT)
        # Equality residual makes DC states actual equilibria, independently of measurement noise.
        eq=(ve-mechanism*cat)/1e-7
        rr=np.r_[eq,(current-curr)/csig,((zp-z)/sig).real.ravel(),((zp-z)/sig).imag.ravel()]
        if return_predictions:return rr,th,np.column_stack([phi,th,c]),current,zp
        return rr
    bounds=np.array([case['parameter_bounds'][key] for key in case['parameter_order']])
    low=np.r_[encode(bounds[:,0]),np.full(n,.002)]
    high=np.r_[encode(bounds[:,1]),np.full(n,.998)]
    rng=np.random.default_rng(86413+mechanism);best=None
    start=time.monotonic()
    for j in range(starts):
        # Data-driven coverage initialization from the DC catalytic law.
        p0=decode((low[:9]+high[:9])/2)
        if j:
            p0=decode(np.clip(encode(p0)+rng.normal(0,.15,9),low[:9]+.01,high[:9]-.01))
        kc,G,km=p0[4:7];tr=km*np.sqrt(rpm/1000);h=1000*G/film
        tpower=curr/(F*G*mechanism*kc*bulk-curr*h*kc/tr)
        theta=np.clip(np.maximum(tpower,.001)**(1/mechanism),.005,.99)
        ans=least_squares(residual,np.r_[encode(p0),theta],bounds=(low,high),max_nfev=2500,
                          ftol=1e-10,xtol=1e-10,gtol=1e-6,x_scale='jac',diff_step=2e-6)
        score=float(np.sum(ans.fun**2))
        print('fit',mechanism,j,'chi2',score,'eval',ans.nfev,'seconds',time.monotonic()-start,flush=True)
        if best is None or score<best[0]:best=(score,ans)
        if score<2*n*(1+2*len(f)):break
    score,ans=best
    rr,th,x,ip,zp=residual(ans.x,True)
    return {'chi2':score,'parameters':decode(ans.x[:9]).tolist(),'theta':th.tolist(),
            'states':x.tolist(),'jacobian_singular_values':np.linalg.svd(ans.jac,compute_uv=False).tolist()}
