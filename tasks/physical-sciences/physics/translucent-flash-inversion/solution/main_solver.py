import json, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.special import jv, jnp_zeros, ive
from numpy.polynomial.legendre import leggauss

DATA=Path('/app/data'); OUT=Path('/app/output'); OUT.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(DATA/'flash_runs.csv'); st=pd.read_csv(DATA/'detector_step.csv')

# ---------- standalone calibrations ----------
def fit_step(y):
    t=st.time_s.values; y=np.asarray(y)
    def f(p):
        tau,t0,A,b=p; u=np.clip(t-t0,0,None); return b+A*np.where(t>=t0,1-np.exp(-u/tau),0)-y
    return least_squares(f,[.02,-.0004,1,0],bounds=([.001,-.003,.8,-.02],[.1,.003,1.2,.02]),x_scale='jac').x
tau_m=float(fit_step(st.monitor_normalized_signal.values)[0])     # monitor chain: no range switch, applies to the runs
tau_cal=float(fit_step(st.detector_normalized_signal.values)[0])  # rear chain on the 1e5 V/A range only; the runs used 1e6 V/A
ff=float(st.flash_blocked_beam_ratio.mean()); fi=float(st.ir_blocked_beam_ratio.mean())

def slab(mu,Rr,L):
    """incoherent slab: transmission, absorbed fraction"""
    e=np.exp(-mu*L); return (1-Rr)**2*e/(1-Rr**2*e*e),(1-Rr)*(1-e)/(1-Rr*e)
meta=df.groupby('run_id',sort=False).first().reset_index(); Lm=meta.thickness_mm.values; wl=meta.laser_wavelength_nm.values
Tf=np.clip((meta.flash_transmission_ratio.values-ff)/(1-ff),1e-10,1); Ti=np.clip((meta.ir_transmission_ratio.values-fi)/(1-fi),1e-12,1)
mu0,pexp,Rf=least_squares(lambda p:np.log(slab(p[0]*(940./wl)**p[1],p[2],Lm)[0])-np.log(Tf),[1.25,.6,.1],bounds=([.1,-1.5,0],[10,2.5,.5]),x_scale='jac').x
RI=.09
kap=float(least_squares(lambda q:np.log(slab(q[0],RI,Lm)[0])-np.log(Ti),[6],bounds=([.1],[50])).x[0])

# ---------- in-plane: Neumann-Bessel modes of the insulated disk ----------
R=float(df.disk_radius_mm.iloc[0]); alpha_r=1.31; M=16; P=16; NQ=500
def gl(a,b,n):
    x,w=leggauss(n); return .5*(b-a)*x+.5*(b+a),.5*(b-a)*w
ms=[0]; ks=[0.]
for m in range(M+1):
    for root in jnp_zeros(m,P):
        if m==0: ms.append(0); ks.append(root/R)
        else: ms+=[m,-m]; ks+=[root/R]*2
ms=np.array(ms); ks=np.array(ks); am=np.abs(ms)
norm=np.where(ks==0,np.pi*R*R,np.pi*R*R*(1-(am/np.where(ks==0,1,ks*R))**2)*jv(am,ks*R)**2)
mask_defs={'M1':(.75,3.20,8.,72.),'M2':(1.15,4.10,132.,66.),'M3':(1.45,4.45,244.,61.)}
def mask_avg(ri,ro,cen,wid):
    r,w=gl(ri,ro,NQ); rad=(jv(am[:,None],ks[:,None]*r)*r*w).sum(1); h=np.deg2rad(wid/2)
    ang=np.where(ms==0,2*h,2*np.sin(ms*h)/np.where(ms==0,1,ms))*np.exp(1j*ms*np.deg2rad(cen))
    return rad*ang   # each detector integrates over its aperture
ALLMASKS=('M1','M2','M3'); DAQS=('D1','D2','D3')
DP=np.stack([mask_avg(*mask_defs[k]) for k in ALLMASKS],axis=1)
prog={}
for p,g0 in meta.groupby('pulse_program'):
    row=g0.iloc[0]; d=float(row.source_offset_mm); wb=float(row.beam_1e2_radius_mm)
    r,w=gl(0,R,2*NQ); g=np.exp(-2*(r-d)**2/wb**2); x=4*r*d/wb**2; Z=2*np.pi*np.sum(g*ive(0,x)*r*w)
    rad=np.array([np.sum(g*ive(m,x)*jv(m,k*r)*r*w) for m,k in zip(am,ks)])
    prog[p]=dict(src=2*np.pi*rad/(Z*norm),angle=float(row.source_angle_camera_deg),mu=mu0*(940./float(row.laser_wavelength_nm))**pexp)

# ---------- through-thickness: closed form with coated, exchanging faces ----------
def axial(L,mu,Rf,kap,Ri,alpha,sv,b):
    """emission-weighted temperature per unit incident flash power for one in-plane mode:
    alpha T''-sv T=-S, alpha T'(0)=b T(0), -alpha T'(L)=b T(L), b=h+i omega c"""
    g=np.sqrt(sv/alpha); eL=np.exp(-mu*L); emL=np.exp(mu*L)
    a1=mu*(1-Rf)/(1-Rf**2*eL**2); a2=a1*Rf*eL**2; den=alpha*(g*g-mu*mu); A1=a1/den; A2=a2/den
    th0=A1+A2; dth0=mu*(A2-A1); thL=A1*eL+A2*emL; dthL=mu*(A2*emL-A1*eL)
    E=np.exp(-g*L); Pp=-(alpha*g+b); Qq=(alpha*g-b)*E; r1=b*th0-alpha*dth0; r2=b*thL+alpha*dthL; det=Pp*Pp-Qq*Qq
    C1=(Pp*r1-Qq*r2)/det; C2=(Pp*r2-Qq*r1)/det
    ek=np.exp(-kap*L); K=kap*(1-Ri)/(1-Ri**2*ek**2)
    T1=(eL-ek)/(kap-mu)+Ri*ek*(1-np.exp(-(kap+mu)*L))/(kap+mu)
    T2=(emL-ek)/(kap+mu)+Ri*ek*(1-np.exp((mu-kap)*L))/(kap-mu)
    T3=(E-ek)/(kap-g)+Ri*ek*(1-np.exp(-(kap+g)*L))/(kap+g)
    T4=(1-np.exp(-(kap+g)*L))/(kap+g)+Ri*ek*(ek-E)/(g-kap)
    return K*(A1*T1+A2*T2+C1*T3+C2*T4)

# ---------- data ----------
freq=df.modulation_hz.values; om=2*np.pi*freq
mon=df.incident_monitor_re.values+1j*df.incident_monitor_im.values; rear=df.rear_detector_re.values+1j*df.rear_detector_im.values; sig=df.measurement_sigma.values
drive=mon*(1+1j*om*tau_m)
sess=df.session_id.values.astype(str); off=df.measured_clock_offset_s.values; leakfac=df.leak_channel_factor.values
Tfrow=np.clip((df.flash_transmission_ratio.values-ff)/(1-ff),0,1)
Ls,iL=np.unique(df.thickness_mm.values,return_inverse=True); DLs,iD=np.unique(df.drive_level.values,return_inverse=True)
Fs,iF=np.unique(freq,return_inverse=True); PS,iP=np.unique(df.pulse_program.values,return_inverse=True)
SS,iS=np.unique(sess,return_inverse=True); iQ=np.array([DAQS.index(v) for v in df.daq_channel.values.astype(str)])
RUNS,iR=np.unique(df.run_id.values.astype(str),return_inverse=True)
og=2*np.pi*Fs[None,None,:,None]; Lg=Ls[:,None,None,None]; kk=alpha_r*ks**2

def thermal_all(alpha0,loss0,h,c,aslope,lheat,phis):
    out=np.empty((len(SS),len(PS),len(Ls),len(DLs),len(Fs),3),complex)
    sv0=1j*og+kk; bnd=h+1j*og*c
    for ip,p in enumerate(PS):
        mu=prog[p]['mu']; Tp,Ap=slab(mu,Rf,Lg)
        q=DLs[None,:,None,None]*Ap
        al=alpha0*(1+aslope*q); lv=loss0+lheat*q*q
        sv=lv+sv0
        TH=axial(Lg,mu,Rf,kap,RI,al,sv,bnd)
        for js in range(len(SS)):
            rot=prog[p]['src']*np.exp(-1j*ms*np.deg2rad(prog[p]['angle']+phis[js]))
            out[js,ip]=TH@(rot[:,None]*DP)
    return out

def model(x,maps):
    alpha0,loss0,h,c,aslope,lheat,phA,phB,phC,bA,bB,bC,leak,tau_run=x
    th=thermal_all(alpha0,loss0,h,c,aslope,lheat,(phA,phB,phC))
    mk=np.array([[ALLMASKS.index(maps[s][d]) for d in DAQS] for s in SS])[iS,iQ]
    therm=th[iS,iP,iL,iD,iF,mk]; bias=np.array([bA,bB,bC])[iS]
    Hr=1/(1+1j*om*tau_run)
    return (drive*therm+leak*leakfac*Tfrow*drive)*Hr*np.exp(-1j*om*(off+bias))

def resid(x,maps,sub=None):
    b=model(x,maps); w=1/sig**2 if sub is None else sub/sig**2
    g=np.maximum(np.bincount(iR,w*np.real(np.conj(b)*rear),len(RUNS))/np.maximum(np.bincount(iR,w*np.abs(b)**2,len(RUNS)),1e-300),1e-8)
    z=(g[iR]*b-rear)/sig
    if sub is not None: z=z[sub]
    return np.r_[z.real,z.imag]

# ---------- B/C cable permutations: coarse search on mid thicknesses at drive 1.0 ----------
mapA={'D1':'M1','D2':'M2','D3':'M3'}
perms=[dict(zip(DAQS,p)) for p in itertools.permutations(ALLMASKS)]
nom=np.array([.75,.11,.015,.05,.15,.02,8.,0.,0.,0.,0.,0.,.004,2*tau_cal])
chosen={'A':mapA}
for S in ('B','C'):
    sub=(sess==S)&(np.abs(df.drive_level.values-1.0)<1e-9)&(df.thickness_mm.values>.65)&(df.thickness_mm.values<1.2)&(freq>.12)
    best=(1e300,None,None)
    for mp in perms:
        for phi in np.arange(-28,29,7):
            maps={'A':mapA,'B':mapA,'C':mapA}; maps[S]=mp; x=nom.copy(); x[7 if S=='B' else 8]=phi
            c=np.mean(resid(x,maps,sub)**2)
            if c<best[0]: best=(c,mp,phi)
    chosen[S]=best[1]; nom[7 if S=='B' else 8]=best[2]

# ---------- global fit ----------
maps={'A':mapA,'B':chosen['B'],'C':chosen['C']}
lo=[.3,0,0,0,-.2,0,-40,-40,-40,-.003,-.003,-.003,0,.005]
hi=[1.3,.4,.08,.2,.8,.15,40,40,40,.003,.003,.003,.02,.1]
xs=[.1,.03,.01,.02,.1,.02,10,10,10,.0005,.0005,.0005,.002,.005]
fit=least_squares(lambda x:resid(x,maps),nom,bounds=(lo,hi),x_scale=xs,max_nfev=300)
for S in ('B','C'):   # mapping refinement with fitted continuous parameters
    sub=sess==S; best=(1e300,None)
    for mp in perms:
        tm=dict(maps); tm[S]=mp; c=np.mean(resid(fit.x,tm,sub)**2)
        if c<best[0]: best=(c,mp)
    maps[S]=best[1]
fit=least_squares(lambda x:resid(x,maps),fit.x,bounds=(lo,hi),x_scale=xs,max_nfev=300)
rms=float(np.sqrt(np.mean(fit.fun**2)))
cov=np.linalg.pinv(fit.jac.T@fit.jac)*np.sum(fit.fun**2)/(len(fit.fun)-len(fit.x)-len(RUNS))
se=np.sqrt(np.diag(cov))
print('maps',maps); print('x',[round(float(v),6) for v in fit.x]); print('se alpha %.5f loss %.6f tau %.6f'%(se[0],se[1],se[13])); print('rms',rms)
out={'alpha_mm2_s':float(fit.x[0]),'optical_penetration_mm':float(1/mu0),'loss_s_inv':float(fit.x[1]),'detector_tau_s':float(fit.x[13])}
with open(OUT/'result.json','w') as f: json.dump(out,f,indent=2)
