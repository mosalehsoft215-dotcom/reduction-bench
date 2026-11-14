import csv
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
HERE=Path(__file__).resolve().parent

def _read_csv(name):
    with open(HERE/name,newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

def _load_basis():
    rows=_read_csv('reference_basis_public.csv')
    xs=np.array(sorted({float(r['x_norm']) for r in rows}));ys=np.array(sorted({float(r['y_norm']) for r in rows}))
    ix={x:i for i,x in enumerate(xs)};iy={y:i for i,y in enumerate(ys)}; interps=[]
    for k in range(1,9):
        z=np.empty((len(ys),len(xs)),float)
        for r in rows:z[iy[float(r['y_norm'])],ix[float(r['x_norm'])]]=float(r[f'basis_{k}'])
        interps.append(RegularGridInterpolator((ys,xs),z,bounds_error=True))
    return interps

def _basis(interps,angle,dx,dy,cx,cy,X,Y,Q):
    th=np.deg2rad(angle);c,s=np.cos(th),np.sin(th)
    xr=cx+c*(X-cx)+s*(Y-cy)+dx;yr=cy-s*(X-cx)+c*(Y-cy)+dy
    M=np.column_stack([f(np.c_[yr,xr]) for f in interps]);return M-Q@(Q.T@M)

def _precompute():
    cand=[(int(r['candidate_id']),float(r['angle_deg']),float(r['dither_x_norm']),float(r['dither_y_norm'])) for r in _read_csv('design_candidates_public.csv')]
    scen=[(float(r['zero_clocking_offset_deg']),float(r['center_x_norm']),float(r['center_y_norm']),float(r['settling_alpha']),float(r['drift_prior'])) for r in _read_csv('design_scenarios_public.csv')]
    u=np.linspace(-1,1,56);XX,YY=np.meshgrid(u,u);X=XX.ravel();Y=YY.ravel();N=X.size
    A=np.c_[np.ones(N),X,Y];Q=np.linalg.qr(A,mode='reduced')[0]
    q=X*X+0.6*X*Y+0.35*Y*Y;q=q-Q@(Q.T@q);q/=np.sqrt(np.mean(q*q))
    qd=X**3-0.45*X*Y*Y+0.25*Y**3;qd=qd-Q@(Q.T@qd);qd/=np.sqrt(np.mean(qd*qd))
    probe=np.array([.55,-.35,.25,.4,-.2,.3,-.15,.45],float);probe/=np.linalg.norm(probe)
    interps=_load_basis();M=np.empty((len(scen),len(cand),8,8));G=np.empty((len(scen),len(cand),8));D=np.empty((len(scen),len(cand)))
    for si,(z,cx,cy,alpha,prior) in enumerate(scen):
        base=_basis(interps,z,0.,0.,cx,cy,X,Y,Q)
        for j,(_,a,dx,dy) in enumerate(cand):
            d=_basis(interps,a+z,dx,dy,cx,cy,X,Y,Q)-base
            M[si,j]=d.T@d/N;G[si,j]=d.T@q/N;D[si,j]=float(np.mean(qd*(d@probe)))
    return cand,scen,M,G,D

def policy_metrics(policy):
    cand,scen,M,G,D=_precompute();id2j={c[0]:j for j,c in enumerate(cand)};I=np.eye(8)*1e-6
    pilot=policy['pilot_plan'];pids=[id2j[int(x['candidate_id'])] for x in pilot];pus=[int(x['exposure_units']) for x in pilot]
    sig=np.array([sum(u*D[s,j] for j,u in zip(pids,pus))/sum(pus) for s in range(len(scen))])
    t1,t2=map(float,policy['branch_thresholds']);assign=np.where(sig<t1,0,np.where(sig<t2,1,2));names=['low','mid','high']
    angles=np.array([c[1] for c in cand]);dxs=np.array([c[2] for c in cand]);dys=np.array([c[3] for c in cand])
    def adiff(a,b):return abs(((a-b+180)%360)-180)
    def seg(si,items):
        alpha,prior=scen[si][3],scen[si][4];F=np.zeros((8,8));g=np.zeros(8);h=prior;pa=pdx=pdy=0.;tt=0.
        for x in items:
            j=id2j[int(x['candidate_id'])];uu=int(x['exposure_units'])
            move=adiff(angles[j],pa)/90.0+0.75*np.hypot(dxs[j]-pdx,dys[j]-pdy)/0.070710678
            eta=0.38+0.62*np.exp(-alpha*move);eu=uu*eta;tm=tt+0.5*uu
            F+=eu*M[si,j];g+=eu*tm*G[si,j];h+=eu*tm*tm;tt+=uu;pa=angles[j];pdx=dxs[j];pdy=dys[j]
        return F-np.outer(g,g)/h
    vals=[]
    for s in range(len(scen)):
        F=I+seg(s,pilot)+seg(s,policy['followup_plans'][names[int(assign[s])]])
        vals.append(float(np.linalg.eigvalsh(F)[0]))
    margin=float(min(np.min(np.abs(sig-t1)),np.min(np.abs(sig-t2))))
    return {'score':float(min(vals)),'scenario_scores':vals,'diagnostics':sig.tolist(),'assignments':assign.tolist(),'threshold_margin':margin}

def plan_score(policy):return policy_metrics(policy)['score']
