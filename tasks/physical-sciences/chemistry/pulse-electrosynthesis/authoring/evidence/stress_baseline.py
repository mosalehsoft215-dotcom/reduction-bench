"""Multistart nonlinear drifting pulses: a stronger low-dimensional competitor."""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'solution'))
from design import responses,information,evaluate,summary
from fast_model import trajectory
c=json.load((ROOT/'environment/data/cell.json').open());s=json.load((ROOT/'environment/data/scenarios.json').open())
p=np.array([[v[k] for k in c['scenario_parameter_order']] for v in s.values()])
baseline=json.load((Path(__file__).parent/'calibration_results.json').open())[-1]
seed=[]
for program in baseline['programs']:
    a=np.array(program);seed.extend([a[0,0]/.25,a[-2,0]/.25,a[0,1],a[-2,1],a[1,1],a[-1,1],1.])
results=[]
for exponent in [.3,1.,3.]:
    start=time.monotonic();x=np.array(seed);x[[6,13]]=exponent
    bounds=[b for _ in range(2) for b in [(.4,.919),(.4,.919),(-1.3,-1.1),(-1.3,-1.1),(-1.12,-.8),(-1.12,-.8),(.15,5.)]]
    def decode(v):
        programs=[]
        for z in v.reshape(2,7):
            a=[]
            for t in np.linspace(0,1,24)**z[6]:
                duty=(1-t)*z[0]+t*z[1];on=(1-t)*z[2]+t*z[3];off=(1-t)*z[4]+t*z[5]
                a.extend([[.25*duty,on],[.25*(1-duty),off]])
            programs.append(np.array(a))
        return programs
    cache={};best=[None,-np.inf]
    def calc(v):
        key=v.tobytes()
        if key in cache:return cache[key]
        programs=decode(v);m,b,j=responses(programs,p[1],24,.006)
        score=np.linalg.eigvalsh(information(b,j,c))[0]
        metrics=np.array([m if i==1 else [trajectory(a,pars,24,.006)[0] for a in programs] for i,pars in enumerate(p)])
        cons=np.r_[metrics[:,:,0].ravel()-.243,.696-metrics[:,:,2].ravel()]
        cache[key]=(score,cons)
        if cons.min()>-1e-6 and score>best[1]:best[:]=[v.copy(),score]
        return cache[key]
    fit=minimize(lambda v:-np.log(max(calc(v)[0],1e-8)),x,method='SLSQP',bounds=bounds,
        constraints=[{'type':'ineq','fun':lambda v:calc(v)[1]}],options={'maxiter':80,'ftol':1e-6,'eps':2e-5})
    v=best[0] if best[0] is not None else fit.x
    records=evaluate(decode(v),p,c,n=64,dt=.001);scores,metrics=summary(records)
    result=dict(start_exponent=exponent,score=float(scores.min()),min_yield=float(metrics[:,:,0].min()),max_charge=float(metrics[:,:,2].max()),programs=[a.tolist() for a in decode(v)],seconds=time.monotonic()-start,termination=fit.message)
    results.append(result);(Path(__file__).parent/'stress_baseline_results.json').write_text(json.dumps(results,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='programs'}),flush=True)
