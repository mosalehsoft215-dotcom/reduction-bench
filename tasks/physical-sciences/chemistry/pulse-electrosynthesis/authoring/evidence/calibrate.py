import sys,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'solution'))
from design import responses,information,evaluate,summary
from fast_model import trajectory
c=json.load((ROOT/'environment/data/cell.json').open());scenarios=json.load((ROOT/'environment/data/scenarios.json').open())
p=np.array([[v[k] for k in c['scenario_parameter_order']] for v in scenarios.values()])
results=[]
for family in ['periodic12','periodic20','periodic24','drift24']:
    start=time.monotonic(); ncycles=int(family[-2:]);drift=family.startswith('drift')
    bounds=[b for _ in range(4 if drift else 2) for b in [(.45,.91),(-1.30,-1.1),(-1.12,-.8)]]
    x=np.array([[.8,-1.25,-1.055]]*(4 if drift else 2)).ravel();x[1]= -1.24;x[-2]=-1.26
    def decode(v):
        z=v.reshape(2,2 if drift else 1,3);out=[]
        for r in range(2):
            rows=[]
            for t in np.linspace(0,1,ncycles):
                duty,on,off=(1-t)*z[r,0]+t*z[r,-1]
                rows.extend([[6/ncycles*duty,on],[6/ncycles*(1-duty),off]])
            out.append(np.array(rows))
        return out
    cache=[None,None];best=[None,-np.inf]
    def calc(v):
        if cache[0] is not None and np.array_equal(v,cache[0]):return cache[1]
        prog=decode(v);m,b,j=responses(prog,p[1],24,.006)
        score=np.linalg.eigvalsh(information(b,j,c))[0]
        metrics=np.array([m if i==1 else [trajectory(a,pars,24,.006)[0] for a in prog] for i,pars in enumerate(p)])
        cons=np.r_[metrics[:,:,0].ravel()-.243,.696-metrics[:,:,2].ravel()]
        cache[:]=[v.copy(),(score,cons)]
        if cons.min()>-1e-6 and score>best[1]:best[:]=[v.copy(),score]
        return cache[1]
    fit=minimize(lambda v:-np.log(max(calc(v)[0],1e-9)),x,method='SLSQP',bounds=bounds,
        constraints=[{'type':'ineq','fun':lambda v:calc(v)[1]}],options={'maxiter':75,'ftol':1e-6,'eps':2e-5})
    v=best[0] if best[0] is not None else fit.x
    records=evaluate(decode(v),p,c,n=64,dt=.001)
    scores,metrics=summary(records)
    item=dict(family=family,score=float(scores.min()),min_yield=float(metrics[:,:,0].min()),max_charge=float(metrics[:,:,2].max()),programs=[a.tolist() for a in decode(v)],seconds=time.monotonic()-start,termination=fit.message)
    results.append(item);print(json.dumps({k:v for k,v in item.items() if k!='programs'}),flush=True)
    (Path(__file__).parent/'calibration_results.json').write_text(json.dumps(results,indent=2))
