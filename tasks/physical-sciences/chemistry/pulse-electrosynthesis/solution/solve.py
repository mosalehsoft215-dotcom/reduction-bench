"""Compute a paired robust design from the public case only."""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from fast_model import trajectory
from design import responses,information,evaluate,summary

DATA=Path('/app/data'); OUT=Path('/app/output/protocol.json'); START=time.monotonic()
def log(s): print(f'[{time.monotonic()-START:.1f}s] {s}',flush=True)

def main():
    c=json.loads((DATA/'cell.json').read_text()); scenarios=json.loads((DATA/'scenarios.json').read_text())
    params=np.array([[p[k] for k in c['scenario_parameter_order']] for p in scenarios.values()])
    # Each of four epochs contains six generation/recovery cycles. Vary the
    # epoch length, duty fraction and both potentials, independently in both runs.
    def decode(x):
        out=[]
        for run in range(2):
            z=x[run*16:(run+1)*16].reshape(4,4)
            lengths=6*z[:,0]/z[:,0].sum()
            schedule=[]
            for length,(_,duty,on,off) in zip(lengths,z):
                schedule.extend([[length*duty/6,on],[length*(1-duty)/6,off]]*6)
            out.append(np.array(schedule))
        return out
    x=np.array([[1,.80,-1.25,-1.055]]*8).ravel()
    # Break the exchange symmetry of the two fresh-cell experiments.
    x[2]=-1.22;x[18]=-1.27;x[1]=.75;x[17]=.85
    bounds=[b for _ in range(8) for b in [(.45,2.0),(.35,.91),(-1.30,-1.10),(-1.12,-.80)]]
    active={1}; cache={}; n=24;dt=.006
    cap=c['maximum_charge_equivalents']-.004
    floor=c['minimum_yield']+.003
    def calc(v):
        key=v.tobytes()
        if key in cache:return cache[key]
        programs=decode(v); metrics=[];scores=[]
        for j,p in enumerate(params):
            if j in active:
                met,bins,jac=responses(programs,p,n,dt)
                scores.append(np.linalg.eigvalsh(information(bins,jac,c))[0])
            else: met=np.array([trajectory(a,p,n,dt)[0] for a in programs])
            metrics.append(met)
        metrics=np.array(metrics)
        constraints=np.r_[metrics[:,:,0].ravel()-floor,cap-metrics[:,:,2].ravel(),
                          np.concatenate([a[:,0] for a in programs])-.0201]
        cache[key]=(min(scores),constraints)
        return cache[key]
    best=[None,-np.inf];count=[0]
    def remember(v):
        score,cons=calc(v)
        if cons.min()>-2e-6 and score>best[1]: best[:]=[v.copy(),score]
        count[0]+=1
        if count[0]%10==0:log(f'iteration {count[0]} information={score:.6f} constraint={cons.min():.4g}')
    for stage in range(3):
        # Objective and constraint finite differences visit the same controls.
        # Reuse these expensive sensitivities within a grid/scenario stage only.
        cache.clear()
        result=minimize(lambda v:-np.log(max(calc(v)[0],1e-10)),x,method='SLSQP',bounds=bounds,
              constraints=[{'type':'ineq','fun':lambda v:calc(v)[1]}],callback=remember,
              options={'maxiter':100,'ftol':2e-7,'eps':2e-5})
        remember(result.x)
        x=best[0] if best[0] is not None else result.x
        all_records=evaluate(decode(x),params,c,n,dt); scores,metrics=summary(all_records)
        worst=int(np.argmin(scores));log(f'stage {stage}: {result.message}; all information {scores.tolist()}')
        if worst not in active:
            active.add(worst);best[:]=[None,-np.inf];continue
        if stage==0:
            n=40;dt=.002;best[:]=[None,-np.inf];continue
        break
    if best[0] is None: raise RuntimeError('No feasible design found')
    programs=decode(best[0])
    # Final predictions on a fine grid: at n=120, dt=2.5e-4 the charge bins sit up to
    # 9.5e-5 from the converged physics; at n=360, dt=2.5e-5 within 1.0e-5
    # (authoring/evidence/reference_resolution.txt).
    records=evaluate(programs,params,c,n=360,dt=.000025,h=.0015)
    scores,metrics=summary(records)
    if metrics[:,:,0].min()<c['minimum_yield'] or metrics[:,:,2].max()>c['maximum_charge_equivalents']:
        raise RuntimeError('Refined solution violates operating limits')
    report={'programs':{label:[{'duration_s':float(d),'potential_V':float(e)} for d,e in a]
                        for label,a in zip(c['experiments'],programs)},'predictions':{}}
    for name,(outcomes,bins,info) in zip(scenarios,records):
        report['predictions'][name]={'experiments':{label:{'outcomes':dict(zip(c['metrics_order'],map(float,m))),
                'charge_bins':b.tolist()} for label,m,b in zip(c['experiments'],outcomes,bins)},
                'profile_information':info.tolist()}
    OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    log(f'final worst information={scores.min():.7f}; min yield={metrics[:,:,0].min():.7f}; max charge={metrics[:,:,2].max():.7f}')

if __name__=='__main__':main()
