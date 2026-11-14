"""Reproduce a synthetic shortlist of experimentally admissible heating programs.

The shortlist spans a preliminary near-equal-yield screen of four different
ramp/hold histories. Selection here uses generating parameters solely to author
the synthetic case. Inference and the final minimax-regret choice use QMS data.
"""
from pathlib import Path
import os,sys,json,time,copy
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'solution'))
from heterogeneity import KEYS,predict

def main():
    started=time.perf_counter();case=json.loads((ROOT/'environment/data/case.json').read_text());truth=json.loads((Path(__file__).parent/'generating_parameters.json').read_text());x=np.array([truth['parameters'][k] for k in KEYS])
    base=copy.deepcopy(case['forecasts'][0]);base['id']='early_hot';base['times']=[0,10,15,30,55,75,100]
    def value(run):return float(predict(x,run,rtol=1e-9)['yields'][-1,1])
    Y=value(base);programs=[base];records=[dict(id=base['id'],yield_at_generating_parameters=Y)]
    templates=[('prehold',lambda T:[[0,370],[20,T],[65,T],[80,475],[100,475]],2e-4),('double_hold',lambda T:[[0,370],[10,465],[25,465],[40,T],[85,T],[100,465]],-1e-4),('late_ramp',lambda T:[[0,370],[55,T],[80,T],[100,495]],1e-4)]
    for name,make,offset in templates:
        def run(T):
            p=make(float(T));times=sorted(set([0.,10.,30.,60.,90.,100.]+[float(k[0]) for k in p]));return dict(id=name,counts=[4,4],initial=dict(counts=[4,4],preparation_temperature_K=360.),program=p,times=times)
        target=Y+offset;grid=np.linspace(380,485,12);ys=[value(run(T)) for T in grid];bracket=None
        for l,h,yl,yh in zip(grid,grid[1:],ys,ys[1:]):
            if (yl-target)*(yh-target)<=0:bracket=(l,h);break
        if bracket is None:raise RuntimeError(str((name,target,list(zip(grid,ys)))))
        T=round(float(brentq(lambda T:value(run(T))-target,*bracket,xtol=2e-6)),3);selected=run(T);programs.append(selected);records.append(dict(id=name,screened_temperature_K=T,target_yield=target,grid=grid.tolist(),grid_yields=ys,yield_at_generating_parameters=value(selected)))
        print(name,T,records[-1]['yield_at_generating_parameters'],flush=True)
    case['forecasts']=programs
    Path(sys.argv[1]).write_bytes((json.dumps(case,indent=2)+'\n').encode())
    Path(sys.argv[2]).write_bytes((json.dumps(dict(provenance=__doc__,screen=records,elapsed_s=time.perf_counter()-started),indent=2)+'\n').encode())

if __name__=='__main__':main()
