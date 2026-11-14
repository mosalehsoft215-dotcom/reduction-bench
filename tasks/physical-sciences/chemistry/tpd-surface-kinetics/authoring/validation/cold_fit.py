"""Public-input fit from the midpoint or a deterministic displaced start."""
from pathlib import Path
import os,sys,json,time
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'solution'))
from fit import fit
from heterogeneity import KEYS

if __name__=='__main__':
    case=json.loads((ROOT/'environment/data/case.json').read_text());data=json.loads((ROOT/'environment/data/measurements.json').read_text());seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
    box=np.asarray([case['parameter_bounds'][k] for k in KEYS]);start=box.mean(axis=1)
    if seed:
        rng=np.random.default_rng(seed);start=box[:,0]+(.2+.6*rng.random(len(KEYS)))*(box[:,1]-box[:,0])
    t=time.perf_counter();answer,obj,result=fit(case,data,start=start,max_nfev=100)
    out=dict(seed=seed,start=dict(zip(KEYS,start.tolist())),result=answer,success=bool(result.success),message=result.message,elapsed_s=time.perf_counter()-t)
    Path(sys.argv[1]).write_bytes((json.dumps(out,indent=2)+'\n').encode());print('FIT_COMPLETE',out['success'],answer['deviance'],out['elapsed_s'],flush=True)
