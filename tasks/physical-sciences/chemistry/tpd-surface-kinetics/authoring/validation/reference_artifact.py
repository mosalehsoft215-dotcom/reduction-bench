"""Serialize an executed oracle calibration as a calibration artifact.

This does not claim to be a fresh cold solution execution. Its provenance
names the report from which the point and witnesses were copied.
"""
from pathlib import Path
import os,sys,json,time,hashlib
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'solution'))
from heterogeneity import KEYS
from forecast import forecast_all
from solve import robust_decision

def main():
    started=time.perf_counter();source=Path(sys.argv[1]);output=Path(sys.argv[2]);a=json.loads(source.read_text());case=json.loads((ROOT/'environment/data/case.json').read_text());p=a['point']
    rows,_,_=forecast_all(np.array([p['parameters'][k] for k in KEYS]),case)
    comparisons=[{k:w[k] for k in ('pair','lower','upper')} for w in a['executions']]
    answer={k:p[k] for k in ('parameters','run_nuisance','deviance')}
    answer.update(forecasts=rows,comparisons=comparisons,decision=robust_decision(rows,comparisons))
    output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes((json.dumps(answer,indent=2,allow_nan=False)+'\n').encode())
    report=dict(provenance=__doc__,input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),artifact_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),elapsed_s=time.perf_counter()-started)
    output.with_suffix('.provenance.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(report),flush=True)
if __name__=='__main__':main()
