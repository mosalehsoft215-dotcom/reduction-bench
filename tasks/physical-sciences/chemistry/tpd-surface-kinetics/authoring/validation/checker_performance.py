"""Measure the compiled checker against the captured pre-optimization engine.

Both use the same independent event graph and RK45 tolerances. Rates are
evaluated once per distinct event feature, and fluxes are accumulated in a
compiled loop. No model or acceptance allowance is changed.
"""
from pathlib import Path
import os,sys,json,time,hashlib,importlib.util
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tests'))
from experiment_check import experiment
saved=ROOT/'authoring/validation/independent_inputs/experiment_check.py'
spec=importlib.util.spec_from_file_location('captured_checker',saved);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def main():
    a=json.loads(Path(sys.argv[1]).read_text());point=a.get('result',a);case=json.loads((ROOT/'environment/data/case.json').read_text());data=json.loads((ROOT/'environment/data/measurements.json').read_text());records=[]
    for run in (data['runs'][0],case['forecasts'][0]):
        start=time.perf_counter();baseline=old.experiment(point['parameters'],run,full=True);middle=time.perf_counter();fast=experiment(point['parameters'],run,full=True);end=time.perf_counter()
        errors={k:float(np.max(abs(np.asarray(fast[k])-np.asarray(baseline[k])))) for k in baseline}
        assert max(errors.values())<2e-8,errors
        records.append(dict(run=run['id'],baseline_s=middle-start,compiled_s=end-middle,absolute_differences=errors))
    files=[ROOT/'tests/checker_fast.py',ROOT/'tests/experiment_check.py',ROOT/'tests/lattice_check.py',saved]
    out=dict(provenance=__doc__,measurements=records,source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    Path(sys.argv[2]).write_bytes((json.dumps(out,indent=2)+'\n').encode());print(json.dumps(out),flush=True)
if __name__=='__main__':main()
