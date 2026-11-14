"""Time independent likelihood and four-program verification in isolation.

The report includes the actual point and numerical comparison with its source
fit. Runtime is a forward-cost measurement, not a frontier-agent trial.
"""
from pathlib import Path
import os,sys,json,time,resource
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tests'))
from experiment_check import deviance,forecast

def main():
    source=Path(sys.argv[1]);target=Path(sys.argv[2]);a=json.loads(source.read_text());p=a.get('result',a);case=json.loads((ROOT/'environment/data/case.json').read_text());data=json.loads((ROOT/'environment/data/measurements.json').read_text());start=time.perf_counter()
    D,parts=deviance(p['parameters'],p['run_nuisance'],data);middle=time.perf_counter();rows,_=forecast(p['parameters'],case);end=time.perf_counter()
    assert abs(D-p['deviance'])<.03
    report=dict(provenance=__doc__,source_deviance=p['deviance'],checker_deviance=D,absolute_deviance_difference=abs(D-p['deviance']),likelihood_s=middle-start,four_forecasts_s=end-middle,total_s=end-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,likelihood_parts=parts,final_CO2_yields=[v['snapshots'][-1]['gas_yields'][1] for v in rows])
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(report),flush=True)
if __name__=='__main__':main()
