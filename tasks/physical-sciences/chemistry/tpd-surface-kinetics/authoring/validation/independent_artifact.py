"""Serialize a completed independent cold fit and its six prediction profiles."""
from pathlib import Path
import os,sys,json,hashlib,time
from itertools import combinations
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tests'))
from experiment_check import forecast

def main():
    started=time.perf_counter();source=Path(sys.argv[1]);output=Path(sys.argv[2]);a=json.loads(source.read_text());case=json.loads((ROOT/'environment/data/case.json').read_text());pairs=list(combinations(range(4),2))
    assert a['fit']['success'] and set(a['endpoints'])=={str(i)+'_'+str(j) for i,j in pairs}
    p=a['point'];rows,_=forecast(p['parameters'],case);comparisons=[];regret=np.zeros(4);flags={}
    for i,j in pairs:
        record=a['endpoints'][str(i)+'_'+str(j)];assert record['pair']==[i,j]
        item=dict(pair=[i,j]);flags[str(i)+'_'+str(j)]={}
        for name in ('lower','upper'):
            w=record[name];item[name]={k:w[k] for k in ('parameters','run_nuisance','deviance','yield_difference')};flags[str(i)+'_'+str(j)][name]=dict(success=w['success'],message=w['message'],iterations=w['iterations'])
        assert item['lower']['yield_difference']<=item['upper']['yield_difference']
        comparisons.append(item);regret[i]=max(regret[i],-item['lower']['yield_difference']);regret[j]=max(regret[j],item['upper']['yield_difference'])
    out=dict(parameters=p['parameters'],run_nuisance=p['run_nuisance'],deviance=p['deviance'],forecasts=rows,comparisons=comparisons,decision=dict(worst_case_regret=regret.tolist(),chosen_id=rows[int(np.argmin(regret))]['id']))
    output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    report=dict(provenance=__doc__,snapshot_method='Separate checker forward code at independently fitted shared parameters.',input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),artifact_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),endpoint_optimizer_status=flags,elapsed_s=time.perf_counter()-started)
    output.with_suffix('.provenance.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(report),flush=True)
if __name__=='__main__':main()
