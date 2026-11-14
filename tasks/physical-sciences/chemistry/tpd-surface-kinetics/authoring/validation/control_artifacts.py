"""Build schema-valid scientific controls around a validated full artifact.

The unmodified comparisons remain as measured witnesses. Changed point or
pair predictions are recomputed, and regret is recomputed from all witnesses.
"""
from pathlib import Path
import os,sys,json,copy,argparse
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tests'))
from experiment_check import forecast

def main():
    parser=argparse.ArgumentParser();parser.add_argument('baseline',type=Path);parser.add_argument('controls',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--point',action='store_true');args=parser.parse_args()
    base=json.loads(args.baseline.read_text());report=json.loads(args.controls.read_text());case=json.loads((ROOT/'environment/data/case.json').read_text());args.output.mkdir(parents=True,exist_ok=True)
    variants={report['mode']:report} if args.point else report['alternatives']
    for name,v in variants.items():
        a=copy.deepcopy(base)
        if args.point:
            for k in ('parameters','run_nuisance','deviance'):a[k]=v[k]
            a['forecasts']=forecast(a['parameters'],case)[0]
        else:
            item=next(w for w in a['comparisons'] if w['pair']==report['pair'])
            for endpoint in ('lower','upper'):item[endpoint]={k:v[endpoint][k] for k in ('parameters','run_nuisance','deviance','yield_difference')}
        regret=np.zeros(4)
        for item in a['comparisons']:
            i,j=item['pair'];regret[i]=max(regret[i],-item['lower']['yield_difference']);regret[j]=max(regret[j],item['upper']['yield_difference'])
        a['decision']=dict(worst_case_regret=regret.tolist(),chosen_id=case['forecasts'][int(regret.argmin())]['id'])
        (args.output/(name+'.json')).write_bytes((json.dumps(a,indent=2,allow_nan=False)+'\n').encode())
if __name__=='__main__':main()
