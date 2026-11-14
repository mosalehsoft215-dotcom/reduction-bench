"""Fit explicit scientific approximations using the oracle forward engine.

These are method-discrimination controls, not independent calibration. Both
start at the full-model fit. Uniform preparation omits Gibbs weighting; the
weak-heterogeneity control fixes four coefficients to the declared values.
Each fitted point is then evaluated with the original full model and data.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import argparse,copy,json,sys,time,hashlib
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'solution'))
from heterogeneity import KEYS,SCALE
from fit import Objective

def main():
    parser=argparse.ArgumentParser();parser.add_argument('point',type=Path);parser.add_argument('output',type=Path);parser.add_argument('--mode',choices=['uniform_preparation','weak_heterogeneity'],required=True);args=parser.parse_args()
    started=time.perf_counter();case=json.loads((ROOT/'environment/data/case.json').read_text());data=json.loads((ROOT/'environment/data/measurements.json').read_text());source=json.loads(args.point.read_text());point=source.get('result',source)
    altered=copy.deepcopy(data);z0=np.array([point['parameters'][k] for k in KEYS])/SCALE;bounds=np.array([case['parameter_bounds'][k] for k in KEYS])/SCALE[:,None]
    if args.mode=='uniform_preparation':
        for run in altered['runs']:run['initial'].pop('preparation_temperature_K')
        free=np.arange(len(KEYS));fixed={}
    else:
        fixed=dict(fraction_B=.15,delta_E_des_B=.3,delta_E_ox_B=0.,delta_E_hop_B=0.)
        for k,v in fixed.items():z0[KEYS.index(k)]=v/SCALE[KEYS.index(k)]
        free=np.arange(15)
    model=Objective(case,altered,verbose=True)
    def unpack(v):
        z=z0.copy();z[free]=v;return z
    def fun(v):return model.evaluate(unpack(v))[0]
    def jac(v):return model.evaluate(unpack(v))[1][:,free]
    result=least_squares(fun,z0[free],jac=jac,bounds=(bounds[free,0],bounds[free,1]),ftol=2e-7,xtol=2e-7,gtol=1e-5,max_nfev=100,x_scale='jac')
    z=unpack(result.x);r,_,n=model.evaluate(z)
    original=Objective(case,data,verbose=True);rr,_,nn=original.evaluate(z)
    # For omitted preparation weighting the approximation's run coefficients
    # are retained in the complete artifact; also report optimal true-model
    # run refitting separately, to give the approximation its best chance.
    sys.path.insert(0,str(ROOT/'tests'));from experiment_check import deviance
    p=dict(zip(KEYS,(z*SCALE).tolist()));true_D=deviance(p,n,data)[0]
    report=dict(provenance=__doc__,mode=args.mode,fixed_coefficients=fixed,start=point['parameters'],input_sha256=hashlib.sha256(args.point.read_bytes()).hexdigest(),optimizer=dict(success=bool(result.success),message=str(result.message),evaluations=model.calls),approximation_deviance=float(r@r),full_model_deviance_retained_nuisance=float(true_D),full_model_deviance_refitted_runs=float(rr@rr),parameters=p,run_nuisance=n,deviance=float(true_D),elapsed_s=time.perf_counter()-started)
    args.output.write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps({k:report[k] for k in ('mode','approximation_deviance','full_model_deviance_retained_nuisance','full_model_deviance_refitted_runs','elapsed_s')}),flush=True)
if __name__=='__main__':main()
