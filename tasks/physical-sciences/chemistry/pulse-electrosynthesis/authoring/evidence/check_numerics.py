import sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'solution'),str(ROOT/'tests')]
from fast_model import trajectory as split
from physics import trajectory as bdf
from design import responses,information
c=json.load((ROOT/'environment/data/cell.json').open());s=json.load((ROOT/'environment/data/scenarios.json').open())
p=np.array([s['nominal'][k] for k in c['scenario_parameter_order']])
results={}
a=np.array([[6.,-1.15]]);b=np.array([[.07,-1.15],[2.13,-1.15],[3.8,-1.15]])
m1,q1=bdf(a,p,n=32);m2,q2=bdf(b,p,n=32)
results['equivalent_subdivision_outcomes']=float(abs(m1-m2).max())
results['equivalent_subdivision_bins']=float(abs(q1-q2).max())
assert abs(m1-m2).max()<1e-6 and abs(q1-q2).max()<1e-7
zero=p.copy();zero[6:8]=0
for label,fun in [('split',split),('bdf',bdf)]:
    m,q=fun(a,zero,n=32)
    assert np.max(abs(m[:3]))<1e-12 and abs(m[3]-1)<1e-12 and np.max(abs(q))<1e-12
    results[label+'_zero_flux']=True
programs=[np.array([[.24,-1.25],[.06,-1.04]]*20),np.array([[.20,-1.23],[.05,-1.055]]*24)]
m,q,j=responses(programs,p,40,.002)
h=information(q,j,c);hswap=information(q[::-1],j[::-1],c)
results['experiment_exchange_information_error']=float(abs(h-hswap).max())
assert abs(h-hswap).max()<1e-8
rho=c['charge_bin_noise_rho'];indices=np.arange(120)
v=np.zeros((122,122));v[:120,:120]=c['charge_bin_noise_sd']**2*rho**abs(indices[:,None]-indices[None,:]);v[120:,120:]=c['assay_noise_sd']**2*np.eye(2)
f=np.zeros((12,12))
for run in range(2):
    raw=np.zeros((122,12));raw[:,:8]=j[run];raw[:120,8+2*run]=q[run];raw[:120,9+2*run]=1
    f+=raw.T@np.linalg.solve(v,raw)
sd=np.array(c['nuisance_prior_sd']+[c['gain_prior_sd'],c['offset_prior_sd']]*2)
f[2:,2:]+=np.diag(1/sd**2)
dense=f[:2,:2]-f[:2,2:]@np.linalg.solve(f[2:,2:],f[2:,:2])
results['dense_covariance_vs_whitened_projection']=float(abs(h-dense).max())
assert abs(h-dense).max()<1e-8
_,_,j2=responses(programs,p,40,.002,h=.0015)
results['halved_log_sensitivity_step_max_difference']=float(abs(j-j2).max())
assert abs(j-j2).max()<1e-6
results['mass_residual']=float(abs(m[:,[0,1,3,4]].sum(axis=1)-1).max())
results['charge_bin_sum_residual']=float(abs(q.sum(axis=1)-m[:,2]).max())
assert results['mass_residual']<1e-8 and results['charge_bin_sum_residual']<1e-10
(Path(__file__).parent/'numerical_invariants.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2),flush=True)
