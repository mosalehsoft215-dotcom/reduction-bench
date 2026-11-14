"""Independent physical derivation and generic counting accuracy evidence."""
from pathlib import Path
import sys, os, json, time, hashlib
import numpy as np

R=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(R/'solution'),str(R/'tests'),str(R/'environment/data')]
os.environ['RAMAN_DATA']=str(R/'environment/data')
from physical_model import derive
from physical_model import predictions
import optics
from independent_dynamics import independent_rates
from count_engine import comparisons
from independent_contour_kernel import count_comparator_kernel

start=time.perf_counter()
catalog,lam,Q,limits=derive()
derived_s=time.perf_counter()-start
catalog2,lam2,Q2,cfg=independent_rates()
assert np.array_equal(catalog,catalog2)
mask=np.isfinite(lam)
assert np.array_equal(mask,np.isfinite(lam2))
rate_difference=float(np.max(abs(lam[mask]-lam2[mask])/np.maximum(1.,abs(lam2[mask]))))
generator_difference=float(np.max(abs(Q-Q2)))
assert rate_difference<1e-11 and generator_difference<1e-12
checks=[]
for setting,state,channel,dwell,threshold in ((0,0,0,4,25),(14,87,2,8,60),(27,17,3,4,35),(40,11,1,12,100),(53,142,0,20,150)):
    a,b=comparisons(Q[setting,state:state+1],lam[setting,state:state+1,:,channel],dwell*cfg['exposure_unit_s'],[threshold])
    c,d=count_comparator_kernel(Q2[setting,state],lam2[setting,state,:,channel],dwell*cfg['exposure_unit_s'],threshold)
    discrepancy=float(max(np.max(abs(a[0,0]-c)),np.max(abs(b[0,0]-d))))
    assert discrepancy<1e-10
    checks.append(dict(setting_id=setting,initial_state=state,channel=channel,dwell_units=dwell,threshold=threshold,maximum_absolute_difference=discrepancy))
start=time.perf_counter()
comparisons(Q[27],lam[27,:,:,3],2.8,list(range(250)))
batch_s=time.perf_counter()-start
saved_tensors=optics.T
try:
    optics.T={w:t.real.astype(complex) for w,t in saved_tensors.items()}
    _,real_only_rates,real_only_Q,_=derive()
finally:
    optics.T=saved_tensors
correct_rows=predictions(catalog,lam,Q)
real_only_rows=predictions(catalog,real_only_rates,real_only_Q)
phase_control=[]
for correct,wrong in zip(correct_rows,real_only_rows):
    truth=np.asarray(correct['photon_rates_per_s']);approx=np.asarray(wrong['photon_rates_per_s'])
    allowance=1e-8+1e-7*abs(truth)
    phase_control.append(dict(query_id=correct['query_id'],entries_outside_photon_reporting_tolerance=int(np.sum(abs(truth-approx)>allowance)),maximum_absolute_difference=float(np.max(abs(truth-approx)))))
assert sum(row['entries_outside_photon_reporting_tolerance'] for row in phase_control)>0
public=R/'environment/data'
assert not list(public.glob('*.npz'))
assert {p.name for p in public.glob('*.py')}=={'count_engine.py'}
result=dict(raw_physical_derivation_seconds=derived_s,scaled_photon_rate_difference=rate_difference,generator_absolute_difference=generator_difference,kernel_checks=checks,generic_250_threshold_batch_seconds=batch_s,public_python_files=['count_engine.py'],public_derived_response_tables=0,real_tensor_only_scientific_control=phase_control)
(R/'authoring/evidence/physical_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
