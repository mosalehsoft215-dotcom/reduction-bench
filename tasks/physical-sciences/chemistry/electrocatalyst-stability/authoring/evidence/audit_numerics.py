from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
import json,sys,time
import numpy as np
sys.path.insert(0,str(ROOT/'tests'));sys.path.insert(0,str(ROOT/'solution'))
import physics as independent
import model as reference
from validate import CASE,RECORDS,CONDITIONS
from scipy.optimize import brentq
from fit import encode,decode

artifact=json.load(open(ROOT/'authoring/evidence/reference-output.json'));p=artifact['parameters'];m=artifact['mechanism']
pv=np.array([p[k] for k in CASE['parameter_order']])
report={}
derivative_errors=[];tensor_errors=[];hopf=[]
for r,q in zip(artifact['queries'],CASE['hopf_queries']):
    s=CONDITIONS[q['condition_id']];u=q['command_V']
    for h in r['hopf_points']:
        x=np.array(h['state']);cap=h['total_capacitance_F_cm2']
        a=reference.jacobian(x,pv,s,m,cap);b=independent.derivative(x,p,s,m,u,cap)
        derivative_errors.append(float(np.max(abs(a-b))))
        bb,cc=reference.tensors(x,pv,s,m,cap)
        for radius,n in [(.025,16),(.05,20)]:
            b2,c2=independent.fourier_tensors(x,p,s,m,u,cap,radius=radius,n=n)
            tensor_errors.append({'radius':radius,'points_per_axis':n,'second_max_abs':float(np.max(abs(bb-b2))),'third_max_abs':float(np.max(abs(cc-c2)))})
        eig=np.linalg.eigvals(b)
        hopf.append({'id':q['id'],'theta':float(x[1]),'C':cap,'frequency_Hz':h['frequency_Hz'],'l1':h['first_lyapunov'],'eigenvalues':[[float(z.real),float(z.imag)] for z in eig]})
report['jacobian_max_abs_errors']=derivative_errors
report['tensor_comparison']=tensor_errors;report['hopf']=hopf

# Whitened sensitivity matrix after eliminating equilibrium constraints exactly.
base_theta={r['id']:a['state'][1] for r,a in zip(RECORDS,artifact['measurements'])}
def predictions(v):
    pp=dict(zip(CASE['parameter_order'],decode(v)));out=[]
    for r in RECORDS:
        s=CONDITIONS[r['condition_id']];t0=base_theta[r['id']];u=r['command_V']
        t=brentq(lambda th:independent.state_at_coverage(th,pp,s,m)[0]-u,max(1e-7,t0-.005),min(.999999,t0+.005),xtol=1e-14)
        _,x,current=independent.state_at_coverage(t,pp,s,m)
        z=independent.spectrum(x,pp,s,m,u,[a['frequency_Hz'] for a in r['spectrum']])
        out.append(current/r['current_sigma_A_cm2'])
        for val,d in zip(z,r['spectrum']):out.extend([val.real/d['sigma_ohm_cm2'],val.imag/d['sigma_ohm_cm2']])
    return np.array(out)
v=encode(pv);cols=[]
for j in range(9):
    d=np.zeros(9);d[j]=1e-5;cols.append((predictions(v+d)-predictions(v-d))/(2e-5))
J=np.column_stack(cols);sv=np.linalg.svd(J,compute_uv=False);cov=np.linalg.inv(J.T@J)
truth=json.load(open(ROOT/'authoring/provenance/generation_truth.json'))['parameters']
report['local_identifiability']={'singular_values':sv.tolist(),'condition_number':float(sv[0]/sv[-1]),'coordinate_standard_errors':np.sqrt(np.diag(cov)).tolist(),'fitted_parameters':p,'generation_parameters':truth,'coordinate_errors_in_standard_deviations':((v-encode(np.array([truth[k] for k in CASE['parameter_order']])))/np.sqrt(np.diag(cov))).tolist()}
dcJ=J[::51];report['dc_only_singular_values']=np.linalg.svd(dcJ,compute_uv=False).tolist()
report['orbits']=[{'id':q['id'],'period':r['orbits'][0]['period_s'],'C':r['orbits'][0]['total_capacitance_F_cm2'],'multipliers':r['orbits'][0]['nontrivial_floquet_multipliers']} for r,q in zip(artifact['queries'],CASE['hopf_queries'])]
report['chi_square']=artifact['fit_chi_square']
json.dump(report,open(ROOT/'authoring/evidence/science-checks.json','w'),indent=2)
print(json.dumps(report,indent=2))
