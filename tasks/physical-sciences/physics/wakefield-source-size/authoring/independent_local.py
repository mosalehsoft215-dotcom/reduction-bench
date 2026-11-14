"""Physical-coordinate quadrature action and finite-difference Jordan coefficients."""
import json,time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from independent_lattice import Independent
B=Path(__file__).resolve().parent.parent;E=B/'authoring/evidence';t=time.monotonic()
ind=json.loads((E/'independent_lattice.json').read_text());r=json.loads((E/'reference_output.json').read_text())
case=json.loads((B/'environment/data/lattice_case.json').read_text());q=np.array(ind['parameters']);M=Independent(q,case)
records=[]
for query,a,b in zip(case['orbit_queries'],ind['orbit_sets'],r['lattice']['orbit_sets']):
    for o,truth in zip(a['cycles'],b['cycles']):
        action=0.;g=query['strength']
        for x,u in o['states']:
            for j,d,k in [(0,M.L,g),(1,M.H,M.r*g)]:
                y=x-q[4*j+1];c,f=q[4*j+2:4*j+4]
                fun=lambda y:y+c*y**3/M.R**2+f*y**5/M.R**4
                u-=k*q[4*j]*fun(y)
                action+=d*u*u/2-k*q[4*j]*quad(fun,0,y,epsabs=1e-14,epsrel=1e-12)[0];x+=d*u
        _,J=M.iterate(np.array(o['states'][0]),g,query['period'],True);eig=np.linalg.eigvals(J);eig=sorted(eig,key=lambda e:(e.real,e.imag))
        ee=np.array([[e.real,e.imag] for e in eig]);records.append(dict(id=query['id'],action_absolute_difference=abs(action-truth['action_mm_mrad']),multiplier_max_absolute_difference=float(np.max(abs(ee-truth['multipliers'])))))
z=np.array(ind['fold'][:2]);g=ind['fold'][2];S=np.diag([M.R,M.U]);Si=np.linalg.inv(S)
def G(v,h):return Si@M.iterate(S@v,h,3)
zn=Si@z;rows=[]
for h in [5e-4,2e-4,4e-5]:
    J=np.column_stack([(G(zn+h*np.eye(2)[i],g)-G(zn-h*np.eye(2)[i],g))/(2*h) for i in range(2)])
    _,_,V=np.linalg.svd(J-np.eye(2));v=V[-1];v*=1 if v[0]>0 else -1;w=np.array([-v[1],v[0]])
    row=dict(step=h,critical_direction=v.tolist(),shear=float(v@(J-np.eye(2))@w),parameter_coefficient=float(w@(G(zn,g+h)-G(zn,g-h))/(2*h)),
        quadratic_coefficient=float(w@(G(zn+h*v,g)-2*G(zn,g)+G(zn-h*v,g))/(2*h*h)))
    truth=r['lattice']['saddle_center'];row['relative_differences']={k:abs(row[k]/truth[k]-1) for k in ['shear','parameter_coefficient','quadratic_coefficient']};rows.append(row)
out=dict(orbit_checks=records,jordan_checks=rows,runtime_sec=time.monotonic()-t)
(E/'independent_local.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
