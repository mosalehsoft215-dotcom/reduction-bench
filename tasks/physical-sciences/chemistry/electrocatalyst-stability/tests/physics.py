"""Independent verification physics: complex-step and Fourier differentiation."""
from itertools import product
from math import factorial
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

F=96485.33212;BETA=38.9217445

def field(y,p,s,m,u,cap):
    phi,theta,c=y
    eta=phi-p['interaction']*(2*theta-1)
    reduction=p['k0_s_inv']*((1-theta)*np.exp(-p['alpha']*eta)-theta*np.exp((1-p['alpha'])*eta))
    chemical=p['kcat_M_inv_s_inv']*theta**m*c
    return np.array([(BETA*(u-p['E0_V'])-phi)/p['Ru_ohm_cm2']+BETA*F*p['Gamma_mol_cm2']*reduction,
          reduction-m*s['bulk_M']*chemical,
          p['km_s_inv']*np.sqrt(s['rpm']/1000)*(1-c)-1000*p['Gamma_mol_cm2']/s['film_cm']*chemical])/np.array([cap,1,1]).reshape((3,)+(1,)*(np.ndim(y)-1))

def derivative(y,p,s,m,u,cap):
    return np.column_stack([field(y.astype(complex)+1e-25j*np.eye(3)[j],p,s,m,u,cap).imag/1e-25 for j in range(3)])

def state_at_coverage(t,p,s,m):
    supply=p['km_s_inv']*np.sqrt(s['rpm']/1000)
    c=supply/(supply+1000*p['Gamma_mol_cm2']/s['film_cm']*p['kcat_M_inv_s_inv']*t**m)
    reaction=p['kcat_M_inv_s_inv']*s['bulk_M']*t**m*c
    def error(phi):
        eta=phi-p['interaction']*(2*t-1)
        return p['k0_s_inv']*((1-t)*np.exp(-p['alpha']*eta)-t*np.exp((1-p['alpha'])*eta))-m*reaction
    phi=brentq(error,-80,80,xtol=1e-13)
    current=F*p['Gamma_mol_cm2']*m*reaction
    u=p['E0_V']+phi/BETA-p['Ru_ohm_cm2']*current
    return u,np.array([phi,t,c]),current

def roots_on_grid(fn,grid):
    values=np.array([fn(float(v)) for v in grid]);out=[]
    for k in range(len(grid)-1):
        if values[k]==0:out.append(float(grid[k]))
        if values[k]*values[k+1]<0:out.append(float(brentq(fn,grid[k],grid[k+1],xtol=1e-13)))
    if values[-1]==0:out.append(float(grid[-1]))
    return sorted(set(out))

def folds(p,s,m):
    def slope(t):
        h=1e-5
        return sum(weight*state_at_coverage(t+offset*h,p,s,m)[0] for offset,weight in [(-2,1),(-1,-8),(1,8),(2,-1)])/(12*h)
    points=roots_on_grid(slope,np.linspace(.005,.995,1001))
    return [{'theta':t,'command_V':state_at_coverage(t,p,s,m)[0]} for t in points]

def equilibria(p,s,m,u):
    ts=roots_on_grid(lambda t:state_at_coverage(t,p,s,m)[0]-u,np.r_[np.geomspace(1e-8,.005,100),np.linspace(.005,.995,1201),1-np.geomspace(.005,1e-8,100)])
    return [state_at_coverage(t,p,s,m)[1] for t in ts]

def spectrum(y,p,s,m,u,freq):
    R=p['Ru_ohm_cm2'];cap=p['Cdl_F_cm2'];a=derivative(y,p,s,m,u,cap)
    return np.array([R/(1-np.linalg.solve(2j*np.pi*f*np.eye(3)-a,np.array([1/(R*cap),0,0]))[0]) for f in freq])

def hopf_points(p,s,m,u,interval):
    result=[]
    for x in equilibria(p,s,m,u):
        def coeff(cap):
            a=derivative(x,p,s,m,u,cap)
            return np.array([-np.trace(a),sum(a[i,i]*a[j,j]-a[i,j]*a[j,i] for i,j in [(0,1),(0,2),(1,2)]),-np.linalg.det(a)])
        def residual(cap):
            a,b,c=coeff(cap);return a*b-c
        for cap in roots_on_grid(residual,np.geomspace(*interval,801)):
            a,b,c=coeff(cap)
            if min(a,b,c)>0:
                result.append((x,cap,np.sqrt(b)/(2*np.pi)))
    return sorted(result,key=lambda a:(a[0][1],a[1]))

def fourier_tensors(y,p,s,m,u,cap,radius=.035,n=16):
    angles=np.exp(2j*np.pi*np.arange(n)/n)
    z=np.array(np.broadcast_arrays(y[0]+radius*angles[:,None,None],y[1]+radius*angles[None,:,None],y[2]+radius*angles[None,None,:]))
    coeff=np.fft.fftn(field(z,p,s,m,u,cap),axes=(1,2,3))/n**3
    tensors=[]
    for order in [2,3]:
        a=np.zeros((3,)*(order+1),complex)
        for ix in product(range(3),repeat=order):
            powers=tuple(ix.count(j) for j in range(3))
            a[(slice(None),)+ix]=coeff[(slice(None),)+powers]*np.prod([factorial(k) for k in powers])/radius**order
        assert np.max(abs(a.imag))<1e-6
        tensors.append(a.real)
    return tensors

def normal_form(y,p,s,m,u,cap):
    a=derivative(y,p,s,m,u,cap);ev,v=np.linalg.eig(a)
    k=np.argmax(ev.imag);w=ev[k].imag;q=v[:,k]/np.linalg.norm(v[:,k])
    # Left eigenvector from the inverse right-eigenvector matrix (independent normalization).
    left=np.linalg.inv(v)[k,:]*np.linalg.norm(v[:,k])
    b,c=fourier_tensors(y,p,s,m,u,cap)
    B=lambda v,w:np.einsum('ijk,j,k->i',b,v,w)
    C=lambda v,w,z:np.einsum('ijkl,j,k,l->i',c,v,w,z)
    h20=np.linalg.solve(2j*w*np.eye(3)-a,B(q,q))
    h11=np.linalg.solve(a,B(q,q.conj()))
    ell=float((left@(C(q,q,q.conj())+B(q.conj(),h20)-2*B(q,h11))).real/(2*w))
    h=cap*2e-5
    def realpart(c):
        vals=np.linalg.eigvals(derivative(y,p,s,m,u,c));return vals[np.argmax(vals.imag)].real
    slope=(realpart(cap+h)-realpart(cap-h))/(2*h)
    return ell,slope

def integrate_orbit(start,T,cap,p,s,m,u,samples=64):
    def extended(t,y):
        x=y[:3];phi=y[3:].reshape(3,3)
        return np.r_[field(x,p,s,m,u,cap),(derivative(x,p,s,m,u,cap)@phi).ravel()]
    def section(t,y):return y[1]-start[1]
    section.direction=-1
    ans=solve_ivp(extended,[0,T],np.r_[start,np.eye(3).ravel()],method='DOP853',
                  rtol=1e-10,atol=2e-12,dense_output=True,events=section,max_step=T/150)
    if not ans.success:raise ValueError('Orbit integration failed')
    monodromy=ans.y[3:,-1].reshape(3,3)
    multipliers=np.linalg.eigvals(monodromy)
    ix=np.argmin(abs(multipliers-1));unit=multipliers[ix]
    velocity=field(start,p,s,m,u,cap)
    projected=(np.eye(3)-np.outer(velocity,np.array([0.,1.,0.]))/velocity[1])@monodromy
    others=np.linalg.eigvals(projected[np.ix_([0,2],[0,2])])
    if max(abs(others.imag))>1e-5:raise ValueError('Non-real Floquet multipliers')
    internal_crossings=[t for t in ans.t_events[0] if T*1e-4<t<T*(1-1e-4)]
    return {'closure':ans.y[:3,-1]-start,'states':ans.sol(np.arange(samples)*T/samples)[:3].T,
            'dense_states':ans.sol(np.linspace(0,T,1001))[:3].T,
            'unit':float(unit.real),'multipliers':np.sort(others.real),'extra_crossings':len(internal_crossings)}
