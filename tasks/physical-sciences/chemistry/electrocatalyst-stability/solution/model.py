"""Dimensionless-potential model of a cooperative redox catalyst film."""
import numpy as np
from scipy.optimize import brentq
from math import factorial
from itertools import product

F=96485.33212
FRT=38.9217445
NAMES=['E0_V','k0_s_inv','alpha','interaction','kcat_M_inv_s_inv','Gamma_mol_cm2','km_s_inv','Ru_ohm_cm2','Cdl_F_cm2']

def rates(x,p,condition,m):
    phi,theta,c=x
    E0,k0,alpha,g,kcat,gamma,km,ru,cdl=p
    eta=phi-g*(2*theta-1)
    forward=k0*(1-theta)*np.exp(-alpha*eta)
    reverse=k0*theta*np.exp((1-alpha)*eta)
    cat=kcat*condition['bulk_M']*c*theta**m
    return forward-reverse,cat

def rhs(x,p,condition,m,command,cap=None):
    E0,k0,alpha,g,kcat,gamma,km,ru,cdl=p
    cap=cdl if cap is None else cap
    v,cat=rates(x,p,condition,m)
    transfer=km*np.sqrt(condition['rpm']/1000.)
    h=1000*gamma/condition['film_cm']
    return np.array([(FRT*(command-E0)-x[0])/(ru*cap)+FRT*F*gamma*v/cap,
                     v-m*cat,transfer*(1-x[2])-h*kcat*x[1]**m*x[2]])

def jacobian(x,p,condition,m,cap=None):
    E0,k0,alpha,g,kcat,gamma,km,ru,cdl=p
    cap=cdl if cap is None else cap
    phi,th,c=x;eta=phi-g*(2*th-1)
    ef=np.exp(-alpha*eta);er=np.exp((1-alpha)*eta)
    vp=k0*(-alpha*(1-th)*ef-(1-alpha)*th*er)
    vt=k0*((-1+2*alpha*g*(1-th))*ef-(1-2*(1-alpha)*g*th)*er)
    ct=kcat*condition['bulk_M']*c*m*th**(m-1)
    cc=kcat*condition['bulk_M']*th**m
    trans=km*np.sqrt(condition['rpm']/1000.);h=1000*gamma/condition['film_cm']
    a=np.array([[-1/(ru*cap)+FRT*F*gamma*vp/cap,FRT*F*gamma*vt/cap,0],
                [vp,vt-m*ct,-m*cc],
                [0,-h*kcat*c*m*th**(m-1),-trans-h*kcat*th**m]])
    return a

def steady_theta(th,p,condition,m):
    E0,k0,alpha,g,kcat,gamma,km,ru,cdl=p
    trans=km*np.sqrt(condition['rpm']/1000.);h=1000*gamma/condition['film_cm']
    c=trans/(trans+h*kcat*th**m)
    cat=kcat*condition['bulk_M']*c*th**m
    eta=brentq(lambda v:k0*((1-th)*np.exp(-alpha*v)-th*np.exp((1-alpha)*v))-m*cat,-100,100,xtol=1e-13)
    phi=eta+g*(2*th-1)
    e=E0+phi/FRT;command=e-ru*F*gamma*m*cat
    return command,np.array([phi,th,c]),F*gamma*m*cat

def impedance(x,p,condition,m,frequencies):
    a=jacobian(x,p,condition,m);ru=p[7];cap=p[8]
    b=np.array([FRT/(ru*cap),0.,0.])
    mat=2j*np.pi*np.asarray(frequencies)[:,None,None]*np.eye(3)-a
    response=np.linalg.solve(mat,np.broadcast_to(b,(len(frequencies),3))[...,None])[...,0]
    return ru/(1-response[:,0]/FRT)

def hopf_at_theta(th,p,condition,m):
    command,x,current=steady_theta(th,p,condition,m)
    # Routh-Hurwitz determinant is quadratic in inverse capacitance.
    vals=[]
    for inv in [0.,1.,2.]:
        a=jacobian(x,p,condition,m,cap=1/max(inv,1e-30))
        coeff=np.poly(a).real;vals.append(coeff[1:])
    b1,b2,_=vals[0];a1,a2,a3=vals[1]-vals[0]
    poly=[a1*a2,a1*b2+b1*a2-a3,b1*b2]
    out=[]
    for inv in np.roots(poly):
        if abs(inv.imag)>1e-8 or inv.real<=0:continue
        cap=1/inv.real;a=jacobian(x,p,condition,m,cap)
        coeff=np.poly(a).real
        if min(coeff[1:])>0 and np.linalg.eigvals(a).real.max()<1e-6:
            out.append((cap,np.sqrt(coeff[2])/(2*np.pi)))
    return command,out

def tensors(x,p,condition,m,cap):
    """Analytic second and third derivatives in (phi, theta, c)."""
    _,k0,alpha,g,kcat,gamma,km,ru,_=p
    phi,th,c=x;eta=phi-g*(2*th-1)
    def deriv(indices):
        ip=indices.count(0);it=indices.count(1);ic=indices.count(2)
        dv=0.
        if ic==0:
            for amp, slope, power in [(1-th,-1.,-alpha),(-th,-1.,1-alpha)]:
                z=-2*g*power
                dv+=k0*np.exp(power*eta)*power**ip*(amp*z**it+(it*slope*z**(it-1) if it else 0))
        dt=0.
        if ip==0 and it<=m and ic<=1:
            dt=kcat*factorial(m)/factorial(m-it)*th**(m-it)*(c if ic==0 else 1.)
        return [FRT*F*gamma*dv/cap,dv-m*condition['bulk_M']*dt,-1000*gamma/condition['film_cm']*dt]
    arrays=[]
    for order in [2,3]:
        a=np.empty((3,)*(order+1))
        for ix in product(range(3),repeat=order):a[(slice(None),)+ix]=deriv(list(ix))
        arrays.append(a)
    return arrays

def lyapunov(x,p,condition,m,cap):
    a=jacobian(x,p,condition,m,cap)
    eig,vec=np.linalg.eig(a);j=np.argmax(eig.imag);omega=eig[j].imag
    q=vec[:,j];q=q/np.linalg.norm(q)
    ev,v=np.linalg.eig(a.T);adj=v[:,np.argmin(abs(ev+1j*omega))]
    adj=adj/np.conj(np.vdot(adj,q))
    b,c=tensors(x,p,condition,m,cap)
    B=lambda u,v:np.einsum('ijk,j,k->i',b,u,v)
    C=lambda u,v,w:np.einsum('ijkl,j,k,l->i',c,u,v,w)
    h20=np.linalg.solve(2j*omega*np.eye(3)-a,B(q,q))
    h11=np.linalg.solve(a,B(q,q.conj()))
    l1=np.vdot(adj,C(q,q,q.conj())-2*B(q,h11)+B(q.conj(),h20)).real/(2*omega)
    return l1,omega,q
