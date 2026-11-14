import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_bvp,solve_ivp
from model import *

def crossings(fun,grid):
    y=[fun(t) for t in grid];roots=[]
    for a,b,fa,fb in zip(grid[:-1],grid[1:],y[:-1],y[1:]):
        if fa*fb<0:roots.append(brentq(fun,a,b,xtol=2e-13))
    return roots

def fold_slope(th,p,cond,m):
    _,x,_=steady_theta(th,p,cond,m)
    _,k,a,g,kc,G,km,R,_=p;phi,_,c=x;eta=phi-g*(2*th-1)
    ef=np.exp(-a*eta);er=np.exp((1-a)*eta)
    catprime=kc*cond['bulk_M']*m*th**(m-1)*c*c
    vp=-k*(a*(1-th)*ef+(1-a)*th*er)
    etap=(m*catprime+k*(ef+er))/vp
    return (etap+2*g)/FRT-R*F*G*m*catprime

def folds(p,cond,m):
    roots=crossings(lambda t:fold_slope(t,p,cond,m),np.linspace(.005,.995,601))
    return [{'theta':t,'command_V':steady_theta(t,p,cond,m)[0]} for t in roots]

def equilibria(p,cond,m,u):
    cuts=[.000001]+[s['theta'] for s in folds(p,cond,m)]+[.999999]
    return [steady_theta(t,p,cond,m)[1] for t in crossings(lambda t:steady_theta(t,p,cond,m)[0]-u,cuts)]

def hopfs(p,cond,m,u):
    out=[]
    for x in equilibria(p,cond,m,u):
        _,hs=hopf_at_theta(x[1],p,cond,m)
        for cap,freq in hs:
            if not .0005<=cap<=.03:continue
            l1,w,q=lyapunov(x,p,cond,m,cap)
            def eigpart(c):
                es=np.linalg.eigvals(jacobian(x,p,cond,m,c))
                return es[np.argmax(es.imag)].real
            slope=(eigpart(cap*1.00001)-eigpart(cap*.99999))/(cap*.00002)
            out.append({'state':x.tolist(),'total_capacitance_F_cm2':cap,
                        'frequency_Hz':freq,'first_lyapunov':float(l1),
                        'eigenvalue_slope_per_F_cm2_s':float(slope)})
    return sorted(out,key=lambda a:(a['state'][1],a['total_capacitance_F_cm2']))

def orbit(p,cond,m,u,hopf,delta=.05,samples=64):
    xeq=np.array(hopf['state']);cap=hopf['total_capacitance_F_cm2']
    l1,omega,q=lyapunov(xeq,p,cond,m,cap)
    q*=np.exp(1j*(np.pi/2-np.angle(q[1])))
    if q[0].real<0:q=-q
    previous=None
    for d in np.linspace(.005,delta,10):
        mesh=np.linspace(0,1,181)
        if previous is None:
            y=xeq[:,None]+d/q[0].real*np.real(q[:,None]*np.exp(2j*np.pi*mesh))
            pars=np.log([2*np.pi/omega,cap])
        else:
            y=xeq[:,None]+(previous.sol(mesh)-xeq[:,None])*d/old_delta
            pars=previous.p
        def fun(t,y,z):return np.exp(z[0])*rhs(y,p,cond,m,u,np.exp(z[1]))
        def bc(a,b,z):return np.r_[a-b,a[0]-xeq[0]-d,a[1]-xeq[1]]
        sol=solve_bvp(fun,bc,mesh,y,p=pars,tol=1e-9,max_nodes=15000)
        if not sol.success:raise RuntimeError(sol.message)
        previous=sol;old_delta=d
    T,C=np.exp(sol.p)
    def variational(t,flat):return (jacobian(sol.sol(t/T),p,cond,m,C)@flat.reshape(3,3)).ravel()
    result=solve_ivp(variational,[0,T],np.eye(3).ravel(),method='DOP853',rtol=2e-11,atol=1e-13)
    mu=np.linalg.eigvals(result.y[:,-1].reshape(3,3))
    assert max(abs(mu.imag))<1e-7
    mu=mu.real;mu=np.delete(mu,np.argmin(abs(mu-1)))
    return {'period_s':float(T),'total_capacitance_F_cm2':float(C),
            'states':sol.sol(np.arange(samples)/samples).T.tolist(),
            'nontrivial_floquet_multipliers':sorted(mu.tolist())}
