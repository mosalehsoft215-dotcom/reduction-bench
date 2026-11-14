"""Criticality change, periodic saddle node, and two coexisting cycles."""
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.integrate import solve_bvp,solve_ivp
from model import *
from analyze import hopfs

def criticality(p,s,m,interval):
    def point(u):
        hs=hopfs(p,s,m,u)
        assert len(hs)==1
        return hs[0]
    u=brentq(lambda u:point(u)['first_lyapunov'],*interval,xtol=1e-12)
    h=point(u);du=1e-6
    slope=(point(u+du)['first_lyapunov']-point(u-du)['first_lyapunov'])/(2*du)
    return {'command_V':u,'state':h['state'],'total_capacitance_F_cm2':h['total_capacitance_F_cm2'],
            'frequency_Hz':h['frequency_Hz'],'lyapunov_slope_per_V':slope}

class Branch:
    def __init__(self,p,s,m,u):
        self.p=p;self.s=s;self.m=m;self.u=u
        hs=hopfs(p,s,m,u);assert len(hs)==1
        self.h=hs[0];self.eq=np.array(self.h['state'])
        self.cap=self.h['total_capacitance_F_cm2']
        _,self.w,self.q=lyapunov(self.eq,p,s,m,self.cap)
        self.q*=np.exp(1j*(np.pi/2-np.angle(self.q[1])))
        if self.q[0].real<0:self.q=-self.q
        self.cache={}
    def solve(self,d):
        if d in self.cache:return self.cache[d]
        if not self.cache and d>.005:self.solve(.005)
        if self.cache:
            closest=min(self.cache,key=lambda v:abs(v-d))
            if abs(closest-d)>.006:
                for v in np.linspace(closest,d,int(abs(closest-d)/.005)+2)[1:-1]:self.solve(float(v))
                closest=min(self.cache,key=lambda v:abs(v-d))
            previous=self.cache[closest]
        mesh=np.linspace(0,1,241)
        if self.cache:
            y=self.eq[:,None]+(previous.sol(mesh)-self.eq[:,None])*d/closest
            pars=previous.p
        else:
            y=self.eq[:,None]+d/self.q[0].real*np.real(self.q[:,None]*np.exp(2j*np.pi*mesh))
            pars=np.log([2*np.pi/self.w,self.cap])
        def fun(t,y,z):return np.exp(z[0])*rhs(y,self.p,self.s,self.m,self.u,np.exp(z[1]))
        def bc(a,b,z):return np.r_[a-b,a[0]-self.eq[0]-d,a[1]-self.eq[1]]
        sol=solve_bvp(fun,bc,mesh,y,p=pars,tol=3e-10,max_nodes=20000)
        if not sol.success:raise RuntimeError(sol.message)
        self.cache[d]=sol
        return sol
    def cap_at(self,d):return float(np.exp(self.solve(d).p[1]))
    def artifact(self,d):
        sol=self.solve(d);T,C=np.exp(sol.p)
        def var(t,z):return (jacobian(sol.sol(t/T),self.p,self.s,self.m,C)@z.reshape(3,3)).ravel()
        ans=solve_ivp(var,[0,T],np.eye(3).ravel(),method='DOP853',rtol=2e-12,atol=2e-14)
        mon=ans.y[:,-1].reshape(3,3);vel=rhs(sol.sol(0),self.p,self.s,self.m,self.u,C)
        section=(np.eye(3)-np.outer(vel,np.array([0.,1.,0.]))/vel[1])@mon
        mu=np.linalg.eigvals(section[np.ix_([0,2],[0,2])])
        assert np.max(abs(mu.imag))<1e-8
        return {'period_s':float(T),'total_capacitance_F_cm2':float(C),
                'states':sol.sol(np.arange(64)/64).T.tolist(),
                'nontrivial_floquet_multipliers':sorted(mu.real.tolist())}

def bistability(p,s,m,u,interval):
    b=Branch(p,s,m,u)
    # Capacitance minimum along the section-amplitude continuation.
    opt=minimize_scalar(b.cap_at,bounds=interval,method='bounded',options={'xatol':2e-9})
    assert opt.success
    df=float(opt.x);cf=b.cap_at(df);mid=(cf+b.cap)/2
    small=brentq(lambda d:b.cap_at(d)-mid,.005,df,xtol=2e-10)
    large=brentq(lambda d:b.cap_at(d)-mid,df,interval[1],xtol=2e-10)
    return {'equilibrium':b.eq.tolist(),'hopf_capacitance_F_cm2':b.cap,
            'coexistence_capacitance_F_cm2':mid,'fold_orbit':b.artifact(df),
            'unstable_orbit':b.artifact(small),'stable_orbit':b.artifact(large)}
