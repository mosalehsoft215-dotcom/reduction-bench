"""Independent physical-coordinate shooting, orbit frequencies and manifold seeds.

No import of the resonance reference or reference artifacts. The calibration
operator is separately integrated with finite-range Legendre quadrature.
"""
import json,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import root,least_squares,brentq
from scipy.special import ndtr
from numba import njit
import mpmath as mp

NAMES=['scale_1','center_1_mm','cubic_1','quintic_1','scale_2','center_2_mm','cubic_2','quintic_2']

class Independent:
    def __init__(self,q,case):
        self.q=np.array(q);self.case=case;self.R=case['radius_mm'];self.L=case['between_lenses_m'];self.H=case['return_drift_m'];self.U=self.R/self.L;self.r=case['strength_ratio']
    def map(self,z,g,jac=False):
        x,u=z;A=np.eye(2)
        for j,L,k in [(0,self.L,g),(1,self.H,self.r*g)]:
            y=x-self.q[4*j+1];c,d=self.q[4*j+2:4*j+4]
            u-=k*self.q[4*j]*(y+c*y**3/self.R**2+d*y**5/self.R**4)
            s=-k*self.q[4*j]*(1+3*c*y*y/self.R**2+5*d*y**4/self.R**4)
            A=np.array([[1+L*s,L],[s,1.]])@A;x+=L*u
        return (np.array([x,u]),A) if jac else np.array([x,u])
    def iterate(self,z,g,n,jac=False):
        A=np.eye(2)
        for _ in range(n):z,J=self.map(z,g,True);A=J@A
        return (z,A) if jac else z
    def physical(self,z,g,n):
        for _ in range(n):
            box=np.array(self.case.get('survey_box',[[-1,1],[-1.35,1.35]]));phase=z/[self.R,self.U]
            if np.any(phase<box[:,0]-1e-10) or np.any(phase>box[:,1]+1e-10):return False
            x,u=z;y=x-self.q[1]
            u-=g*self.q[0]*(y+self.q[2]*y**3/self.R**2+self.q[3]*y**5/self.R**4)
            if max(abs(x),abs(x+self.L*u))>self.R+1e-10:return False
            z=self.map(z,g)
        return True
    def center(self,g):return root(lambda z:self.map(z,g)-z,[0,0],tol=1e-12).x
    def mode(self,z,g):
        _,A=self.map(z,g,True);S=np.diag([self.R,self.U]);A=np.linalg.solve(S,A@S)
        vals,V=np.linalg.eig(A);J=np.array([[0,1],[-1,0]])
        for i in range(2):
            v=V[:,i].astype(complex);s=(v.conj()@J@v).imag
            if s>0:
                v/=np.sqrt(s);v*=np.exp(-1j*np.angle(v[0]));p=np.linalg.inv(np.stack([v,v.conj()],axis=1))[0]
                return vals[i],v,p
    def shooting(self,g,n,amplitudes=15,phases=18):
        center=self.center(g);lam,v,p=self.mode(center,g);found=[]
        winding=1 if n==3 else 2
        for amp in np.linspace(.025,.9,amplitudes):
            for phase in np.linspace(0,2*np.pi,phases,endpoint=False):
                ws=amp*np.exp(1j*(phase+2*np.pi*winding*np.arange(n)/n))
                initial=center+2*np.real(ws[:,None]*v)*[self.R,self.U]
                def fun(a):
                    z=a.reshape(n,2);res=np.empty((n,2));J=np.zeros((2*n,2*n))
                    for i in range(n):
                        y,A=self.map(z[i],g,True);j=(i+1)%n;res[i]=y-z[j]
                        J[2*i:2*i+2,2*i:2*i+2]+=A;J[2*i:2*i+2,2*j:2*j+2]-=np.eye(2)
                    return res.ravel(),J
                with np.errstate(over='ignore',invalid='ignore'):
                    f=root(fun,initial.ravel(),jac=True,tol=2e-11)
                z=f.x.reshape(n,2)
                if not np.isfinite(z).all() or np.linalg.norm(f.fun)>1e-8 or not self.physical(z[0],g,n):continue
                if any(n%d==0 and np.linalg.norm(self.iterate(z[0],g,d)-z[0])<1e-6 for d in range(1,n)):continue
                z=np.roll(z,-np.lexsort((z[:,1],z[:,0]))[0],axis=0)
                if any(np.linalg.norm(z[0]-v['states'][0])<1e-6 for v in found):continue
                _,A=self.iterate(z[0],g,n,True)
                found.append(dict(states=z.tolist(),trace=float(np.trace(A))))
        return sorted(found,key=lambda v:tuple(v['states'][0]))
    def fold(self,gfixed=None,start=None,scale1=False):
        box=np.array(self.case['saddle_center']['representative_box']);c=box.mean(axis=1)*[self.R,self.U]
        if gfixed is None:
            initial=np.r_[c,np.mean(self.case['saddle_center']['strength_interval'])] if start is None else start
            def fun(a):
                y,A=self.iterate(a[:2],a[2],3,True);return np.r_[y-a[:2],np.trace(A)-2]
            f=root(fun,initial,tol=1e-11)
        else:
            initial=np.r_[c,self.q[0]] if start is None else start
            def fun(a):
                q=self.q.copy();q[0]=a[2];M=Independent(q,self.case);y,A=M.iterate(a[:2],gfixed,3,True)
                return np.r_[y-a[:2],np.trace(A)-2]
            f=root(fun,initial,tol=1e-11)
        if np.linalg.norm(f.fun)>2e-8:raise RuntimeError('independent fold failed')
        return f.x
    def forms(self):
        lo,hi=self.case['resonance_interval']
        g=brentq(lambda k:np.trace(self.map(self.center(k),k,True)[1])+1,lo,hi,xtol=2e-13)
        z=self.center(g);lam,v,p=self.mode(z,g);vals=[]
        phases=2*np.pi*np.arange(128)/128
        for amp in [.002,.004,.008]:
            w=amp*np.exp(1j*phases);xp=z+2*np.real(w[:,None]*v)*[self.R,self.U]
            wp=np.array([p@((self.map(x,g)-z)/[self.R,self.U]) for x in xp])
            vals.append(np.mean((wp-lam*w)*np.exp(2j*phases))/amp**2)
        X=np.array([.002,.004,.008])**2
        chi=np.polynomial.polynomial.polyfit(X,np.array(vals),2)[0]
        k=self.case['twist_strength'];zc=self.center(k);la,vc,pc=self.mode(zc,k)
        freq=[];actions=[]
        for amp in [.014,.022,.032]:
            initial=zc+2*np.real(amp*vc)*[self.R,self.U]
            states,theta=trajectory(self.q,self.R,self.L,self.H,self.r,k,initial,zc,pc,220000)
            coord=(states-zc)/[self.R,self.U];w=coord@pc
            si=np.argsort(np.angle(w));phasecurve=coord[si]
            area=.5*np.sum(phasecurve[:,0]*np.roll(phasecurve[:,1],-1)-phasecurve[:,1]*np.roll(phasecurve[:,0],-1))
            actions.append(abs(area)/(2*np.pi));freq.append(theta-np.angle(la))
        coeff=np.polynomial.polynomial.polyfit(actions,freq,2)
        return dict(resonance_strength=g,resonance_coefficient=[chi.real,chi.imag],twist_coefficient=coeff[1],frequency_actions=actions,frequency_shifts=freq)

@njit
def trajectory(q,R,L,H,ratio,g,z,center,p,n):
    states=np.empty((n,2));sumphase=0.
    w=p[0]*(z[0]-center[0])/R+p[1]*(z[1]-center[1])/(R/L)
    x,u=z
    for i in range(n):
        for j,dist,k in [(0,L,g),(1,H,ratio*g)]:
            y=x-q[4*j+1];u-=k*q[4*j]*(y+q[4*j+2]*y**3/R**2+q[4*j+3]*y**5/R**4);x+=dist*u
        wp=p[0]*(x-center[0])/R+p[1]*(u-center[1])/(R/L)
        r=wp*np.conj(w);sumphase+=np.arctan2(r.imag,r.real);w=wp
        states[i,0]=x;states[i,1]=u
    return states,sumphase/n

def calibration(q,runs,model):
    t,w=np.polynomial.legendre.leggauss(64);t=t*8;w=w*8*np.exp(-t*t/2)/np.sqrt(2*np.pi);w/=w.sum()
    out=[]
    for run in runs:
        j=run['lens']-1;x=run['input_mean_mm']+run['input_sigma_mm']*t;y=x-q[4*j+1]
        E=run['energy_MeV'];p=np.sqrt(E*(E+2*.51099895))/1000
        k=.299792458*model['gradient_T_per_m_per_A']*model['effective_length_m']*run['current_A']*q[4*j]/p
        R=model['radius_mm'];xx=x-run['drift_m']*k*(y+q[4*j+2]*y**3/R**2+q[4*j+3]*y**5/R**4)
        sd=np.hypot(run['drift_m']*run['input_angle_sigma_mrad'],run['resolution_mm'])
        cdf=ndtr((np.array(run['edges_mm'])[:,None]-xx)/sd)@w;prob=np.maximum(np.diff(np.r_[0,cdf,1]),1e-300)
        n=np.array(run['counts']);mu=n.sum()*prob;v=mu-n;use=n>0;v[use]+=n[use]*np.log(n[use]/mu[use])
        out.extend(np.sign(mu-n)*np.sqrt(np.maximum(2*v,0)))
    return np.array(out)

def independent_profile(qhat,dhat,runs,calcase,case):
    lo,hi=np.array(calcase['parameter_bounds']).T;span=hi-lo
    qcenter=qhat.copy();base=Independent(qhat,case);z0=base.fold()[:2];foldstart=np.r_[z0,qhat[0]]
    cache={};witness={}
    def objective(g):
        def residual(v):
            q=qcenter.copy();q[1:]=lo[1:]+span[1:]*v
            inner=Independent(q,case).fold(g,foldstart);q[0]=inner[2]
            r=calibration(q,runs,calcase['model'])
            penalty=10000*max(0,lo[0]-q[0],q[0]-hi[0])
            return np.r_[r,penalty]
        f=least_squares(residual,(qhat[1:]-lo[1:])/span[1:],bounds=(np.zeros(7),np.ones(7)),
            x_scale='jac',diff_step=2e-5,ftol=2e-10,xtol=2e-10,gtol=2e-6,max_nfev=80)
        q=qcenter.copy();q[1:]=lo[1:]+span[1:]*f.x;q[0]=Independent(q,case).fold(g,foldstart)[2]
        d=float(np.sum(calibration(q,runs,calcase['model'])**2));witness[g]=dict(strength=g,parameters=q.tolist(),deviance=d)
        return d-dhat-calcase['likelihood_delta']
    gc=base.fold()[2];out={}
    for sign,side in [(-1,'minimum'),(1,'maximum')]:
        a,b=sorted([gc,gc+sign*.019]);g=brentq(objective,a,b,xtol=2e-7)
        objective(g);out[side]=witness[g]
    return out

class SeedChart:
    def __init__(self,q,case,z,g,stable):
        mp.mp.dps=70;self.q=[mp.mpf(str(v)) for v in q];self.R=mp.mpf(str(case['radius_mm']));self.L=mp.mpf(str(case['between_lenses_m']));self.H=mp.mpf(str(case['return_drift_m']));self.r=mp.mpf(str(case['strength_ratio']));self.U=self.R/self.L;self.g=mp.mpf(str(g));self.stable=stable
        def diff(x,y):
            p=mp.matrix([x,y]);return self.iterate(p,False)-p
        self.z=mp.findroot(diff,(mp.mpf(str(z[0])),mp.mpf(str(z[1]))),tol=mp.mpf('1e-60'))
        A=mp.matrix([[mp.diff(lambda t:self.iterate(self.z+mp.matrix([t,0]),False)[i],0),mp.diff(lambda t:self.iterate(self.z+mp.matrix([0,t]),False)[i],0)] for i in range(2)])
        tr=A[0,0]+A[1,1];self.lam=(tr-mp.sqrt(tr*tr-4))/2 if stable else (tr+mp.sqrt(tr*tr-4))/2
        v=mp.matrix([1,(self.lam-A[0,0])/A[0,1]]);self.v=v/mp.sqrt((v.T*v)[0])
        B=mp.matrix([mp.diff(lambda t:self.iterate(self.z+t*self.v,False)[i],0,2)/2 for i in range(2)])
        self.c2=mp.lu_solve(self.lam**2*mp.eye(2)-A,B)
        self.grow=1/self.lam if stable else self.lam
    def iterate(self,z,inverse):
        x,u=z[0]*self.R,z[1]*self.U
        ops=[(0,self.L,self.g),(1,self.H,self.r*self.g)]
        for _ in range(3):
            for j,L,k in (ops[::-1] if inverse else ops):
                if inverse:x-=L*u
                y=x-self.q[4*j+1];kick=k*self.q[4*j]*(y+self.q[4*j+2]*y**3/self.R**2+self.q[4*j+3]*y**5/self.R**4)
                u+=kick if inverse else -kick
                if not inverse:x+=L*u
        return mp.matrix([x/self.R,u/self.U])
    def eval(self,t):
        t=mp.mpf(str(t));m=max(0,int(mp.ceil(mp.log(abs(t)/mp.mpf('1e-8'))/mp.log(self.grow))))
        h=t/self.grow**m;z=self.z+h*self.v+h*h*self.c2
        for _ in range(m):z=self.iterate(z,self.stable)
        return np.array([float(z[0]),float(z[1])])

def independent_manifolds(q,case,cycles):
    M=Independent(q,case);g=case['manifolds']['strength']
    hyper=[c for c in cycles if c['trace']>2]
    z=np.array(sorted(hyper,key=lambda c:c['states'][0])[0]['states'][0])/[M.R,M.U]
    target=M.map(z*[M.R,M.U],g)/[M.R,M.U]
    U=SeedChart(q,case,z,g,False);S=SeedChart(q,case,target,g,True)
    cfg=case['manifolds'];tu=np.linspace(*cfg['unstable_parameter_interval'],1000);ts=np.linspace(*cfg['stable_parameter_interval'],1000)
    A=np.array([U.eval(t) for t in tu]);B=np.array([S.eval(t) for t in ts]);found=[]
    # Separately bracket the two-dimensional segment intersection determinant.
    for i in range(len(tu)-1):
        l=np.minimum(A[i],A[i+1]);h=np.maximum(A[i],A[i+1]);inds=np.where((np.minimum(B[:-1],B[1:])<=h).all(axis=1)&(np.maximum(B[:-1],B[1:])>=l).all(axis=1))[0]
        for j in inds:
            mat=np.stack([A[i+1]-A[i],B[j]-B[j+1]],axis=1)
            if abs(np.linalg.det(mat))<1e-18:continue
            w=np.linalg.solve(mat,B[j]-A[i])
            if w.min()<0 or w.max()>1:continue
            f=root(lambda a:U.eval(a[0])-S.eval(a[1]),[tu[i]+w[0]*(tu[i+1]-tu[i]),ts[j]+w[1]*(ts[j+1]-ts[j])],tol=1e-10)
            if np.linalg.norm(f.fun)<1e-9 and not any(np.linalg.norm(f.x-v['parameters'])<1e-5 for v in found):found.append(dict(parameters=f.x.tolist(),point=(U.eval(f.x[0])*[M.R,M.U]).tolist()))
    found=sorted(found,key=lambda v:v['parameters'][0]);areas=[]
    for a,b in zip(found[:-1],found[1:]):
        up=np.array([U.eval(t) for t in np.linspace(a['parameters'][0],b['parameters'][0],700)])
        sp=np.array([S.eval(t) for t in np.linspace(b['parameters'][1],a['parameters'][1],700)])
        curve=np.r_[up,sp]
        area=.5*np.sum(curve[:,0]*np.roll(curve[:,1],-1)-curve[:,1]*np.roll(curve[:,0],-1))*M.R*M.U
        areas.append(float(area))
    return dict(intersections=found,signed_areas=areas,method='70-digit quadratic local seeds, repeated inverse/forward maps, independent segment bracketing and polygon area')

def main(data,evidence):
    t=time.monotonic();D=Path(data);case=json.loads((D/'lattice_case.json').read_text());cc=json.loads((D/'delivery_case.json').read_text());runs=json.loads((D/'delivery_calibration.json').read_text())['runs']
    lo,hi=np.array(cc['parameter_bounds']).T
    f=least_squares(lambda q:calibration(q,runs,cc['model']),(lo+hi)/2,bounds=(lo,hi),x_scale=hi-lo,ftol=1e-12,xtol=1e-12,gtol=1e-8)
    q=f.x;d=float(np.sum(f.fun**2));M=Independent(q,case);out=dict(parameters=q.tolist(),deviance=d,orbit_sets=[])
    for query in case['orbit_queries']:
        cycles=M.shooting(query['strength'],query['period']);out['orbit_sets'].append(dict(id=query['id'],cycles=cycles));print('independent cycles',query['id'],len(cycles),flush=True)
    out['fold']=M.fold().tolist();out['forms']=M.forms();print('independent forms',out['forms'],flush=True)
    out['confidence']=independent_profile(q,d,runs,cc,case);print('independent endpoints',out['confidence'],flush=True)
    cyc=next(v['cycles'] for v,query in zip(out['orbit_sets'],case['orbit_queries']) if query['period']==3 and query['strength']==case['manifolds']['strength'])
    out['manifolds']=independent_manifolds(q,case,cyc);print('independent manifolds',out['manifolds'],flush=True)
    out['runtime_sec']=time.monotonic()-t;Path(evidence).write_text(json.dumps(out,indent=2)+'\n')
    print('INDEPENDENT DONE',out['runtime_sec'],flush=True)

if __name__=='__main__':main(sys.argv[1],sys.argv[2])
