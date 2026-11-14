"""Symplectic plasma-cell resonances and parameterized invariant manifolds."""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import root,brentq,minimize
from scipy.spatial import cKDTree
from delivery import NAMES,deviance,deviance_residual,fit

class Lattice:
    def __init__(self,q,case):
        self.q=np.asarray(q,float);self.case=case
        self.R=case['radius_mm'];self.L=case['between_lenses_m'];self.H=case['return_drift_m']
        self.U=self.R/self.L;self.ratio=case['strength_ratio']
    def step(self,z,g,inverse=False,jac=False):
        z=np.asarray(z);scalar=z.ndim==1;z=np.atleast_2d(z)
        x=z[:,0]*self.R;u=z[:,1]*self.U
        J=np.broadcast_to(np.eye(2),(len(z),2,2)).copy()
        ops=[('kick',0,g),('drift',self.L,0),('kick',1,self.ratio*g),('drift',self.H,0)]
        if inverse:ops=ops[::-1]
        sign=-1 if inverse else 1
        for kind,a,b in ops:
            A=np.broadcast_to(np.eye(2),(len(z),2,2)).copy()
            if kind=='drift':
                x=x+sign*a*u;A[:,0,1]=sign*a*self.U/self.R
            else:
                j=a;y=(x-self.q[4*j+1])/self.R
                C,Q=self.q[4*j+2:4*j+4]
                u=u-sign*b*self.q[4*j]*self.R*(y+C*y**3+Q*y**5)
                A[:,1,0]=-sign*b*self.q[4*j]*(1+3*C*y*y+5*Q*y**4)*self.R/self.U
            J=A@J
        ans=np.stack([x/self.R,u/self.U],axis=1)
        if scalar:ans,J=ans[0],J[0]
        return (ans,J) if jac else ans
    def iterate(self,z,g,n,inverse=False,jac=False):
        z=np.asarray(z);J=np.broadcast_to(np.eye(2),z.shape[:-1]+(2,2)).copy()
        for _ in range(n):z,A=self.step(z,g,inverse,True);J=A@J
        return (z,J) if jac else z
    def physical(self,z,g,n):
        z=np.asarray(z);scalar=z.ndim==1;z=np.atleast_2d(z);margin=np.full(len(z),np.inf)
        for _ in range(n):
            box=np.array(self.case.get('survey_box',[[-1,1],[-1.35,1.35]]))
            margin=np.minimum(margin,np.minimum((z-box[:,0]).min(axis=1),(box[:,1]-z).min(axis=1))*self.R)
            x=z[:,0]*self.R;u=z[:,1]*self.U;y=(x-self.q[1])/self.R
            u=u-g*self.q[0]*self.R*(y+self.q[2]*y**3+self.q[3]*y**5)
            margin=np.minimum(margin,np.minimum(self.R-abs(x),self.R-abs(x+self.L*u)))
            z=self.step(z,g)
        return float(margin[0]) if scalar else margin
    def poly(self,w,g,nperiod=1,inverse=False):
        n=w.shape[1]-1;x=w[0].copy()*self.R;u=w[1].copy()*self.U
        def mul(a,b):return np.convolve(a,b)[:n+1]
        ops=[('kick',0,g),('drift',self.L,0),('kick',1,self.ratio*g),('drift',self.H,0)]
        if inverse:ops=ops[::-1]
        sign=-1 if inverse else 1
        for _ in range(nperiod):
            for kind,a,b in ops:
                if kind=='drift':x=x+sign*a*u
                else:
                    j=a;y=x.copy()/self.R;y[0]-=self.q[4*j+1]/self.R
                    yy=mul(y,y);yyy=mul(yy,y);yyyyy=mul(yyy,yy)
                    u=u-sign*b*self.q[4*j]*self.R*(y+self.q[4*j+2]*yyy+self.q[4*j+3]*yyyyy)
        return np.stack([x/self.R,u/self.U])
    def B(self,z,g,a,b,n=1):
        def quad(v):
            w=np.zeros((2,3),dtype=np.result_type(v,float));w[:,0]=z;w[:,1]=v
            return 2*self.poly(w,g,n)[:,2]
        return (quad(a+b)-quad(a-b))/4
    def C(self,z,g,a,b,c,n=1):
        ans=np.zeros(2,dtype=np.result_type(a,b,c,float))
        for s in [-1,1]:
            for t in [-1,1]:
                for u in [-1,1]:
                    w=np.zeros((2,4),dtype=ans.dtype);w[:,0]=z;w[:,1]=s*a+t*b+u*c
                    ans+=s*t*u*self.poly(w,g,n)[:,3]
        return ans/8
    def center(self,g):
        f=root(lambda z:(self.step(z,g,jac=True)[0]-z,self.step(z,g,jac=True)[1]-np.eye(2)),[0,0],jac=True,tol=1e-11)
        if np.linalg.norm(f.fun)>1e-10:raise RuntimeError('centroid did not converge')
        return f.x
    def orbit(self,z,g,n):
        seq=[];action=0.
        for _ in range(n):
            seq.append(z.copy());x=z[0]*self.R;u=z[1]*self.U
            for j,L,k in [(0,self.L,g),(1,self.H,self.ratio*g)]:
                y=x-self.q[4*j+1];C,Q=self.q[4*j+2:4*j+4]
                potential=y*y/2+C*y**4/(4*self.R**2)+Q*y**6/(6*self.R**4)
                u=u-k*self.q[4*j]*(y+C*y**3/self.R**2+Q*y**5/self.R**4)
                action+=L*u*u/2-k*self.q[4*j]*potential;x+=L*u
            z=np.array([x/self.R,u/self.U])
        seq=np.array(seq);i=np.lexsort((seq[:,1],seq[:,0]))[0];seq=np.roll(seq,-i,axis=0)
        _,M=self.iterate(seq[0],g,n,jac=True);tr=float(np.trace(M))
        eig=np.linalg.eigvals(M);eig=sorted(eig,key=lambda e:(e.real,e.imag))
        return dict(states=(seq*np.array([self.R,self.U])).tolist(),trace=tr,
            residue=(2-tr)/4,class_= 'elliptic' if abs(tr)<2 else 'hyperbolic',
            multipliers=[[float(e.real),float(e.imag)] for e in eig],action_mm_mrad=float(action))
    def cycles(self,g,n,grid=25):
        found=[]
        for x in np.linspace(-1,1,grid):
            for u in np.linspace(-1.5,1.5,grid):
                def fun(z):
                    y,J=self.iterate(z,g,n,jac=True);return y-z,J-np.eye(2)
                with np.errstate(over='ignore',invalid='ignore'):
                    f=root(fun,[x,u],jac=True,tol=3e-11)
                z=f.x
                if not np.isfinite(z).all() or np.linalg.norm(fun(z)[0])>1e-8 or abs(z).max()>2:continue
                if self.physical(z,g,n)<-1e-10:continue
                if any(n%d==0 and np.linalg.norm(self.iterate(z,g,d)-z)<1e-6 for d in range(1,n)):continue
                o=self.orbit(z,g,n);first=np.array(o['states'][0])/[self.R,self.U]
                if any(np.linalg.norm(first-np.array(v['states'][0])/[self.R,self.U])<1e-6 for v in found):continue
                o['classification']=o.pop('class_');found.append(o)
        return sorted(found,key=lambda o:tuple(o['states'][0]))
    def complex_mode(self,z,g):
        _,A=self.step(z,g,jac=True);v,Q=np.linalg.eig(A)
        cand=[]
        J=np.array([[0,1],[-1,0]])
        for i in range(2):
            q=Q[:,i].astype(complex);s=(q.conj()@J@q).imag
            if s>0:
                q/=np.sqrt(s);q*=np.exp(-1j*np.angle(q[0]));lam=v[i]
                p=np.linalg.inv(np.stack([q,q.conj()],axis=1))[0].conj()
                return A,lam,q,p
        raise RuntimeError('not an ordinary elliptic point')
    def normal_forms(self):
        lo,hi=self.case['resonance_interval'];g=brentq(lambda g:np.trace(self.step(self.center(g),g,jac=True)[1])+1,lo,hi,xtol=1e-12)
        z=self.center(g);A,lam,q,p=self.complex_mode(z,g)
        chi=np.vdot(p,self.B(z,g,q.conj(),q.conj()))/2
        h=self.case['twist_strength'];z0=self.center(h);A0,l0,q0,p0=self.complex_mode(z0,h)
        h20=np.linalg.solve(l0*l0*np.eye(2)-A0,self.B(z0,h,q0,q0))
        h11=np.linalg.solve(np.eye(2)-A0,self.B(z0,h,q0,q0.conj()))
        c=np.vdot(p0,self.C(z0,h,q0,q0,q0.conj())/2+self.B(z0,h,q0,h11)+self.B(z0,h,q0.conj(),h20)/2)
        return dict(resonance_strength=float(g),resonance_center=(z*[self.R,self.U]).tolist(),
            resonance_coefficient=[float(chi.real),float(chi.imag)],twist_coefficient=float((c/l0).imag))
    def saddle_center(self,start=None):
        cfg=self.case['saddle_center'];box=np.array(cfg['representative_box'])
        start=np.r_[box.mean(axis=1),np.mean(cfg['strength_interval'])] if start is None else start
        def fun(v):
            z,g=v[:2],v[2];zz,M=self.iterate(z,g,3,jac=True)
            return np.r_[zz-z,np.trace(M)-2]
        f=root(fun,start,tol=1e-10)
        z,g=f.x[:2],float(f.x[2])
        if np.linalg.norm(f.fun)>3e-8 or not cfg['strength_interval'][0]<g<cfg['strength_interval'][1] or np.any(z<box[:,0]) or np.any(z>box[:,1]):raise RuntimeError('wrong saddle-center branch')
        _,M=self.iterate(z,g,3,jac=True);_,_,V=np.linalg.svd(M-np.eye(2));v=V[-1];v*=1 if v[0]>0 else -1
        w=np.array([-v[1],v[0]]);eps=2e-5
        fg=(self.iterate(z,g+eps,3)-self.iterate(z,g-eps,3))/(2*eps)
        o=self.orbit(z,g,3);o['classification']='parabolic';o.pop('class_',None)
        sigma=float(v@(M-np.eye(2))@w);mu=float(w@fg);nu=float(w@self.B(z,g,v,v,3)/2)
        return dict(strength=g,orbit=o,critical_direction=v.tolist(),shear=sigma,parameter_coefficient=mu,quadratic_coefficient=nu),np.r_[z,g]

class Chart:
    def __init__(self,lattice,z,g,stable=False,order=18,local=.002):
        self.lattice=lattice;self.z=np.asarray(z);self.g=g;self.stable=stable;self.local=local
        _,A=lattice.iterate(z,g,3,jac=True);v,V=np.linalg.eig(A)
        i=np.argmin(v.real) if stable else np.argmax(v.real)
        self.lam=float(v[i].real);vec=V[:,i].real;vec/=np.linalg.norm(vec)
        if vec[0]<0:vec=-vec
        self.c=np.zeros((2,order+1));self.c[:,0]=z;self.c[:,1]=vec
        for n in range(2,order+1):
            rhs=lattice.poly(self.c,g,3)[:,n]
            self.c[:,n]=np.linalg.solve(self.lam**n*np.eye(2)-A,rhs)
    def eval(self,t,derivative=False):
        scalar=np.ndim(t)==0;t=np.atleast_1d(t).astype(float)
        grow=1/self.lam if self.stable else self.lam
        m=np.maximum(0,np.ceil(np.log(np.maximum(abs(t),1e-300)/self.local)/np.log(grow)).astype(int))
        local=t/grow**m;z=np.polynomial.polynomial.polyval(local,self.c.T).T
        d=np.polynomial.polynomial.polyval(local,(self.c[:,1:]*np.arange(1,self.c.shape[1])).T).T/(grow**m)[:,None]
        for n in range(m.max(initial=0)):
            mask=m>n;zz,J=self.lattice.iterate(z[mask],self.g,3,self.stable,True)
            z[mask]=zz;d[mask]=np.einsum('nij,nj->ni',J,d[mask])
        if scalar:z,d=z[0],d[0]
        return (z,d) if derivative else z

def intersections(U,S,boundsU,boundsS,grid=1600):
    tu=np.linspace(*boundsU,grid);ts=np.linspace(*boundsS,grid)
    A=U.eval(tu);B=S.eval(ts);ac=(A[:-1]+A[1:])/2;bc=(B[:-1]+B[1:])/2
    al=np.linalg.norm(np.diff(A,axis=0),axis=1);bl=np.linalg.norm(np.diff(B,axis=0),axis=1)
    tree=cKDTree(bc);out=[]
    for i in range(len(ac)):
        for j in tree.query_ball_point(ac[i],(al[i]+bl.max())/2):
            matrix=np.stack([A[i+1]-A[i],B[j]-B[j+1]],axis=1)
            if abs(np.linalg.det(matrix))<1e-18:continue
            c=np.linalg.solve(matrix,B[j]-A[i])
            if c.min()<0 or c.max()>1:continue
            start=[tu[i]+c[0]*(tu[i+1]-tu[i]),ts[j]+c[1]*(ts[j+1]-ts[j])]
            def fun(v):
                a,da=U.eval(v[0],True);b,db=S.eval(v[1],True)
                return a-b,np.stack([da,-db],axis=1)
            f=root(fun,start,jac=True,tol=1e-10)
            if np.linalg.norm(f.fun)>2e-9 or not boundsU[0]<f.x[0]<boundsU[1] or not boundsS[0]<f.x[1]<boundsS[1]:continue
            if any(np.linalg.norm(f.x-v['parameters'])<1e-6 for v in out):continue
            z,du=U.eval(f.x[0],True);_,ds=S.eval(f.x[1],True)
            sin=float(np.linalg.det(np.stack([du,ds]))/(np.linalg.norm(du)*np.linalg.norm(ds)))
            out.append(dict(parameters=f.x.tolist(),point=(z*[U.lattice.R,U.lattice.U]).tolist(),crossing_sine=sin))
    return sorted(out,key=lambda v:v['parameters'][0])

def arc_integral(chart,a,b,order=96):
    x,w=np.polynomial.legendre.leggauss(order);t=(a+b)/2+(b-a)*x/2
    z,d=chart.eval(t,True)
    return float((b-a)/4*np.dot(w,z[:,0]*d[:,1]-z[:,1]*d[:,0]))

def manifolds(lattice,cycles=None,grid=1600,order=18,local=.002,quadrature=96):
    cfg=lattice.case['manifolds'];g=cfg['strength']
    cycles=lattice.cycles(g,3) if cycles is None else cycles
    first=sorted([o for o in cycles if o['classification']=='hyperbolic'],key=lambda o:tuple(o['states'][0]))[0]
    z=np.array(first['states'][0])/[lattice.R,lattice.U];target=lattice.step(z,g)
    U=Chart(lattice,z,g,False,order,local);S=Chart(lattice,target,g,True,order,local)
    rows=intersections(U,S,cfg['unstable_parameter_interval'],cfg['stable_parameter_interval'],grid)
    lobes=[]
    for i in range(len(rows)-1):
        a,b=rows[i]['parameters'],rows[i+1]['parameters']
        if any(min(a[1],b[1])<v['parameters'][1]<max(a[1],b[1]) for j,v in enumerate(rows) if j not in [i,i+1]):continue
        area=(arc_integral(U,a[0],b[0],quadrature)+arc_integral(S,b[1],a[1],quadrature))*lattice.R*lattice.U
        lobes.append(dict(intersection_indices=[i,i+1],signed_area_mm_mrad=area))
    t=np.array(cfg['sample_parameters']);up=U.eval(t[:,0])*[lattice.R,lattice.U];sp=S.eval(t[:,1])*[lattice.R,lattice.U]
    return dict(intersections=rows,lobes=lobes,unstable_samples=up.tolist(),stable_samples=sp.tolist())

def calculate(q,case,grid=25,manifold_grid=1600):
    L=Lattice(q,case);sets=[]
    for query in case['orbit_queries']:
        cycles=L.cycles(query['strength'],query['period'],grid)
        sets.append(dict(id=query['id'],cycles=cycles))
    forms=L.normal_forms();fold,seed=L.saddle_center()
    cfg=case['saddle_center'];near=[]
    z,g=seed[:2],seed[2];v=np.array(fold['critical_direction'])
    for delta in cfg['continuation_offsets']:
        eps=np.sqrt(abs(fold['parameter_coefficient']*delta/fold['quadratic_coefficient']))
        pair=[]
        for sign in [-1,1]:
            h=g+delta
            f=root(lambda y:(L.iterate(y,h,3,jac=True)[0]-y,L.iterate(y,h,3,jac=True)[1]-np.eye(2)),z+sign*eps*v,jac=True,tol=1e-11)
            if np.linalg.norm(f.fun)>1e-8:raise RuntimeError('secondary orbit did not converge')
            o=L.orbit(f.x,h,3);o['classification']=o.pop('class_');pair.append(o)
        near.append(dict(offset=delta,cycles=sorted(pair,key=lambda o:tuple(o['states'][0]))))
    matching=next((v['cycles'] for qry,v in zip(case['orbit_queries'],sets) if qry['strength']==case['manifolds']['strength'] and qry['period']==3),None)
    return dict(orbit_sets=sets,normal_forms=forms,saddle_center=fold,secondary_branches=near,
                manifolds=manifolds(L,matching,manifold_grid))

def packet_capture(q,case,source,manifold_result,order_x=96,order_y=32):
    """Change of variables through the symplectic map, with exact wall history."""
    from scipy.special import logsumexp
    L=Lattice(q,case);cfg=case['packet'];g=case['manifolds']['strength']
    cycles=L.cycles(g,3,21)
    o=sorted([v for v in cycles if v['classification']=='hyperbolic'],key=lambda v:tuple(v['states'][0]))[0]
    z=np.array(o['states'][0])/[L.R,L.U];target=L.step(z,g)
    U=Chart(L,z,g);S=Chart(L,target,g,True)
    rows=manifold_result['intersections'];a,b=[rows[i]['parameters'] for i in cfg['lobe_intersection_indices']]
    tx,wx=np.polynomial.legendre.leggauss(order_x);ty,wy=np.polynomial.legendre.leggauss(order_y)
    t=(a[0]+b[0])/2+(b[0]-a[0])*tx/2
    up,du=U.eval(t,True);stable=[]
    for xx in up[:,0]:
        h=brentq(lambda v:S.eval(v)[0]-xx,min(a[1],b[1]),max(a[1],b[1]),xtol=1e-12)
        stable.append(S.eval(h)[1])
    height=np.array(stable)-up[:,1]
    vertices=np.stack([np.repeat(up[:,0],order_y),np.repeat(up[:,1],order_y)+np.repeat(height,order_y)*np.tile((ty+1)/2,order_x)],axis=1)
    weights=(np.abs(du[:,0]*height)*wx*(b[0]-a[0])/2)[:,None]*wy[None,:]/2*L.R*L.U
    weights=weights.ravel()
    sx=source['source_size_um']*.001;su=source['divergence_mrad'];d=source['source_distance_mm']*.001+cfg['upstream_drift_m']
    C0=np.diag([sx*sx,su*su,cfg['position_jitter_mm']**2,cfg['angle_jitter_mrad']**2])
    B=np.array([[1,d,1,0],[0,1,0,1]],float);cov=B@C0@B.T;inv=np.linalg.inv(cov)
    gain=C0@B.T@inv;conditional_cov=C0-gain@B@C0
    center=np.array(cfg['centroid_mm_mrad']);normal=2*np.pi*np.sqrt(np.linalg.det(cov))
    outputs=[]
    for query in cfg['queries']:
        n=query['cells'];checkpoint=query['checkpoint_cell'];z=vertices.copy();live=np.ones(len(z),bool);mid=None
        for i in range(n):
            z=L.step(z,g,True)
            x=z[:,0]*L.R;u=z[:,1]*L.U;y=(x-L.q[1])/L.R
            after=u-g*L.q[0]*L.R*(y+L.q[2]*y**3+L.q[3]*y**5)
            live&=(abs(x)<=L.R)&(abs(x+L.L*after)<=L.R)
            if n-i-1==checkpoint:mid=z*np.array([L.R,L.U])
        initial=z*np.array([L.R,L.U]);v=initial-center
        logs=-.5*np.einsum('ni,ij,nj->n',v,inv,v)-np.log(normal)+np.log(weights)
        logs[~live]=-np.inf;logp=float(logsumexp(logs));w=np.exp(logs-logp)
        latent=v@gain.T;mean=w@latent;cv=conditional_cov+(latent-mean).T@(w[:,None]*(latent-mean))
        mm=w@mid;mc=(mid-mm).T@(w[:,None]*(mid-mm))
        outputs.append(dict(id=query['id'],probability=float(np.exp(logp)),latent_mean=mean.tolist(),latent_covariance=cv.tolist(),checkpoint_mean=mm.tolist(),checkpoint_covariance=mc.tolist()))
    return outputs

def confidence(qhat,dhat,runs,calcase,case):
    from scipy.optimize import least_squares
    lo,hi=np.array(calcase['parameter_bounds']).T
    jac=least_squares(deviance_residual,qhat,args=(runs,calcase['model']),max_nfev=2).jac
    _,s,V=np.linalg.svd(jac,full_matrices=False);W=V.T/s*np.sqrt(calcase['likelihood_delta'])
    limit=dhat+calcase['likelihood_delta']-.015
    initial=Lattice(qhat,case).saddle_center()[1];cache={}
    def evaluate(y):
        q=qhat+W@y;key=q.tobytes()
        if key not in cache:
            fold,_=Lattice(q,case).saddle_center(initial)
            cache[key]=(fold['strength'],deviance(q,runs,calcase['model']),q)
        return cache[key]
    def box(y):
        q=qhat+W@y;return np.r_[(q-lo)/(hi-lo),(hi-q)/(hi-lo)]
    results={}
    for sign,label in [(1,'minimum'),(-1,'maximum')]:
        candidates=[]
        for i in range(2):
            f=minimize(lambda y:sign*evaluate(y)[0],np.zeros(8) if i==0 else .15*np.sin(np.arange(8)+2*i),method='SLSQP',
                constraints=[dict(type='ineq',fun=lambda y:(limit-evaluate(y)[1])/20),dict(type='ineq',fun=box)],
                options=dict(ftol=2e-10,eps=2e-5,maxiter=100))
            val,d,q=evaluate(f.x)
            if d<=limit+.02 and box(f.x).min()>-1e-8:candidates.append((val,d,q))
        if not candidates:raise RuntimeError('no feasible saddle-center endpoint')
        val,d,q=min(candidates,key=lambda c:sign*c[0])
        results[label]=dict(strength=float(val),parameters=dict(zip(NAMES,map(float,q))),deviance=float(d))
    return results

def main(data,output):
    data=Path(data);t=time.monotonic();out=json.loads(Path(output).read_text())
    case=json.loads((data/'lattice_case.json').read_text());calcase=json.loads((data/'delivery_case.json').read_text())
    runs=json.loads((data/'delivery_calibration.json').read_text())['runs']
    q,d,_=fit(runs,calcase)
    out['delivery_parameters']=dict(zip(NAMES,map(float,q)));out['calibration_deviance']=d
    out['lattice']=calculate(q,case)
    if 'packet' in case:out['lobe_capture']=packet_capture(q,case,out,out['lattice']['manifolds'])
    out['lattice_confidence']=confidence(q,d,runs,calcase,case)
    p=np.sqrt(out['energy_MeV']*(out['energy_MeV']+2*.51099895))/1000
    out['resonance_currents_A']=(out['lattice']['normal_forms']['resonance_strength']*p/(.299792458*calcase['model']['gradient_T_per_m_per_A']*calcase['model']['effective_length_m'])*np.array([1,case['strength_ratio']])).tolist()
    Path(output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('LATTICE DONE',time.monotonic()-t,flush=True)

if __name__=='__main__':
    import sys
    main(sys.argv[1],sys.argv[2])
