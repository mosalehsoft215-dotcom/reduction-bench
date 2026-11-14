"""Independent lobe integration in physical x, with high-precision boundary samples."""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from independent_lattice import Independent,SeedChart,NAMES

def captures(data,source_file,independent_file,evidence_file):
    start=time.monotonic();D=Path(data);case=json.loads((D/'lattice_case.json').read_text())
    src=json.loads(Path(source_file).read_text());ind=json.loads(Path(independent_file).read_text());q=np.array(ind['parameters'])
    M=Independent(q,case);g=case['manifolds']['strength'];cfg=case['packet']
    cycles=next(v['cycles'] for v in ind['orbit_sets'] if v['id']=='three-strong')
    first=sorted([c for c in cycles if c['trace']>2],key=lambda c:tuple(c['states'][0]))[0]
    z=np.array(first['states'][0])/[M.R,M.U];target=M.map(z*[M.R,M.U],g)/[M.R,M.U]
    U=SeedChart(q,case,z,g,False);S=SeedChart(q,case,target,g,True)
    a,b=[ind['manifolds']['intersections'][i]['parameters'] for i in cfg['lobe_intersection_indices']]
    sx=src['source_size_um']*.001;su=src['divergence_mrad'];d=src['source_distance_mm']*.001+cfg['upstream_drift_m']
    C0=np.diag([sx*sx,su*su,cfg['position_jitter_mm']**2,cfg['angle_jitter_mrad']**2]);H=np.array([[1,d,1,0],[0,1,0,1]])
    Sigma=H@C0@H.T;inverse=np.linalg.inv(Sigma);gain=C0@H.T@inverse;extra=C0-gain@H@C0;offset=np.array(cfg['centroid_mm_mrad'])
    report=[]
    for ns,nx,ny in [(192,96,24),(384,160,48)]:
        up=np.array([U.eval(t) for t in np.linspace(a[0],b[0],ns)])*[M.R,M.U]
        sp=np.array([S.eval(t) for t in np.linspace(a[1],b[1],ns)])*[M.R,M.U]
        sort=np.argsort(sp[:,0]);sp=sp[sort]
        lower=PchipInterpolator(up[:,0],up[:,1]);upper=PchipInterpolator(sp[:,0],sp[:,1])
        tx,wx=np.polynomial.legendre.leggauss(nx);ty,wy=np.polynomial.legendre.leggauss(ny)
        xa,xb=up[0,0],up[-1,0];x=(xa+xb)/2+(xb-xa)*tx/2;lo=lower(x);hi=upper(x)
        xs=np.repeat(x,ny);us=np.repeat(lo,ny)+np.repeat(hi-lo,ny)*np.tile((ty+1)/2,nx)
        weights=((hi-lo)*wx*(xb-xa)/2)[:,None]*wy[None,:]/2;weights=weights.ravel();out=[]
        for query in cfg['queries']:
            x=xs.copy();u=us.copy();live=np.ones(len(x),bool);checkpoint=None
            for i in range(query['cells']):
                x-=M.H*u;y=x-q[5];u+=M.r*g*q[4]*(y+q[6]*y**3/M.R**2+q[7]*y**5/M.R**4)
                x-=M.L*u;y=x-q[1];u+=g*q[0]*(y+q[2]*y**3/M.R**2+q[3]*y**5/M.R**4)
                after=u-g*q[0]*(y+q[2]*y**3/M.R**2+q[3]*y**5/M.R**4)
                live&=(abs(x)<=M.R)&(abs(x+M.L*after)<=M.R)
                if query['cells']-i-1==query['checkpoint_cell']:checkpoint=np.stack([x,u],axis=1)
            delta=np.stack([x,u],axis=1)-offset
            weight=weights*np.exp(-.5*np.einsum('ni,ij,nj->n',delta,inverse,delta))/(2*np.pi*np.sqrt(np.linalg.det(Sigma)))
            weight[~live]=0.;prob=weight.sum();weight/=prob;latent=delta@gain.T;mean=weight@latent;cm=weight@checkpoint
            cv=extra+(latent-mean).T@(weight[:,None]*(latent-mean));cc=(checkpoint-cm).T@(weight[:,None]*(checkpoint-cm))
            out.append(dict(id=query['id'],probability=float(prob),latent_mean=mean.tolist(),latent_covariance=cv.tolist(),checkpoint_mean=cm.tolist(),checkpoint_covariance=cc.tolist(),surviving_quadrature_fraction=float(live.mean())))
        record=dict(boundary_samples=ns,x_order=nx,y_order=ny,values=out);report.append(record);print(json.dumps(record),flush=True)
    Path(evidence_file).write_text(json.dumps(dict(method='Independent physical-x PCHIP boundary quadrature, direct physical inverse maps and Gaussian latent conditioning',results=report,runtime_sec=time.monotonic()-start),indent=2)+'\n')
    print('INDEPENDENT CAPTURE DONE',time.monotonic()-start,flush=True)

if __name__=='__main__':captures(*sys.argv[1:])
