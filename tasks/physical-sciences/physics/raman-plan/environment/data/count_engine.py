"""Matrix-valued photon probabilities from the specified Markov count law."""
import numpy as np
from scipy.linalg import expm
from scipy.special import pdtr,pdtrc
from scipy.stats import poisson

def distributions(q,rates,duration,tolerance=2e-15):
    count_bound=int(poisson.isf(tolerance,float(np.max(rates))*duration))+1
    N=1<<max(4,int(count_bound).bit_length())
    z=np.exp(2j*np.pi*np.arange(N)/N)[None,:]
    a=duration*(q[:,0,0,None]+(z-1)*rates[:,0,None]);e=duration*(q[:,1,1,None]+(z-1)*rates[:,1,None]);g=duration*(z-1)*rates[:,2,None]
    b=duration*q[:,0,1,None];d=duration*q[:,1,0,None];c=duration*q[:,0,2,None];f=duration*q[:,1,2,None]
    h=(a+e)/2;v=(a-e)/2;delta=np.sqrt(v*v+b*d)
    plus=np.exp(h+delta);minus=np.exp(h-delta);C=(plus+minus)/2
    D=np.divide(plus-minus,2*delta,out=np.zeros_like(delta),where=abs(delta)>1e-8)
    small=abs(delta)<=1e-8;D[small]=np.exp(h[small])*(1+delta[small]**2/6)
    E00=C+D*v;E01=D*b;E10=D*d;E11=C-D*v;Eg=np.exp(g)
    det=(a-g)*(e-g)-b*d
    va=((e-g)*c-b*f)/det;vb=((a-g)*f-d*c)/det
    E02=(E00-Eg)*va+E01*vb;E12=E10*va+(E11-Eg)*vb
    tilted=np.zeros((len(q),N,3,3),complex)
    tilted[:,:,0,0]=E00;tilted[:,:,0,1]=E01;tilted[:,:,0,2]=E02;tilted[:,:,1,0]=E10;tilted[:,:,1,1]=E11;tilted[:,:,1,2]=E12;tilted[:,:,2,2]=Eg
    bad=np.abs(det)<1e-11*np.maximum(1,abs(a-g)*abs(e-g)+abs(b*d))
    for s,k in zip(*np.where(bad)):
        matrix=q[s].astype(complex)+(z[0,k]-1)*np.diag(rates[s]);tilted[s,k]=expm(matrix*duration)
    pmf=(np.fft.fft(tilted,axis=1)/N).real
    assert pmf.min()>-1e-10,float(pmf.min())
    pmf=np.maximum(pmf,0.)
    material=np.asarray([expm(x*duration) for x in q]);total=pmf.sum(axis=1)
    pmf*=np.divide(material,total,out=np.zeros_like(material),where=total>0)[:,None]
    lo=pmf.cumsum(axis=1);hi=pmf[:,::-1].cumsum(axis=1)[:,::-1]
    return lo,hi

def comparisons(q,rates,duration,thresholds):
    q=np.asarray(q,float);rates=np.asarray(rates,float);ts=np.asarray(thresholds,int)
    if not np.any(q):
        lo=np.zeros((len(ts),len(q),3,3));hi=np.zeros_like(lo)
        for f in range(3):
            lo[:,:,f,f]=pdtr(ts[:,None],duration*rates[None,:,f]);hi[:,:,f,f]=pdtrc(ts[:,None],duration*rates[None,:,f])
        return lo,hi
    return tails(distributions(q,rates,duration),rates,duration,ts)

def tails(distributions,rates,duration,thresholds):
    low,high=distributions;ts=np.asarray(thresholds,int);N=low.shape[1]
    lo=np.empty((len(ts),len(rates),3,3));hi=np.empty_like(lo)
    for j,t in enumerate(ts):
        lo[j]=low[:,min(t,N-1)];hi[j]=high[:,t+1] if t<N-1 else 0.
    lo[:,:,2,:]=0.;hi[:,:,2,:]=0.
    lo[:,:,2,2]=pdtr(ts[:,None],duration*rates[None,:,2]);hi[:,:,2,2]=pdtrc(ts[:,None],duration*rates[None,:,2])
    return lo,hi
