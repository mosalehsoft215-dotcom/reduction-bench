"""Finite-cell reaction rates, thermal contact and species-specific QMS response."""
from functools import lru_cache
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from dynamics import construct,subset,initial,rates,rate_jac

KEYS=('E_des','log10_k_des_450','epsilon_CC','epsilon_CO','E_ox',
      'log10_k_ox_450','E_hop','log10_k_hop_450','epsilon_OO','E_hop_O',
      'tau_CO','tau_CO2','tau_T','temperature_offset','cracking_ratio')
SCALE=np.array([20.,1.,2.,2.,20.,1.,20.,1.,2.,30.,2.,3.,3.,3.,.2])
NP=len(KEYS)

@lru_cache(maxsize=1)
def whole_space():return construct()

@lru_cache(maxsize=12)
def cell(c,o):return subset(whole_space(),c,o)

def thermal_knots(program,tau,offset):
    program=np.asarray(program,float);values=[program[0,1]+offset];derivs=[0.]
    for (a,T0),(b,T1) in zip(program,program[1:]):
        dt=b-a;slope=(T1-T0)/dt;e=np.exp(-dt/tau)
        tail=values[-1]-T0-offset+slope*tau
        deriv=-slope+(derivs[-1]+slope)*e+tail*e*dt/tau**2
        value=T1+offset-slope*tau+tail*e
        values.append(value);derivs.append(deriv)
    return program,np.asarray(values),np.asarray(derivs)

def temperature(t,knots,tau,offset):
    program,values,derivs=knots
    i=min(max(int(np.searchsorted(program[:,0],t,side='right'))-1,0),len(program)-2)
    t0,T0=program[i];t1,T1=program[i+1];dt=t-t0
    slope=(T1-T0)/(t1-t0);e=np.exp(-dt/tau);tail=values[i]-T0-offset+slope*tau
    value=T0+offset+slope*(dt-tau)+tail*e
    derivative=-slope+(derivs[i]+slope)*e+tail*e*dt/tau**2
    return value,derivative

def rate_temperature_derivative(x,T,types,v):
    tag=types[:,0];E=np.where(tag==0,x[0],np.where(tag==1,x[4],np.where(tag==2,x[6],x[9])))
    lateral=x[2]*types[:,1]+x[3]*types[:,2]+x[8]*types[:,3]
    energy=np.where(tag==0,E+lateral,np.where(tag>=2,E-lateral/2,E))
    return v*energy/(.00831446261815324*T*T)

@njit
def probability_rhs(p,v,dv,src,dst,labels,tags):
    out=np.zeros_like(p);flux=np.zeros((p.shape[0],2))
    for e in range(len(src)):
        i=src[e];j=dst[e];k=labels[e];tag=tags[k]
        q=v[k]*p[0,i];out[0,i]-=q;out[0,j]+=q
        if tag<2:flux[0,tag]+=q/12.
        for a in range(p.shape[0]-1):
            q=v[k]*p[a+1,i]+dv[k,a]*p[0,i]
            out[a+1,i]-=q;out[a+1,j]+=q
            if tag<2:flux[a+1,tag]+=q/12.
    return out,flux

def predict(x,run,sens=False,rtol=3e-9,atol=2e-13,full=False):
    x=np.asarray(x,float);s=cell(*run['counts']);n=len(s['reps']);m=NP+1 if sens else 1
    y=np.zeros((m,n+2));y[0,:n]=initial(s,run.get('initial',{'counts':run['counts']}))
    prep=run.get('initial',{}).get('preparation_temperature_K') if isinstance(run.get('initial'),dict) else None
    if prep is not None:
        lateral=s['obs'][:,2:]@x[[2,3,8]]/(.00831446261815324*prep)
        support=y[0,:n]>0;weight=np.zeros(n);weight[support]=y[0,:n][support]*np.exp(lateral[support]-lateral[support].max());weight/=weight.sum();y[0,:n]=weight
        if sens:
            for k,j in enumerate((2,3,8)):y[j+1,:n]=weight*(s['obs'][:,2+k]-weight@s['obs'][:,2+k])/(.00831446261815324*prep)
    knots=thermal_knots(run['program'],x[12],x[13]);tags=s['types'][:,0]
    def rhs(t,flat):
        state=flat.reshape(m,n+2);T,dT=temperature(t,knots,x[12],x[13]);v=rates(x,T,s['types'])
        dv=np.zeros((len(v),NP))
        if sens:
            dv[:,:10]=rate_jac(x,T,s['types']);dtemp=rate_temperature_derivative(x,T,s['types'],v)
            dv[:,12]=dtemp*dT;dv[:,13]=dtemp
        dp,flux=probability_rhs(state[:,:n],v,dv,s['src'],s['dst'],s['labels'],tags)
        out=np.empty_like(state);out[:,:n]=dp;out[:,n:]=(flux-state[:,n:])/x[10:12]
        if sens:
            out[11,n]-=(flux[0,0]-state[0,n])/x[10]**2
            out[12,n+1]-=(flux[0,1]-state[0,n+1])/x[11]**2
        return out.ravel()
    times=np.asarray(run['times'],float);bytime={};flat=y.ravel();edges=sorted(set([0.,float(times[-1])]+[float(t) for t,T in run['program'] if 0<t<times[-1]]))
    if times[0]==0:bytime[0.]=flat.copy()
    for a,b in zip(edges,edges[1:]):
        requested=sorted(set(times[(times>a)&(times<=b)].tolist()+[b]))
        sol=solve_ivp(rhs,(a,b),flat,method='DOP853',t_eval=requested,rtol=rtol,atol=atol)
        if not sol.success:raise RuntimeError(sol.message)
        for t,col in zip(sol.t,sol.y.T):bytime[float(t)]=col
        flat=sol.y[:,-1]
    values=np.stack([bytime[float(t)].reshape(m,n+2) for t in times],axis=-1)
    q=values[0,n:,:].T;dq=values[1:,n:,:].transpose(2,1,0)
    c0,o0=run['counts'];nc=s['obs'][:,0]/12.;no=s['obs'][:,1]/12.
    counts=np.stack([nc@values[0,:n],no@values[0,:n]],axis=1)
    yields=np.stack([c0/12.-counts[:,0]-(o0/12.-counts[:,1]),o0/12.-counts[:,1]],axis=1)
    dyields=np.zeros((len(times),2,NP))
    if sens:
        dc=np.einsum('n,ant->ta',nc,values[1:,:n]);do=np.einsum('n,ant->ta',no,values[1:,:n])
        dyields[:,0]=-dc+do;dyields[:,1]=-do
    result=dict(q=q,dq=dq,yields=yields,dyields=dyields,coverages=counts,
                sample_temperature=np.asarray([temperature(t,knots,x[12],x[13])[0] for t in times]))
    if full:
        pairs=(s['obs'][:,2:]/24.).T@values[0,:n];flux=[]
        for j,t in enumerate(times):
            v=rates(x,result['sample_temperature'][j],s['types'])
            _,f=probability_rhs(values[:,:n,j],v,np.zeros((len(v),NP)),s['src'],s['dst'],s['labels'],tags)
            flux.append(f[0])
        result.update(pair_fractions=pairs.T,flux=np.asarray(flux),mass_error=float(np.max(abs(values[0,:n].sum(axis=0)-1))))
    return result

def whiten(values,times,covariance,correlation_time):
    values=np.asarray(values);shape=values.shape
    temp=np.linalg.solve(np.linalg.cholesky(covariance),values.reshape(len(times),3,-1).transpose(1,0,2).reshape(3,-1))
    temp=temp.reshape(3,len(times),-1).transpose(1,0,2)
    out=temp.copy();rho=np.exp(-np.diff(times)/correlation_time)
    out[1:]=(temp[1:]-rho[:,None,None]*temp[:-1])/np.sqrt(1-rho[:,None,None]**2)
    return out.reshape(shape)

def design(x,run,pred):
    t=np.asarray(run['times']);v=t/t[-1]-.5;q=pred['q'];dq=pred['dq'];N=len(t)
    A=np.zeros((N,3,5));A[:,0,0]=1200*q[:,0]+x[14]*850*q[:,1];A[:,1,0]=850*q[:,1];A[:,2,0]=run['argon_reference_pA']
    A[:,0,1]=1;A[:,0,2]=v;A[:,1,3]=1;A[:,1,4]=v
    dA=np.zeros((N,3,5,NP))
    if dq.shape[-1]:
        dA[:,0,0]=1200*dq[:,0]+x[14]*850*dq[:,1];dA[:,0,0,14]+=850*q[:,1];dA[:,1,0]=850*dq[:,1]
    return A,dA
