"""Joint generalized least squares with exact linear-nuisance profiling."""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares,lsq_linear
from heterogeneity import KEYS,SCALE,NP,predict,design,whiten,thermal_knots,temperature

def calibrations(x,data,sens=True):
    residuals=[];jac=[]
    c=data['thermal_calibration'];knots=thermal_knots(c['program'],x[12],x[13])
    for t,y in zip(c['times'],c['sample_temperature_K']):
        T,d=temperature(t,knots,x[12],x[13]);row=np.zeros(NP);row[12]=d/c['sigma_K'];row[13]=1/c['sigma_K']
        residuals.append((T-y)/c['sigma_K']);jac.append(row)
    c=data['CO_pulse'];t=np.asarray(c['times']);f=c['initial_current_pA']*np.exp(-t/x[10]);r=(f-c['current_28_pA'])/c['sigma_pA']
    J=np.zeros((len(t),NP));J[:,10]=f*t/x[10]**2/c['sigma_pA'];residuals.extend(r);jac.extend(J)
    c=data['CO2_pulse'];t=np.asarray(c['times']);f=c['initial_current_44_pA']*np.exp(-t/x[11]);model=f[:,None]*np.array([x[14],1.])
    residual=(model-c['current_28_44_pA'])/np.asarray(c['sigma_pA']);J=np.zeros((len(t),2,NP));J[:,:,11]=model*t[:,None]/x[11]**2
    J[:,0,14]=f;J/=np.asarray(c['sigma_pA'])[None,:,None]
    residuals.extend(residual.ravel());jac.extend(J.reshape(-1,NP))
    return np.asarray(residuals),np.asarray(jac)

def project_run(x,run,cov,taucorr,sens=True,rtol=3e-9):
    pred=predict(x,run,sens=sens,rtol=rtol);A,dA=design(x,run,pred);times=np.asarray(run['times'])
    B=whiten(A,times,cov,taucorr).reshape(-1,5);Y=whiten(np.asarray(run['ion_current_pA']),times,cov,taucorr).ravel()
    lower=np.array([.6,-.25,-.25,-.25,-.25]);upper=np.array([1.7,.25,.25,.25,.25])
    coeff=np.linalg.lstsq(B,Y,rcond=None)[0]
    if np.any(coeff<lower) or np.any(coeff>upper):
        sol=lsq_linear(B,Y,bounds=(lower,upper),tol=1e-12);coeff=sol.x;free=sol.active_mask==0
    else:free=np.ones(5,bool)
    residual=B@coeff-Y
    J=np.zeros((len(Y),NP))
    if sens:
        dB=whiten(dA,times,cov,taucorr).reshape(-1,5,NP)
        direct=np.einsum('nkp,k->np',dB,coeff);BF=B[:,free]
        normal=BF.T@BF
        rhs=BF.T@direct+np.einsum('nkp,n->kp',dB[:,free],residual)
        dc=-np.linalg.solve(normal,rhs);J=direct+BF@dc
    return residual,J,coeff,pred

class Objective:
    def __init__(self,case,data,rtol=3e-9,verbose=True):
        self.case=case;self.data=data;self.rtol=rtol;self.cached=None;self.calls=0;self.verbose=verbose
    def evaluate(self,z):
        if self.cached is not None and np.array_equal(z,self.cached[0]):return self.cached[1:]
        started=time.perf_counter();x=np.asarray(z)*SCALE;res=[];jac=[];nuisance={}
        r,J=calibrations(x,self.data);res.append(r);jac.append(J)
        for run in self.data['runs']:
            r,J,c,_=project_run(x,run,self.data['noise']['channel_covariance_pA2'],self.data['noise']['correlation_time_s'],True,self.rtol)
            res.append(r);jac.append(J);nuisance[run['id']]=c.tolist()
        residual=np.concatenate(res);J=np.vstack(jac)*SCALE[None,:];self.calls+=1
        self.cached=(np.asarray(z).copy(),residual,J,nuisance)
        if self.verbose:print('GLS',self.calls,'D',float(residual@residual),'seconds',round(time.perf_counter()-started,2),'x',x.tolist(),flush=True)
        return residual,J,nuisance
    def fun(self,z):return self.evaluate(z)[0]
    def jac(self,z):return self.evaluate(z)[1]

def fit(case,data,start=None,max_nfev=90,rtol=3e-9):
    bounds=np.asarray([case['parameter_bounds'][key] for key in KEYS]);low,high=bounds.T
    x0=(low+high)/2 if start is None else np.asarray(start)
    obj=Objective(case,data,rtol=rtol)
    result=least_squares(obj.fun,x0/SCALE,jac=obj.jac,bounds=(low/SCALE,high/SCALE),
                         ftol=2e-8,xtol=2e-8,gtol=1e-6,max_nfev=max_nfev,x_scale='jac')
    residual,J,nuisance=obj.evaluate(result.x)
    return dict(parameters=dict(zip(KEYS,(result.x*SCALE).tolist())),deviance=float(residual@residual),
                run_nuisance=nuisance,optimality=float(result.optimality),evaluations=obj.calls),obj,result
