"""Nuisance-refitted profile limits for a forecast yield difference."""
import numpy as np
from scipy.optimize import minimize,LinearConstraint
from heterogeneity import KEYS,SCALE
from forecast import forecast

def profile_limits(case,data,fit_answer,obj,start_u=None,maxiter=65,ftol=2e-4,diagnostics=None):
    z0=np.array([fit_answer['parameters'][k] for k in KEYS])/SCALE
    r,J,_=obj.evaluate(z0);minimum=float(r@r);delta=case['profile_delta']
    _,singular,Vt=np.linalg.svd(J,full_matrices=False)
    transform=Vt.T/np.maximum(singular,.1)
    bounds=np.asarray([case['parameter_bounds'][k] for k in KEYS])/SCALE[:,None]
    linear=LinearConstraint(transform,bounds[:,0]-z0,bounds[:,1]-z0)
    cache={}
    def evaluate(u):
        key=np.asarray(u).tobytes()
        if key not in cache:
            z=z0+transform@u;res,jac,nuisance=obj.evaluate(z)
            _,value,gradient=forecast(z*SCALE,case,sens=True)
            cache[key]=dict(parameters=dict(zip(KEYS,(z*SCALE).tolist())),deviance=float(res@res),run_nuisance=nuisance,
                            yield_difference=value,gradient=gradient@transform,deviance_gradient=(2*res@jac)@transform)
        return cache[key]
    answers={}
    for name,sign in (('lower',1.),('upper',-1.)):
        def fun(u):return sign*1000*evaluate(u)['yield_difference']
        def jac(u):return sign*1000*evaluate(u)['gradient']
        def constraint(u):return (minimum+delta-evaluate(u)['deviance'])/delta
        def cjac(u):return -evaluate(u)['deviance_gradient']/delta
        iterations=[0]
        def callback(u):
            e=evaluate(u);iterations[0]+=1
            print('Profile',name,iterations[0],'contrast',e['yield_difference'],'D',e['deviance'],flush=True)
        u0=np.zeros(len(z0)) if start_u is None else np.asarray(start_u[name])
        result=minimize(fun,u0,jac=jac,method='SLSQP',constraints=[linear,dict(type='ineq',fun=constraint,jac=cjac)],
                        callback=callback,options=dict(ftol=ftol,maxiter=maxiter,disp=False))
        e=evaluate(result.x)
        if diagnostics is not None:
            diagnostics[name]=dict(success=bool(result.success),message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),deviance=e['deviance'],yield_difference=e['yield_difference'])
        if e['deviance']>minimum+delta+2e-3:raise RuntimeError('Infeasible profile endpoint: '+str(result.message))
        answers[name]={k:v for k,v in e.items() if k not in ('gradient','deviance_gradient')}
        print('Profile complete',name,result.success,str(result.message),flush=True)
    return answers
