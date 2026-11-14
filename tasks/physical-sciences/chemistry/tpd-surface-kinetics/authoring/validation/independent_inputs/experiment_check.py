"""Independent sparse generator and coupled temperature/filter integration."""
from functools import lru_cache
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import cholesky,solve_triangular
from lattice_check import build

PARAMETERS=('E_des','log10_k_des_450','epsilon_CC','epsilon_CO','E_ox','log10_k_ox_450','E_hop','log10_k_hop_450','epsilon_OO','E_hop_O','tau_CO','tau_CO2','tau_T','temperature_offset','cracking_ratio','fraction_B','delta_E_des_B','delta_E_ox_B','delta_E_hop_B')

@lru_cache(maxsize=10)
def space(c,o):return build(c,o)

def integrate(rhs,y,program,times,rtol=5e-10,atol=1e-12):
    program=np.asarray(program);times=np.asarray(times,float)
    knots=sorted(set([0.,float(times[-1])]+[float(t) for t,T in program if 0<t<times[-1]]));out={}
    if times[0]==0:out[0.]=y.copy()
    for left,right in zip(knots,knots[1:]):
        requested=sorted(set([float(t) for t in times if left<t<=right]+[right]))
        sol=solve_ivp(rhs,(left,right),y,method='RK45',rtol=rtol,atol=atol,t_eval=requested)
        if not sol.success:raise ArithmeticError(sol.message)
        for t,v in zip(sol.t,sol.y.T):out[float(t)]=v.copy()
        y=sol.y[:,-1]
    return np.asarray([out[float(t)] for t in times])

def component(parameters,run,population,full=False):
    x=np.array([parameters[k] for k in PARAMETERS]) if isinstance(parameters,dict) else np.asarray(parameters)
    s=space(*run['counts']);program=np.asarray(run['program']);initial=s.initial(run['initial']);prep=run['initial'].get('preparation_temperature_K') if isinstance(run['initial'],dict) else None
    if prep is not None:
        selected=initial>0;energy=s.obs[:,2:]@x[[2,3,8]]/(.00831446261815324*prep);initial[selected]*=np.exp(energy[selected]-energy[selected].max());initial/=initial.sum()
    y=np.r_[initial,0.,0.,program[0,1]+x[13]]
    def rates(T):
        v=s.rates(x,T)
        if population:
            shifts=np.where(s.feat[:,0]==0,x[16],np.where(s.feat[:,0]==1,x[17],x[18]));v*=np.exp(-shifts/(.00831446261815324*T))
        return v
    def rhs(t,y):
        T=y[-1];v=rates(T);p=y[:s.n]
        # Matrix-free action of the independently enumerated sparse generator.
        # Sum incoming and outgoing event flux separately, retaining self-events.
        # Rebuilding and sorting COO/CSR matrices at every ODE call is unnecessary.
        flow=v*p[s.i]
        dp=np.bincount(s.j,weights=flow,minlength=s.n)-np.bincount(s.i,weights=flow,minlength=s.n)
        CO=np.dot(v[s.des],p[s.i[s.des]])/12.;CO2=np.dot(v[s.ox],p[s.i[s.ox]])/12.
        return np.r_[dp,(CO-y[-3])/x[10],(CO2-y[-2])/x[11],(np.interp(t,program[:,0],program[:,1])+x[13]-T)/x[12]]
    ys=integrate(rhs,y,program,run['times']);p=ys[:,:s.n]
    coverage=p@s.obs[:,:2]/12.;yield2=run['counts'][1]/12.-coverage[:,1]
    yield1=run['counts'][0]/12.-coverage[:,0]-yield2
    result=dict(q=ys[:,-3:-1],sample_temperature=ys[:,-1],coverages=coverage,yields=np.c_[yield1,yield2])
    if full:
        flux=[]
        for T,pi in zip(ys[:,-1],p):
            v=rates(T);flux.append([np.dot(v[s.des],pi[s.i[s.des]])/12.,np.dot(v[s.ox],pi[s.i[s.ox]])/12.])
        result.update(pair_fractions=p@s.obs[:,2:]/24.,flux=np.asarray(flux),mass_error=float(np.max(abs(p.sum(axis=1)-1))))
    return result

def experiment(parameters,run,full=False):
    x=np.array([parameters[k] for k in PARAMETERS]) if isinstance(parameters,dict) else np.asarray(parameters)
    a=component(x,run,0,full);b=component(x,run,1,full);f=x[15]
    out={k:(1-f)*v+f*b[k] for k,v in a.items() if k!='mass_error'}
    if full:out['mass_error']=max(a['mass_error'],b['mass_error'])
    return out

def thermal_prediction(x,cal):
    program=np.asarray(cal['program']);y=np.array([program[0,1]+x[13]])
    rhs=lambda t,y: (np.interp(t,program[:,0],program[:,1])+x[13]-y)/x[12]
    return integrate(rhs,y,program,cal['times'],rtol=2e-11,atol=1e-10)[:,0]

def deviance(parameters,nuisance,data):
    x=np.array([parameters[k] for k in PARAMETERS]);total=0.;parts={}
    cal=data['thermal_calibration'];r=(thermal_prediction(x,cal)-cal['sample_temperature_K'])/cal['sigma_K'];parts['thermal']=float(r@r)
    cal=data['CO_pulse'];t=np.asarray(cal['times']);r=(cal['initial_current_pA']*np.exp(-t/x[10])-cal['current_28_pA'])/cal['sigma_pA'];parts['CO_pulse']=float(r@r)
    cal=data['CO2_pulse'];t=np.asarray(cal['times']);v=cal['initial_current_44_pA']*np.exp(-t/x[11]);r=(v[:,None]*[x[14],1.]-cal['current_28_44_pA'])/cal['sigma_pA'];parts['CO2_pulse']=float(np.sum(r*r))
    covariance=np.asarray(data['noise']['channel_covariance_pA2'])
    for run in data['runs']:
        pred=experiment(parameters,run);q=pred['q'];t=np.asarray(run['times']);g,b28,s28,b44,s44=nuisance[run['id']];u=t/t[-1]-.5
        model=np.c_[g*(1200*q[:,0]+x[14]*850*q[:,1])+b28+s28*u,g*850*q[:,1]+b44+s44*u,np.full(len(t),g*run['argon_reference_pA'])]
        error=(model-run['ion_current_pA']).ravel()
        temporal=np.exp(-abs(t[:,None]-t[None,:])/data['noise']['correlation_time_s'])
        matrix=np.kron(temporal,covariance);L=cholesky(matrix,lower=True)
        white=solve_triangular(L,error,lower=True);parts[run['id']]=float(white@white)
    return float(sum(parts.values())),parts

def forecast(parameters,case):
    rows=[]
    for run in case['forecasts']:
        p=experiment(parameters,run,full=True);snapshots=[]
        for j,t in enumerate(run['times']):
            snapshots.append(dict(time=t,sample_temperature=float(p['sample_temperature'][j]),coverages=p['coverages'][j].tolist(),gas_yields=p['yields'][j].tolist(),gas_rates=p['flux'][j].tolist(),pair_fractions=p['pair_fractions'][j].tolist()))
        rows.append(dict(id=run['id'],snapshots=snapshots))
    contrast=rows[0]['snapshots'][-1]['gas_yields'][1]-rows[1]['snapshots'][-1]['gas_yields'][1]
    return rows,float(contrast)
