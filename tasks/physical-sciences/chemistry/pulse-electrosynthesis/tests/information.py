"""Independent finite-difference likelihood curvature and Schur complement."""
import numpy as np
from physics import trajectory

def evaluate_pair(programs, pars, cell, n=80, step=.002):
    metrics=[]; signals=[]; gram=np.zeros((12,12))
    correlation=cell['charge_bin_noise_rho']
    # Cholesky factor of inverse AR(1) covariance; use explicit whitening matrix.
    whitening=np.eye(120)
    for i in range(1,120):
        whitening[i,i-1]=-correlation
    whitening[1:]/=np.sqrt(1-correlation**2)
    whitening/=cell['charge_bin_noise_sd']
    for experiment,program in enumerate(programs):
        metric,signal=trajectory(program,pars,n=n)
        metrics.append(metric);signals.append(signal)
        sensitivity=np.zeros((122,12))
        for column,index in enumerate(cell['sensitivity_parameter_indices']):
            lower=np.array(pars,dtype=float);upper=lower.copy()
            lower[index]*=np.exp(-step);upper[index]*=np.exp(step)
            ml,sl=trajectory(program,lower,n=n)
            mu,su=trajectory(program,upper,n=n)
            sensitivity[:,column]=(np.r_[su,mu[:2]]-np.r_[sl,ml[:2]])/(2*step)
        sensitivity[:120,8+2*experiment]=signal
        sensitivity[:120,9+2*experiment]=1
        sensitivity[:120]=whitening@sensitivity[:120]
        sensitivity[120:]/=cell['assay_noise_sd']
        gram+=sensitivity.T@sensitivity
    sd=np.array(cell['nuisance_prior_sd']+[cell['gain_prior_sd'],cell['offset_prior_sd']]*2)
    gram[2:,2:]+=np.diag(1/sd**2)
    # Diagonal scaling of nuisance coordinates keeps offset units from harming the solve.
    scale=1/np.sqrt(np.diag(gram)[2:])
    cross=gram[:2,2:]*scale
    block=gram[2:,2:]*scale[:,None]*scale[None,:]
    effective=gram[:2,:2]-cross@np.linalg.solve(block,cross.T)
    effective=(effective+effective.T)/2
    assert np.linalg.eigvalsh(effective)[0]>0
    return np.array(metrics),np.array(signals),effective
