"""Local sensitivity and nuisance-profiled information for paired experiments."""
import numpy as np
from numba import njit
from fast_model import trajectory

@njit
def responses(programs, pars, n=32, dt=.004, h=.003):
    indices=np.array([4,5,2,3,6,7,9,10])
    outputs=np.zeros((2,7)); bins=np.zeros((2,120)); jac=np.zeros((2,122,8))
    for run in range(2):
        outputs[run],bins[run]=trajectory(programs[run],pars,n,dt)
        for k in range(8):
            plus=pars.copy(); minus=pars.copy()
            plus[indices[k]]*=np.exp(h); minus[indices[k]]*=np.exp(-h)
            op,bp=trajectory(programs[run],plus,n,dt)
            om,bm=trajectory(programs[run],minus,n,dt)
            jac[run,:120,k]=(bp-bm)/(2*h)
            jac[run,120:,k]=(op[:2]-om[:2])/(2*h)
    return outputs,bins,jac

def information(bins, jac, cell):
    blocks=[]
    rho=cell['charge_bin_noise_rho']; sigma=cell['charge_bin_noise_sd']
    for run in range(2):
        j=np.zeros((122,12)); j[:,:8]=jac[run]
        j[:120,8+2*run]=bins[run]
        j[:120,9+2*run]=1.
        raw=j[:120].copy()
        j[0]=raw[0]/sigma
        j[1:120]=(raw[1:]-rho*raw[:-1])/(sigma*np.sqrt(1-rho*rho))
        j[120:]/=cell['assay_noise_sd']
        blocks.append(j)
    prior=np.zeros((10,12))
    sd=cell['nuisance_prior_sd']+[cell['gain_prior_sd'],cell['offset_prior_sd']]*2
    prior[:,2:]=np.diag(1/np.array(sd))
    full=np.vstack(blocks+[prior])
    # Project in whitened observation space, avoiding cancellation of large Gramians.
    interest=full[:,:2]; nuisance=full[:,2:]
    residual=interest-nuisance@np.linalg.lstsq(nuisance,interest,rcond=None)[0]
    return residual.T@residual

def evaluate(programs, parameters, cell, n=32, dt=.004, h=.003):
    records=[]
    for pars in parameters:
        outcomes,bins,jac=responses(programs,pars,n,dt,h)
        info=information(bins,jac,cell)
        records.append((outcomes,bins,info))
    return records

def summary(records):
    return np.array([min(np.linalg.eigvalsh(r[2])) for r in records]),np.array([r[0] for r in records])
