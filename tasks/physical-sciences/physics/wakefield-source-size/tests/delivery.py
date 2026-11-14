"""Calibrated nonlinear phase-sheet likelihood and bounded parameter fitting."""
import json
from pathlib import Path
import numpy as np
from functools import lru_cache
from scipy.special import ndtr
from scipy.optimize import least_squares

ME = .51099895
NAMES = ['scale_1', 'center_1_mm', 'cubic_1', 'quintic_1',
         'scale_2', 'center_2_mm', 'cubic_2', 'quintic_2']

@lru_cache(None)
def hermite_rule(order):
    t,w=np.polynomial.hermite.hermgauss(order)
    return t,w

def strength(current, energy, scale, model):
    p = np.sqrt(energy * (energy + 2 * ME)) / 1000
    return .299792458 * model['gradient_T_per_m_per_A'] * model['effective_length_m'] * current * scale / p

def field(x, q, j, model):
    y = x - q[4*j+1]
    r = model['radius_mm']
    return y + q[4*j+2]*y**3/r**2 + q[4*j+3]*y**5/r**4

def calibration_prob(q, run, model, order=48):
    t, w = hermite_rule(order)
    x = run['input_mean_mm'] + np.sqrt(2)*run['input_sigma_mm']*t
    j = run['lens'] - 1
    k = strength(run['current_A'], run['energy_MeV'], q[4*j], model)
    y = x - run['drift_m']*k*field(x, q, j, model)
    sig = np.hypot(run['drift_m']*run['input_angle_sigma_mrad'], run['resolution_mm'])
    edges = np.asarray(run['edges_mm'])
    cdf = ndtr((edges[:,None]-y[None,:])/sig) @ (w/np.sqrt(np.pi))
    # Overflow bins are part of the multinomial observation.
    return np.maximum(np.diff(np.r_[0., cdf, 1.]), 1e-300)

def deviance(q, runs, model, order=48):
    total = 0.
    for r in runs:
        n = np.asarray(r['counts'], float)
        expected = n.sum()*calibration_prob(q, r, model, order)
        use = n > 0
        total += 2*np.sum(n[use]*np.log(n[use]/expected[use]))
    return float(total)

def deviance_residual(q, runs, model):
    out = []
    for r in runs:
        n = np.asarray(r['counts'], float)
        mu = n.sum()*calibration_prob(q, r, model)
        v = mu-n
        use = n>0
        v[use] += n[use]*np.log(n[use]/mu[use])
        out.extend(np.sign(n-mu)*np.sqrt(np.maximum(2*v, 0)))
    return np.asarray(out)

def fit(runs, case):
    lo, hi = np.asarray(case['parameter_bounds']).T
    center = (lo+hi)/2
    fits = [least_squares(deviance_residual, center, args=(runs, case['model']),
                         bounds=(lo,hi), x_scale=hi-lo, ftol=1e-11, xtol=1e-11, gtol=1e-7)]
    # Starts differ in the correlated cubic/quintic directions.
    for f in (.25,.75):
        start = center.copy(); start[[2,3,6,7]] = lo[[2,3,6,7]]+f*(hi-lo)[[2,3,6,7]]
        fits.append(least_squares(deviance_residual,start,args=(runs,case['model']),
                    bounds=(lo,hi), x_scale=hi-lo, ftol=1e-11,xtol=1e-11,gtol=1e-7))
    best = min(fits, key=lambda f: f.cost)
    return best.x, deviance(best.x,runs,case['model']), best
