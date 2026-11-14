"""Coherent plane-wave transmission, using reflection recursion.

The amplitudes propagated here are tangential E fields. Returned p amplitude
is converted to the full electric-field amplitude in the active substrate.
No selected schedule or optical-response table is stored.
"""
import numpy as np
from functools import lru_cache

def transmission(ci, substrate_index, wavelength_nm, layers):
    key=tuple((r['thickness_nm'],r['cauchy_A'],r['cauchy_B_um2']) for r in layers)
    return _cached(float(ci),float(substrate_index),float(wavelength_nm),key)

@lru_cache(maxsize=65536)
def _cached(ci, substrate_index, wavelength_nm, key):
    layers=[dict(thickness_nm=d,cauchy_A=a,cauchy_B_um2=b) for d,a,b in key]
    n=np.asarray([1.]+[r['cauchy_A']+r['cauchy_B_um2']/(wavelength_nm/1000.)**2 for r in layers]+[substrate_index])
    c=np.sqrt(1-(1-ci*ci)/n**2)
    phase=np.asarray([0.]+[2*np.pi*n[j+1]*c[j+1]*r['thickness_nm']/wavelength_nm for j,r in enumerate(layers)]+[0.])
    answer=[]
    for admittance in (n*c,n/c):
        reflection=0j; transmitted=1+0j
        for j in range(len(n)-2,-1,-1):
            r=(admittance[j]-admittance[j+1])/(admittance[j]+admittance[j+1])
            t=2*admittance[j]/(admittance[j]+admittance[j+1])
            propagation=np.exp(1j*phase[j+1]);denominator=1+r*reflection*propagation**2
            transmitted=t*transmitted*propagation/denominator
            reflection=(r+reflection*propagation**2)/denominator
        answer.append(transmitted)
    return answer[0],answer[1]*ci/c[-1]
