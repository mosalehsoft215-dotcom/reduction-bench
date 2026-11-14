"""Coating fields from a direct boundary-value system, no recursion.

Every film carries forward/backward tangential E amplitudes. Air has a known
unit incoming wave; the active half-space has only an outgoing wave.
"""
import numpy as np
from functools import lru_cache

def field_amplitudes(cos_air,index_substrate,wavelength_nm,layers):
    key=tuple((r['thickness_nm'],r['cauchy_A'],r['cauchy_B_um2']) for r in layers)
    return _boundary_solution(float(cos_air),float(index_substrate),float(wavelength_nm),key)

@lru_cache(maxsize=65536)
def _boundary_solution(cos_air,index_substrate,wavelength_nm,key):
    layers=[dict(thickness_nm=d,cauchy_A=a,cauchy_B_um2=b) for d,a,b in key]
    if not layers:
        cs=np.sqrt(1-(1-cos_air*cos_air)/index_substrate**2)
        return 2*cos_air/(cos_air+index_substrate*cs),2*cos_air/(index_substrate*cos_air+cs)
    indices=[1.]+[float(r['cauchy_A'])+float(r['cauchy_B_um2'])*(1000./wavelength_nm)**2 for r in layers]+[index_substrate]
    cosines=np.sqrt(1-(1-cos_air*cos_air)/np.square(indices));nfilm=len(layers);count=2*nfilm+2
    phases=[np.exp(2j*np.pi*indices[k+1]*cosines[k+1]*r['thickness_nm']/wavelength_nm) for k,r in enumerate(layers)]
    output=[]
    for p in (False,True):
        y=np.asarray(indices)/cosines if p else np.asarray(indices)*cosines
        matrix=np.zeros((count,count),complex);rhs=np.zeros(count,complex)
        for interface in range(nfilm+1):
            row=2*interface
            if interface==0:
                matrix[row,0]=1;matrix[row+1,0]=-y[0];rhs[row]=-1;rhs[row+1]=-y[0]
            else:
                a=2*interface-1;b=a+1;phase=phases[interface-1]
                matrix[row,a]=phase;matrix[row,b]=1/phase
                matrix[row+1,a]=y[interface]*phase;matrix[row+1,b]=-y[interface]/phase
            if interface==nfilm:
                matrix[row,-1]=-1;matrix[row+1,-1]=-y[-1]
            else:
                a=2*interface+1;b=a+1
                matrix[row,a]=-1;matrix[row,b]=-1
                matrix[row+1,a]=-y[interface+1];matrix[row+1,b]=y[interface+1]
        amplitude=np.linalg.solve(matrix,rhs)[-1]
        output.append(amplitude*cos_air/cosines[-1] if p else amplitude)
    return tuple(output)
