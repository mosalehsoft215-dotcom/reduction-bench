"""Content-checked private caches contain physical kernels, never policies."""
from pathlib import Path
import hashlib,json
import numpy as np
import optics as o
from adaptive_rates import load,signature
from physical_model import derive
from count_kernels import fourier_kernels

def comparison_bank(cache):
    catalog,rates,G,limits=derive()
    rows=[];groups=[]
    for i,row in enumerate(catalog):
        for channel in range(3 if row[0]==0 else 4):
            for n in (2,4,8,12,20,32):
                ordered=np.unique(np.round(rates[i,:,:,channel].ravel()*n*o.CFG['exposure_unit_s'],9))
                thresholds=np.unique(np.maximum(0,np.rint((ordered[:-1]+ordered[1:])/2).astype(int)))
                if len(thresholds)>16:thresholds=np.unique(thresholds[np.linspace(0,len(thresholds)-1,16).astype(int)])
                if not len(thresholds):continue
                start=len(rows);rows.extend([[i,channel,n,int(t)] for t in thresholds]);groups.append((i,channel,n,thresholds,start))
    actions=np.asarray(rows,int)
    key=hashlib.sha256(signature().encode()+actions.tobytes()+Path(__file__).read_bytes()).hexdigest();path=Path(cache).with_suffix('.comparisons.npz')
    if path.exists():
        with np.load(path,allow_pickle=False) as z:
            if str(z['signature'])==key:return catalog,actions,z['zero'],z['one'],rates,G,limits
    zero=np.empty((len(actions),rates.shape[1],3,3));one=np.empty_like(zero)
    for group,(i,ch,n,taus,start) in enumerate(groups):
        lo,hi=fourier_kernels(G[i],rates[i,:,:,ch],n*o.CFG['exposure_unit_s'],taus);end=start+len(taus);zero[start:end]=lo;one[start:end]=hi
        if group%128==0:print('Physical comparison groups',group,'/',len(groups),flush=True)
    np.savez(path,signature=key,zero=zero,one=one)
    return catalog,actions,zero,one,rates,G,limits
