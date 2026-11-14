"""Recompute every physical path, conditional risk and posterior from a tree."""
import math
import numpy as np
from adaptive_rates import leg
from count_kernels import fourier_kernels
from physical_model import predictions

def audit(policy,catalog,rates,G,cfg,limits):
    index={int(row[1]):i for i,row in enumerate(catalog)};cache={};leaves=[];max_depth=0
    groups={}
    def collect(node):
        if 'domain' in node:return
        setting,ch,n,tau=node['action'];groups.setdefault((setting,ch,n),set()).add(tau)
        for child in node['children']:collect(child)
    collect(policy)
    for (setting,ch,n),thresholds in groups.items():
        ts=sorted(thresholds);i=index[setting]
        lo,hi=fourier_kernels(G[i],rates[i,:,:,ch],n*cfg['exposure_unit_s'],ts)
        for j,t in enumerate(ts):cache[(setting,ch,n,t)]=(lo[j],hi[j])
    def visit(node,w,path,previous,clock):
        nonlocal max_depth
        assert len(path)<=limits['max_policy_depth'];max_depth=max(max_depth,len(path))
        if 'domain' in node:
            total=math.fsum(w.ravel());marginal=float(np.clip(total/rates.shape[1],0,1));posterior=None if marginal<=1e-6 else np.clip(w.ravel()/total,0,1).tolist()
            leaves.append(dict(path=path,domain=node['domain'],conditional_probabilities=np.clip(w.sum(axis=1),0,1).tolist(),lost_probabilities=np.clip(w[:,2],0,1).tolist(),posterior=posterior,marginal_probability=marginal,end_time_s=float(clock)));return
        setting,ch,n,tau=node['action'];i=index[setting];row=catalog[i];key=tuple(node['action'])
        if key not in cache:
            a,b=fourier_kernels(G[i],rates[i,:,:,ch],n*cfg['exposure_unit_s'],[tau]);cache[key]=(a[0],b[0])
        end=clock+(leg(catalog[previous],row,cfg) if previous is not None else 0)+n*cfg['exposure_unit_s']+cfg['readout_s_per_spectrum'];assert end<=limits['time_budget_s']+1e-9
        for bit,k in enumerate(cache[key]):visit(node['children'][bit],np.einsum('sf,sfg->sg',w,k),path+str(bit),i,end)
    initial=np.zeros((rates.shape[1],3));initial[:,0]=1.;visit(policy,initial,'',None,0.)
    mass=np.asarray([x['conditional_probabilities'] for x in leaves]);loss=np.asarray([x['lost_probabilities'] for x in leaves]);labels=np.asarray([x['domain'] for x in leaves]);truth=np.repeat(np.arange(12),rates.shape[1]//12);times=np.asarray([x['end_time_s'] for x in leaves]);assert np.max(abs(mass.sum(axis=0)-1))<1e-10
    errors=np.clip((mass*(labels[:,None]!=truth[None,:])).sum(axis=0),0,1);lost=np.clip(loss.sum(axis=0),0,1)
    return dict(optical_predictions=predictions(catalog,rates,G),errors=errors.tolist(),mean_error=float(errors.mean()),worst_error=float(errors.max()),memory_loss=lost.tolist(),mean_memory_loss=float(lost.mean()),worst_memory_loss=float(lost.max()),expected_time_s=(mass*times[:,None]).sum(axis=0).tolist(),max_path_time_s=float(times.max()),leaves=leaves)
