"""Independent whole-policy probabilities with persistent material status."""
import math
from functools import lru_cache
import numpy as np
from optical_check import TRUTH
from independent_dynamics import independent_rates
from independent_contour_kernel import count_comparator_kernel
from policy_validation import validate

@lru_cache(maxsize=1)
def physical_inputs():
    return independent_rates()

@lru_cache(maxsize=1024)
def action_kernel(key):
    catalog,rates,G,cfg=physical_inputs()
    setting,ch,n,tau=key
    i=next(i for i,row in enumerate(catalog) if int(row[1])==setting)
    return count_comparator_kernel(G[i],rates[i,:,:,ch],n*cfg['exposure_unit_s'],tau)

def evaluate(policy):
    catalog,rates,G,cfg=physical_inputs();limits=TRUTH['limits'];structure=validate(policy,catalog,cfg,limits)
    clocks={row['path']:row['end_time_s'] for row in structure['leaves']};index={int(row[1]):i for i,row in enumerate(catalog)};cache={};rows=[]
    def walk(node,w,path):
        if 'domain' in node:
            conditional=w.sum(axis=1);total=math.fsum(w.ravel());marginal=total/144
            posterior=(w.ravel()/total).tolist() if total>0 else None
            rows.append(dict(path=path,domain=node['domain'],conditional_probabilities=conditional.tolist(),lost_probabilities=w[:,2].tolist(),posterior=posterior,marginal_probability=marginal,end_time_s=clocks[path]));return
        key=tuple(node['action'])
        if key not in cache:
            cache[key]=action_kernel(key)
        for bit,kernel in enumerate(cache[key]):walk(node['children'][bit],np.einsum('sf,sfg->sg',w,kernel),path+str(bit))
    initial=np.zeros((144,3));initial[:,0]=1.;walk(policy,initial,'')
    masses=np.asarray([x['conditional_probabilities'] for x in rows]);lost=np.asarray([x['lost_probabilities'] for x in rows]);labels=np.asarray([x['domain'] for x in rows]);truth=np.repeat(np.arange(12),12);times=np.asarray([x['end_time_s'] for x in rows])
    assert np.max(abs(masses.sum(axis=0)-1))<1e-10
    error=(masses*(labels[:,None]!=truth[None,:])).sum(axis=0);loss=lost.sum(axis=0)
    return dict(errors=error.tolist(),mean_error=float(error.mean()),worst_error=float(error.max()),memory_loss=loss.tolist(),mean_memory_loss=float(loss.mean()),worst_memory_loss=float(loss.max()),expected_time_s=(masses*times[:,None]).sum(axis=0).tolist(),max_path_time_s=structure['max_path_time_s'],max_depth=structure['max_depth'],leaves=rows)
