"""Derive coupled optical and source-state-dependent switching responses."""
import json
import numpy as np
import optics as o
from adaptive_rates import rates,inputs

def strengths(catalog,frames):
    model=json.loads((o.DATA/'photodynamics.json').read_text());A=np.asarray(model['local_strength_tensor'])
    _,branches,_=inputs();out=np.zeros((len(catalog),72))
    original={key:[dict(layer) for layer in o.STACK[key]] for key in ('sample_layers','reference_layers')}
    try:
        for i,row in enumerate(catalog):
            if row[0]==0:continue
            wavelength=int(row[2]);n=o.index['sample'][wavelength];mirror=o.reflect[wavelength]
            for b,branch in enumerate(branches):
                for key,layers in original.items():o.STACK[key]=[dict(x,thickness_nm=x['thickness_nm']*branch['common_layer_thickness_scale']) for x in layers]
                command=row[2:].copy();command[0]=int(command[0]);command[1]+=branch['azimuth_registration_deg'];command[2]+=branch['tilt_x_registration_deg'];command[3]+=branch['tilt_y_registration_deg']
                pose,angle,_,_,_=o.geom(command)
                excitation=np.array([np.sqrt(mirror['R_p'])*np.cos(angle),np.sqrt(mirror['R_s'])*np.exp(1j*np.deg2rad(mirror['phase_s_minus_p_deg']))*np.sin(angle),0.])
                E=o.inside(excitation,pose[:,2],n,wavelength,'sample')
                for d,U in enumerate(frames):
                    frame=pose@U;out[i,6*d+b]=np.real(np.vdot(E,(frame@A@frame.T)@E))/(2/(1+n))**2
    finally:o.STACK.update(original)
    return out

def derive():
    catalog,mu,limits=rates();model=json.loads((o.DATA/'photodynamics.json').read_text())
    families=model['kinetic_families'];unit=o.CFG['exposure_unit_s'];c=model['specimen_collection_scale'];s=model['specimen_signal_scale']
    photons=np.repeat(mu[:,:,None,:]/unit,3,axis=2);photons=np.repeat(photons,2,axis=1)
    generators=np.zeros((len(catalog),144,3,3));j0=strengths(catalog,o.H)
    for k,family in enumerate(families):
        V=np.asarray(family['switched_parent_rotation']);frames=[V@U for U in o.H]
        _,switched,_=rates(frames);j1=strengths(catalog,frames)
        for i,row in enumerate(catalog):
            if int(row[0])==0:continue
            wavelength=str(int(row[2]));B=o.sb[int(row[2])][None,:]
            bright=c*(B+s*(mu[i]-B));dark=c*(B+model['dark_signal_fraction'][wavelength]*s*(switched[i]-B))
            mixed=np.tile(bright.reshape(12,6,4).mean(axis=0),(12,1))
            photons[i,k::2]=np.stack([bright,dark,mixed],axis=1)/unit
            co=family['transition_coefficients_per_s_per_strength'][wavelength]
            for a,b,name,j in ((0,1,'bright_to_dark',j0),(0,2,'bright_to_scrambled',j0),(1,0,'dark_to_bright',j1),(1,2,'dark_to_scrambled',j1)):
                generators[i,k::2,a,b]=j[i]*co[name]
            for a in (0,1):generators[i,k::2,a,a]=-generators[i,k::2,a].sum(axis=-1)
    return catalog,photons,generators,limits

def predictions(catalog, photons, generators):
    queries=json.loads((o.DATA/'response_queries.json').read_text())['queries'];positions={int(row[1]):i for i,row in enumerate(catalog)};result=[]
    for q in queries:
        i=positions[q['setting_id']];s=q['initial_state'];ch=3 if int(catalog[i,0])==0 else 4
        result.append(dict(query_id=q['query_id'],photon_rates_per_s=photons[i,s,:,:ch].tolist(),generator_per_s=generators[i,s].tolist()))
    return result
