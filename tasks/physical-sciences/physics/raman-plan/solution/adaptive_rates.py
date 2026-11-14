"""Exact nonlinear optical rates for the public persistent-state bank."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import optics as o

def inputs():
    rows=list(csv.DictReader((o.DATA/'settings.csv').open()))
    catalog=np.asarray([[int(r['target']),int(r['setting_id']),int(r['laser_nm']),float(r['stage_azimuth_deg']),float(r['tilt_x_deg']),float(r['tilt_y_deg']),float(r['half_wave_plate_deg']),int(r['crossed'])] for r in rows])
    bank=json.loads((o.DATA/'optical_states.json').read_text())['branches']
    return catalog,bank,json.loads((o.DATA/'design_limits.json').read_text())

def rates(frames=None):
    original_frames=o.H
    if frames is not None:o.H=frames
    catalog,bank,limits=inputs();base={k:[dict(r) for r in o.STACK[k]] for k in ('sample_layers','reference_layers')};allmeans=[]
    try:
        for target,index,*setting in catalog:
            bybranch=[]
            for branch in bank:
                for crystal in ('sample','reference'):
                    o.STACK[crystal+'_layers']=[dict(r,thickness_nm=r['thickness_nm']*branch['common_layer_thickness_scale']) for r in base[crystal+'_layers']]
                command=setting.copy();command[0]=int(command[0]);command[1]+=branch['azimuth_registration_deg'];command[2]+=branch['tilt_x_registration_deg'];command[3]+=branch['tilt_y_registration_deg']
                value=o.cresp(tuple(command))[None,:].repeat(12,axis=0) if target==0 else o.sresp(tuple(command))
                bybranch.append(value)
            value=np.asarray(bybranch).transpose(1,0,2).reshape(72,-1)
            if target==0:value=np.pad(value,((0,0),(0,1)),constant_values=np.nan)
            allmeans.append(value)
    finally:
        o.STACK.update(base);o.H=original_frames
    return catalog,np.asarray(allmeans),limits

def signature():
    h=hashlib.sha256()
    for p in sorted(o.DATA.iterdir()):
        if p.suffix in ('.csv','.json','.npz','.py'):h.update(p.name.encode());h.update(p.read_bytes())
    for name in ('adaptive_rates.py','optics.py','coating_fields.py','physical_model.py','count_kernels.py'):h.update(Path(__file__).with_name(name).read_bytes())
    return h.hexdigest()

def load(cache):
    key=signature();cache=Path(cache)
    if cache.exists():
        with np.load(cache,allow_pickle=False) as z:
            if str(z['signature'])==key:return z['catalog'],z['means'],json.loads(str(z['limits']))
    catalog,means,limits=rates();cache.parent.mkdir(parents=True,exist_ok=True)
    np.savez(cache,catalog=catalog,means=means,limits=json.dumps(limits),signature=key)
    return catalog,means,limits

def leg(a,b,cfg):
    m=cfg['motion_cost_s'];az=abs(a[3]-b[3])%360;hw=abs(a[6]-b[6])%90
    return m['azimuth_per_deg']*min(az,360-az)+m['tilt_per_deg_each_axis']*(abs(a[4]-b[4])+abs(a[5]-b[5]))+m['half_wave_plate_per_deg_mod90']*min(hw,90-hw)+m['analyzer_toggle']*(a[7]!=b[7])+m['laser_switch']*(a[2]!=b[2])+cfg['target_switch_s']*(a[0]!=b[0])
