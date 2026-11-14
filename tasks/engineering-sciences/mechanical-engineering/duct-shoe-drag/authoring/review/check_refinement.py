"""Check next-cycle row peaks under four linear substeps, preserving brake-command ordering."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
p=ROOT/'solution/shoe_reduction.py'
spec=importlib.util.spec_from_file_location('forward',p);r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
truth=json.loads((ROOT/'authoring/evidence/generation.json').read_text())['physical_inputs']
results=[]
for name,physical in truth.items():
    folder=ROOT/('environment/data/startups' if name.startswith('W') else 'tests/cycles')/name
    d=json.loads((folder/'assembly.json').read_text());a=np.array(d['next_cycle']['rows'])
    lock=[d['next_cycle']['columns'].index(k) for k in ['lock_P','lock_S']]
    refined=[a[0]]
    for prev,current in zip(a[:-1],a[1:]):
        for sub in range(1,5):
            q=prev+(current-prev)*sub/4;q[lock]=current[lock]
            refined.append(q)
    fine=copy.deepcopy(d);fine['next_cycle']['rows']=np.array(refined).tolist()
    args=(np.array(physical['top']),np.array(physical['carrier']),physical['force_span'])
    coarse=r.replay(d,*args);fine_trace=r.replay(fine,*args)[::4]
    coarse_peak=float(np.max(abs(coarse[:,0]-coarse[:,1])))
    fine_peak=float(np.max(abs(fine_trace[:,0]-fine_trace[:,1])))
    row={'campaign':name,'coarse_peak_mm':coarse_peak,'refined_peak_at_original_rows_mm':fine_peak,
         'difference_mm':abs(coarse_peak-fine_peak),'coordinate_max_difference_mm':float(np.max(abs(coarse-fine_trace)))}
    print(row,flush=True);results.append(row)
    assert row['difference_mm']<.035,row
(ROOT/'authoring/evidence/refinement.json').write_text(json.dumps({'substeps':4,'rows':results,
    'reference_source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2)+'\n',encoding='utf-8',newline='\n')
