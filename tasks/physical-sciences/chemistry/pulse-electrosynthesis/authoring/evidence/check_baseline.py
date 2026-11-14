import sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tests'))
from information import evaluate_pair
c=json.load((ROOT/'environment/data/cell.json').open());s=json.load((ROOT/'environment/data/scenarios.json').open())
p=np.array([s['slow_transport'][k] for k in c['scenario_parameter_order']])
entry=json.load((Path(__file__).parent/'calibration_results.json').open())[-1]
programs=[np.array(a) for a in entry['programs']]
start=time.monotonic();m,q,h=evaluate_pair(programs,p,c,n=80)
result=dict(family=entry['family'],scenario='slow_transport',information=h.tolist(),score=float(np.linalg.eigvalsh(h)[0]),min_yield=float(m[:,0].min()),max_charge=float(m[:,2].max()),seconds=time.monotonic()-start)
(Path(__file__).parent/'baseline_independent.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
