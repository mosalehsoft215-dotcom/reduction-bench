"""Additional independent-grid check on the final limiting scenario."""
import sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tests'))
from test_outputs import read_report,parse_programs,parse_predictions
from information import evaluate_pair
c=json.load((ROOT/'environment/data/cell.json').open());s=json.load((ROOT/'environment/data/scenarios.json').open())
report=read_report((Path(__file__).parent/'reference_output.json'));programs=parse_programs(report)
pred=parse_predictions(report)
index=list(s).index('slow_transport');p=np.array([s['slow_transport'][k] for k in c['scenario_parameter_order']])
start=time.monotonic();m,b,h=evaluate_pair(programs,p,c,n=120,step=.001)
val,vec=np.linalg.eigh(h);whiten=(vec/np.sqrt(val))@vec.T
r=dict(grid_cells=120,sensitivity_step=.001,scenario='slow_transport',seconds=time.monotonic()-start,
       score=float(val[0]),information=h.tolist(),outcomes=m.tolist(),
       max_outcome_error=float(abs(pred[index][0]-m).max()),max_bin_error=float(abs(pred[index][1]-b).max()),
       scaled_information_error=float(np.linalg.norm(whiten@(pred[index][2]-h)@whiten,ord=2)))
assert r['max_outcome_error']<.003 and r['max_bin_error']<.00010 and r['scaled_information_error']<.05
assert val[0]>=9
(Path(__file__).parent/'reference_grid_check.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
