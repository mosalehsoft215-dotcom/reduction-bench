"""Generate only the added calibrated delivery experiment, with a fixed seed."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
BUNDLE=HERE.parent.parent
sys.path.insert(0,str(BUNDLE/'solution'))
from delivery import calibration_prob,NAMES

TRUTH=np.array([1.036,.024,-.48,.13,.972,-.018,-.39,.105])
def generate():
    model=dict(radius_mm=.60,gradient_T_per_m_per_A=.8,effective_length_m=.010)
    rng=np.random.default_rng(8041026)
    runs=[]
    for j in (1,2):
        for current in (420.,880.):
            for x in (-.49,-.36,-.22,-.07,.09,.24,.39,.51):
                run=dict(id=f'L{j}-I{int(current)}-x{x:+.2f}',lens=j,current_A=current,
                         energy_MeV=252.,input_mean_mm=x,input_sigma_mm=.004,
                         input_angle_sigma_mrad=.008,drift_m=.50,resolution_mm=.008)
                center=x-.50*.299792458*.8*.010*current/(np.sqrt(252*(252+1.0219979))/1000)*x
                run['edges_mm']=np.linspace(center-.60,center+.60,97).tolist()
                pr=calibration_prob(TRUTH,run,model,96); pr=pr/pr.sum()
                run['counts']=rng.multinomial(400,pr).tolist()
                runs.append(run)
    bounds=[[.90,1.15],[-.06,.06],[-.65,-.20],[.03,.24],
            [.85,1.10],[-.06,.06],[-.62,-.18],[.03,.23]]
    case=dict(model=model,parameter_order=NAMES,parameter_bounds=bounds,likelihood_delta=20.)
    data=BUNDLE/'environment'/'data'
    (data/'delivery_calibration.json').write_bytes((json.dumps(dict(runs=runs),indent=2)+'\n').encode())
    (data/'delivery_case.json').write_bytes((json.dumps(case,indent=2)+'\n').encode())
    record=dict(seed=8041026,parameters=dict(zip(NAMES,TRUTH.tolist())),
                model='nonlinear paraxial phase sheets; independent multinomial probe bins',
                files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [data/'delivery_case.json',data/'delivery_calibration.json']})
    (HERE/'delivery_generation_record.json').write_bytes((json.dumps(record,indent=2)+'\n').encode())
if __name__=='__main__': generate()
