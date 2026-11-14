"""Fast strict-contract, physical-gate, and information-error regressions."""
import sys,json,copy,tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tests'))
import test_outputs as t

passed=[]
def check(name,fun,expected=True):
    try:fun();ok=True
    except Exception:ok=False
    assert ok==expected,name
    passed.append(name)

fake={'programs':{a:[{'duration_s':6.,'potential_V':-1.1}] for a in ['A','B']},'predictions':{}}
for name in t.SCENARIOS:
    fake['predictions'][name]={'experiments':{a:{'outcomes':dict(zip(t.CELL['metrics_order'],[.25,.2,.69,.5,.05,.5,.5])),
        'charge_bins':[.69/120]*120} for a in ['A','B']},'profile_information':[[10.,0.],[0.,12.]]}
check('valid programs',lambda:t.parse_programs(fake));check('valid predictions',lambda:t.parse_predictions(fake))
for field,values in [('duration_s',[True,'6',None,[],0,-1,.019,6.01,float('nan'),float('inf')]),
                     ('potential_V',[True,'-1.1',None,[],0,-1.31,float('nan'),float('inf')])]:
    for i,v in enumerate(values):
        r=copy.deepcopy(fake);r['programs']['A'][0][field]=v
        check(f'{field} invalid {i}',lambda:t.parse_programs(r),False)
for key in ['A','B']:
    r=copy.deepcopy(fake);del r['programs'][key];check('missing '+key,lambda:t.parse_programs(r),False)
for count in [0,49]:
    r=copy.deepcopy(fake);r['programs']['A']=[{'duration_s':6/max(1,count),'potential_V':-1.1}]*count
    check(f'count {count}',lambda:t.parse_programs(r),False)
for residual in [-5e-7,5e-7]:
    r=copy.deepcopy(fake);r['programs']['A'][0]['duration_s']+=residual
    check('duration residue '+str(residual),lambda:t.parse_programs(r))
    assert abs(t.parse_programs(r)[0][:,0].sum()-6)<1e-12
for name in t.SCENARIOS:
    r=copy.deepcopy(fake);del r['predictions'][name];check('missing scenario '+name,lambda:t.parse_predictions(r),False)
for metric in t.CELL['metrics_order']:
    for value in [True,'0.2',None,float('nan'),float('inf')]:
        r=copy.deepcopy(fake);r['predictions']['nominal']['experiments']['A']['outcomes'][metric]=value
        check(f'metric {metric} {value}',lambda:t.parse_predictions(r),False)
for value in [[],[0]*119,[0]*121,[True]*120,[-.01]*120]:
    r=copy.deepcopy(fake);r['predictions']['nominal']['experiments']['A']['charge_bins']=value
    check('invalid bins '+str(len(passed)),lambda:t.parse_predictions(r),False)
r=copy.deepcopy(fake);r['predictions']['nominal']['experiments']['A']['charge_bins']=[.70/120]*120
check('inconsistent total charge',lambda:t.parse_predictions(r),False)
for value in [[[10,0]],[[10,0],[0]],[[10,1],[0,12]],[[0,0],[0,12]],[[10,0],[0,-12]],[[True,0],[0,12]]]:
    r=copy.deepcopy(fake);r['predictions']['nominal']['profile_information']=value
    check('invalid H '+str(len(passed)),lambda:t.parse_predictions(r),False)
with tempfile.TemporaryDirectory() as folder:
    path=Path(folder)/'protocol.json'
    check('missing file',lambda:t.read_report(path),False)
    for text in ['{} garbage','{"a":1,"a":2}','{"a":NaN}','{"a":Infinity}','[]','x'*262144]:
        path.write_text(text);check('invalid JSON '+str(len(passed)),lambda:t.read_report(path),False)
valid=t.parse_predictions(fake)
check('performance gate valid',lambda:t.check_performance(valid))
for name,index,column,value in [('low yield',0,0,.20),('over charge',0,2,.80)]:
    a=copy.deepcopy(valid);a[0][index][0,column]=value
    check(name,lambda:t.check_performance(a),False)
a=copy.deepcopy(valid);a[0][2][0,0]=.1
check('unidentifiable direction',lambda:t.check_performance(a),False)
check('matching predictions',lambda:t.check_predictions(valid,valid))
for kind in ['invented matrix','wrong bins','wrong outcomes','weak direction hidden']:
    a=copy.deepcopy(valid)
    if kind=='invented matrix':a[0][2][:]*=10
    if kind=='wrong bins':a[0][1][0,0]+=.01
    if kind=='wrong outcomes':a[0][0][0,0]+=.1
    if kind=='weak direction hidden':a[0][2][0,0]*=.8
    check(kind,lambda:t.check_predictions(a,valid),False)
(Path(__file__).parent/'contract_results.json').write_text(json.dumps({'passed':len(passed),'cases':passed},indent=2))
print(f'{len(passed)} contract and scoring regressions passed',flush=True)
