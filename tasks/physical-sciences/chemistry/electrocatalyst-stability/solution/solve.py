from pathlib import Path
import argparse,json,time
import numpy as np
from fit import fit
from model import impedance,steady_theta
from analyze import folds,equilibria,hopfs,orbit
from cycles import criticality,bistability

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--data',default='/app/data');parser.add_argument('--output',default='/app/output/stability.json');args=parser.parse_args()
    start=time.monotonic();data=Path(args.data)
    case=json.loads((data/'case.json').read_text());records=json.loads((data/'experiments.json').read_text())
    fits={m:fit(case,records,m) for m in case['mechanisms']}
    m=min(fits,key=lambda m:fits[m]['chi2']);p=np.array(fits[m]['parameters'])
    conditions={c['id']:c for c in case['conditions']}
    artifact={'mechanism':m,'parameters':dict(zip(case['parameter_order'],p.tolist())),
              'fit_chi_square':0.,'measurements':[],'folds':{},'queries':[]}
    chi=0.
    for r in records:
        cond=conditions[r['condition_id']];roots=equilibria(p,cond,m,r['command_V'])
        x=roots[0 if r['branch']=='lowest-theta' else -1]
        _,_,current=steady_theta(x[1],p,cond,m)
        freq=[s['frequency_Hz'] for s in r['spectrum']];z=impedance(x,p,cond,m,freq)
        chi+=((current-r['current_A_cm2'])/r['current_sigma_A_cm2'])**2
        chi+=sum(((zz.real-s['real_ohm_cm2'])/s['sigma_ohm_cm2'])**2+((zz.imag-s['imag_ohm_cm2'])/s['sigma_ohm_cm2'])**2 for zz,s in zip(z,r['spectrum']))
        artifact['measurements'].append({'id':r['id'],'state':x.tolist(),'current_A_cm2':float(current),'impedance_ohm_cm2':[[float(zz.real),float(zz.imag)] for zz in z]})
    artifact['fit_chi_square']=float(chi)
    for name,cond in conditions.items():artifact['folds'][name]=folds(p,cond,m)
    for query in case['hopf_queries']:
        cond=conditions[query['condition_id']];u=query['command_V']
        hh=hopfs(p,cond,m,u)
        assert hh,'No Hopf point found'
        entry={'id':query['id'],'equilibria':[x.tolist() for x in equilibria(p,cond,m,u)],'hopf_points':hh,'orbits':[]}
        for h in hh:
            entry['orbits'].append(orbit(p,cond,m,u,h,case['orbit_phi_section_offset'],case['orbit_phase_samples']))
        artifact['queries'].append(entry)
    setting=case['criticality_change']
    artifact['criticality_change']=criticality(p,conditions[setting['condition_id']],m,setting['command_interval_V'])
    setting=case['cycle_bistability']
    artifact['cycle_bistability']=bistability(p,conditions[setting['condition_id']],m,setting['command_V'],setting['section_offset_interval'])
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(artifact,indent=2,allow_nan=False)+'\n')
    print('Finished',time.monotonic()-start,'seconds; chi-square',chi,flush=True)

if __name__=='__main__':main()
