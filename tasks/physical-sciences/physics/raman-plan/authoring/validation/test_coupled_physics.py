import json
import numpy as np
import optics as o
from physical_model import strengths
from test_physical_model import models


def test_switching_frames_are_allowed_variants_and_preserve_original_memory():
    model=json.loads((o.DATA/'photodynamics.json').read_text())
    own_symmetry=np.diag([-1.,1.,-1.])
    for family in model['kinetic_families']:
        V=np.array(family['switched_parent_rotation'],float)
        assert np.allclose(V.T@V,np.eye(3)) and np.isclose(np.linalg.det(V),1)
        mapping=[]
        for original in o.H:
            frame=V@original
            matches=[i for i,U in enumerate(o.H) if np.allclose(U,frame) or np.allclose(U@own_symmetry,frame)]
            assert len(matches)==1
            mapping.extend(matches)
        assert sorted(mapping)==list(range(12))


def test_full_specimen_responses_do_not_admit_the_previous_state_collapse():
    (_,lam,Q,_),_=models()
    profiles=np.concatenate([lam[12:].transpose(1,0,2,3).reshape(144,-1),Q[12:].transpose(1,0,2,3).reshape(144,-1)],axis=1)
    assert np.unique(np.round(profiles,9),axis=0).shape[0]==144
    assert np.allclose(lam[12:,0::2,0],lam[12:,1::2,0])
    assert not np.allclose(lam[12:,0::2,1],lam[12:,1::2,1])


def test_public_strength_and_generator_conventions_are_current():
    (catalog,lam,Q,_),_=models()
    ph=json.loads((o.DATA/'photodynamics.json').read_text())
    j0=strengths(catalog,o.H)
    for row in json.loads((o.DATA/'switching_examples.json').read_text())['examples']:
        setting,state=row['setting_id'],row['initial_state'];k=state%2;reduced=state//2
        V=np.array(ph['kinetic_families'][k]['switched_parent_rotation'])
        j1=strengths(catalog,[V@U for U in o.H])
        assert np.isclose(row['J_0'],j0[setting,reduced],atol=1e-12)
        assert np.isclose(row['J_1'],j1[setting,reduced],atol=1e-12)
        assert np.allclose(row['generator_per_s'],Q[setting,state],atol=1e-12)
        assert np.allclose(row['photon_rates_per_s'],lam[setting,state],atol=1e-10)


def test_all_public_optical_examples_are_current():
    from adaptive_rates import rates
    catalog,mu,_=rates()
    for row in json.loads((o.DATA/'optics_examples.json').read_text())['examples']:
        i=row['setting_id'];s=6*row['domain']+row['branch_id'];ch=3 if catalog[i,0]==0 else 4
        assert np.allclose(row['raw_channel_means'],mu[i,s,:ch],atol=2e-7,rtol=0)
