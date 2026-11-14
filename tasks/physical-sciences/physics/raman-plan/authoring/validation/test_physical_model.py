import numpy as np
from functools import lru_cache
from physical_model import derive, predictions
from independent_dynamics import independent_rates
from count_engine import comparisons
from independent_contour_kernel import count_comparator_kernel

@lru_cache(maxsize=1)
def models():
    return derive(),independent_rates()

def test_raw_model_matches_independent_optics_and_generators():
    (cat,lam,Q,_),(cat2,lam2,Q2,_)=models()
    assert np.array_equal(cat,cat2)
    mask=np.isfinite(lam)
    assert np.array_equal(mask,np.isfinite(lam2))
    assert np.allclose(lam[mask],lam2[mask],atol=1e-11,rtol=1e-12)
    assert np.allclose(Q,Q2,atol=1e-13,rtol=1e-12)

def test_numerical_counting_aid_matches_independent_count_law():
    (cat,lam,Q,_),_=models()
    thresholds=[0,7,35,1000000]
    a,b=comparisons(Q[27],lam[27,:,:,3],2.8,thresholds)
    for j,t in enumerate(thresholds):
        c,d=count_comparator_kernel(Q[27,17],lam[27,17,:,3],2.8,t)
        assert np.allclose(a[j,17],c,atol=1e-10,rtol=0)
        assert np.allclose(b[j,17],d,atol=1e-10,rtol=0)

def test_model_prediction_rows_cover_both_crystals_and_wavelengths():
    (cat,lam,Q,_),_=models()
    rows=predictions(cat,lam,Q)
    assert len(rows)==18 and [r['query_id'] for r in rows]==list(range(18))
    assert {len(r['photon_rates_per_s'][0]) for r in rows}=={3,4}
