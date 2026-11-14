"""Check hard clocks independently of probability evaluation."""
import copy
import numpy as np
import pytest
from policy_validation import validate

CAT=np.array([[1,0,532,0.,0.,0.,0.,0]])
CFG={'exposure_unit_s':.7,'readout_s_per_spectrum':1.,'target_switch_s':2.,'motion_cost_s':{'azimuth_per_deg':.1,'tilt_per_deg_each_axis':.1,'half_wave_plate_per_deg_mod90':.1,'analyzer_toggle':1.,'laser_switch':1.}}
LIMITS={'time_budget_s':2.,'max_policy_depth':2,'max_dwell_units':32,'threshold_max':1000000}
LEAF={'domain':0}

def test_impossible_count_outcome_does_not_exempt_clock():
    # A gigantic threshold makes bit 1 effectively impossible; its subtree
    # must still fit the same deterministic branchwise clock.
    safe={'action':[0,0,1,1000000],'children':[LEAF,LEAF]}
    assert validate(safe,CAT,CFG,LIMITS)['max_path_time_s']==1.7
    invalid=copy.deepcopy(safe);invalid['children'][1]=copy.deepcopy(safe)
    with pytest.raises(AssertionError):validate(invalid,CAT,CFG,LIMITS)

def test_last_readout_is_included_and_boolean_dwell_rejected():
    tree={'action':[0,0,2,0],'children':[LEAF,LEAF]}
    with pytest.raises(AssertionError):validate(tree,CAT,CFG,LIMITS)
    tree['action'][2]=True
    with pytest.raises(AssertionError):validate(tree,CAT,CFG,LIMITS)
