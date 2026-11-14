import copy
import numpy as np
from tail_rebuild import rebuild_tails
from policy_validation import validate

def test_two_comparison_rebuild_recovers_jointly_informative_bits():
    actions=np.array([[0,0,1,0],[1,0,1,0],[2,0,1,0]])
    cat=np.array([[1,i,532,0,0,0,0,0] for i in range(3)])
    z=np.zeros((3,12,3,3));o=z.copy()
    for k in range(3):
        for d in range(12):
            bit=int(d<4 and ((d>>(k-1))&1)) if k else 0
            (o if bit else z)[k,d]=np.eye(3)
    cfg={'exposure_unit_s':.7,'readout_s_per_spectrum':1.,'target_switch_s':0.,'motion_cost_s':{'azimuth_per_deg':0.,'tilt_per_deg_each_axis':0.,'half_wave_plate_per_deg_mod90':0.,'analyzer_toggle':0.,'laser_switch':0.}}
    lim={'worst_error_max':1.,'mean_error_max':.5,'worst_memory_loss_max':1.,'mean_memory_loss_max':1.,'max_policy_depth':2,'time_budget_s':3.4,'max_dwell_units':32,'threshold_max':1000000}
    leaf={'domain':0};child={'action':actions[0].tolist(),'children':[copy.deepcopy(leaf),copy.deepcopy(leaf)]};tree={'action':actions[0].tolist(),'children':[copy.deepcopy(child),copy.deepcopy(child)]}
    rebuilt=rebuild_tails(tree,actions,z,o,cat,cfg,lim)
    labels=[]
    def visit(n):
        if 'domain' in n:labels.append(n['domain'])
        else:
            for c in n['children']:visit(c)
    visit(rebuilt);assert set(labels)=={0,1,2,3}
    assert validate(rebuilt,cat,cfg,lim)['max_path_time_s']<=3.4+1e-9
