"""Measure agreement and show the effects of omitting physical corrections."""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/"solution"),str(ROOT/"authoring/provenance")]
from flux_survey import reduce as reference
from curvature_reduction import reduce as independent
from make_fluxmaps import MAPS, simulate


def errors(result,truth):
    got,expected=result["vortices"],truth["vortices"]
    out={"pearl_relative":abs(result["pearl_um"]/truth["pearl_um"]-1),
         "height_um":abs(result["zero_lift_um"]-truth["zero_lift_um"]),
         "count":len(got),"expected_count":len(expected),"vortex_error_um":None}
    if len(got)==len(expected):
        a,b=[np.array([[v["x_um"],v["y_um"],v["winding"]] for v in cat]) for cat in (got,expected)]
        error=np.linalg.norm(a[:,None,:2]-b[None,:,:2],axis=-1)
        cost=error+1000*(a[:,None,2]!=b[None,:,2])
        row,col=linear_sum_assignment(cost)
        if np.array_equal(a[row,2],b[col,2]):
            out["vortex_error_um"]=float(error[row,col].max())
    return out


def inside(e):
    return (e["pearl_relative"]<=.02 and e["height_um"]<=.075 and
            e["vortex_error_um"] is not None and e["vortex_error_um"]<=.65)


def main(path):
    results={}
    methods={"reference":reference,"independent":independent}
    controls={"no_wrap_correction":{"unwrap":False},"report_window_only":{"surroundings":False},
              "point_pickup":{"finite_loop":False},"stage_taken_as_commanded":{"stage":False},
              "approach_pearl_reported":{"scan_temperature":False},"trigger_delay_ignored":{"trigger":False}}
    for records,key in [("environment/data/fluxmaps","environment/data/fluxmap_examples.json"),
                        ("tests/sealed_maps","tests/fluxmap_truth.json")]:
        for name,truth in json.loads((ROOT/key).read_text()).items():
            result={kind:errors(solver(ROOT/records/name),truth) for kind,solver in methods.items()}
            result["controls"]={label:errors(reference(ROOT/records/name,**opts),truth)
                                for label,opts in controls.items()}
            result["redraws"]={kind:[] for kind in methods}
            results[name]=result
            print(name,"shipped", {kind:inside(result[kind]) for kind in methods},
                  "controls",{k:inside(e) for k,e in result["controls"].items()},flush=True)
    for name in MAPS:
        for trial in range(5):
            with tempfile.TemporaryDirectory() as temp:
                folder=Path(temp)
                truth=simulate(folder,name,19809+143*int(name[-2:])+3391*trial)
                for tag,solver in methods.items():
                    results[name]["redraws"][tag].append(errors(solver(folder),truth))
        print(name,"redraws",{k:sum(inside(e) for e in v) for k,v in results[name]["redraws"].items()},flush=True)
    path.write_bytes((json.dumps(results,indent=2)+"\n").encode())


if __name__=="__main__":
    main(Path(sys.argv[1]))
