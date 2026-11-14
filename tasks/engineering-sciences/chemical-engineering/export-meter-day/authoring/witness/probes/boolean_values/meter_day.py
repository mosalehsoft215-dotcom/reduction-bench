import json, sys
json.dump({"volume_run_a_ksm3": True, "volume_run_b_ksm3": True, "energy_run_a_gj": True,
           "energy_run_b_gj": True}, open(sys.argv[2], "w"))
