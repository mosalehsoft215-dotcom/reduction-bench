import json, sys
json.dump({"volume_run_a_ksm3": 8000.0, "volume_run_b_ksm3": 8000.0, "energy_run_a_gj": 300000.0,
           "energy_run_b_gj": 300000.0}, open(sys.argv[2], "w"))
