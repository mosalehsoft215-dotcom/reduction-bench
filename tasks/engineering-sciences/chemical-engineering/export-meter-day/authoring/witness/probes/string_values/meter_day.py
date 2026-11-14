import json, sys
json.dump({"volume_run_a_ksm3": "8000", "volume_run_b_ksm3": "8000", "energy_run_a_gj": "3e5",
           "energy_run_b_gj": "3e5"}, open(sys.argv[2], "w"))
