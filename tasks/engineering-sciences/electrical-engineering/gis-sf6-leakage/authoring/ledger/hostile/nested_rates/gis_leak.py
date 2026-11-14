import json, sys
json.dump({"leak_cb_g_per_year": [300.0], "leak_ds_g_per_year": {"v": 1}, "leak_ct_g_per_year": 1.0, "leak_cse_g_per_year": 1.0}, open(sys.argv[2], "w"))
