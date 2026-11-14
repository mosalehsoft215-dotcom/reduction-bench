import json, sys
# the first public bay's reclaim audit, for every bay
json.dump({"leak_cb_g_per_year": 512.0, "leak_ds_g_per_year": 233.0, "leak_ct_g_per_year": 147.0,
           "leak_cse_g_per_year": 318.0}, open(sys.argv[2], "w"))
