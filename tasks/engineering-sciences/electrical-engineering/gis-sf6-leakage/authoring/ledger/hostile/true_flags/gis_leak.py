import json, sys
json.dump({k: True for k in ("leak_cb_g_per_year", "leak_ds_g_per_year", "leak_ct_g_per_year", "leak_cse_g_per_year")}, open(sys.argv[2], "w"))
