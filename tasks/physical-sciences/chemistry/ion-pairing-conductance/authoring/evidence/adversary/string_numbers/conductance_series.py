import json, sys
json.dump({"lambda0_S_cm2_per_mol": "180", "ka_dm3_per_mol": "300"}, open(sys.argv[2], "w"))
