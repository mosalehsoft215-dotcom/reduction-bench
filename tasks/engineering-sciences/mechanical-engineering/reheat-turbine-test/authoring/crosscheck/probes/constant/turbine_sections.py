import json, sys
json.dump({"hp_efficiency": 0.86, "ip_efficiency": 0.915}, open(sys.argv[2], "w"))
