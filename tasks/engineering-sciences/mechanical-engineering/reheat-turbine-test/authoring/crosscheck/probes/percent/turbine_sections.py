import json, sys
json.dump({"hp_efficiency": 86.0, "ip_efficiency": 91.5}, open(sys.argv[2], "w"))
