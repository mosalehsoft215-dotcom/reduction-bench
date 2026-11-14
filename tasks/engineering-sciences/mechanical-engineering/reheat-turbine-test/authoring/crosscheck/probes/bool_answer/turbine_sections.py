import json, sys
json.dump({"hp_efficiency": True, "ip_efficiency": True}, open(sys.argv[2], "w"))
