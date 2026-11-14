import json, sys
json.dump({"set_pressure_barg": "120"}, open(sys.argv[2], "w"))
