import json, sys
json.dump({"set_pressure_barg": 120.0, "extra": {"x": 1}}, open(sys.argv[2], "w"))
