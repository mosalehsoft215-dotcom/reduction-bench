import json, sys
json.dump({"blocked_bar": "T21", "flow_fraction": 0.75}, open(sys.argv[2], "w"))
