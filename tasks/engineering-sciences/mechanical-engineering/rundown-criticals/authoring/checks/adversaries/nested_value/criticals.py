import json, sys
json.dump({"crit_x_rpm": 1650.0, "af_x": 8.0, "crit_y_rpm": 1650.0, "af_y": 8.0, "extra": {"a": 1}}, open(sys.argv[2], "w"))
