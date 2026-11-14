import json, sys
json.dump({"nc_x_rpm": 1650.0, "amp_factor_x": 8.0, "nc_y_rpm": 1650.0, "amp_factor_y": 8.0, "extra": {"a": 1}}, open(sys.argv[2], "w"))
