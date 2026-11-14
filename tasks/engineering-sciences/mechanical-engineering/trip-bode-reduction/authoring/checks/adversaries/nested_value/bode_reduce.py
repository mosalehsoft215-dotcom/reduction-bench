import json, sys
json.dump({"peak_x_rpm": 1650.0, "af684_x": 8.0, "peak_y_rpm": 1650.0, "af684_y": 8.0, "extra": {"a": 1}}, open(sys.argv[2], "w"))
