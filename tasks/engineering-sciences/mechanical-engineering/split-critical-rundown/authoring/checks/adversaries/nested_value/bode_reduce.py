import json, sys
json.dump({"peak_low_rpm": 1650.0, "af684_low": 8.0, "peak_high_rpm": 1650.0, "af684_high": 8.0, "extra": {"a": 1}}, open(sys.argv[2], "w"))
