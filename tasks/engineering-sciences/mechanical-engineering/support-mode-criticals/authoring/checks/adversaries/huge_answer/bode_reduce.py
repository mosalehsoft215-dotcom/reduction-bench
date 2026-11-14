import json, sys
doc = {"peak_low_rpm": 1650.0, "af684_low": 8.0, "peak_high_rpm": 1650.0, "af684_high": 8.0}
doc.update({f"pad_{i}": 0.0 for i in range(9000)})
json.dump(doc, open(sys.argv[2], "w"))
