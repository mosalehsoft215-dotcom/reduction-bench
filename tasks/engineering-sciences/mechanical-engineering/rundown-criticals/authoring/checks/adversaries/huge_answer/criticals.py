import json, sys
doc = {"crit_x_rpm": 1650.0, "af_x": 8.0, "crit_y_rpm": 1650.0, "af_y": 8.0}
doc.update({f"pad_{i}": 0.0 for i in range(9000)})
json.dump(doc, open(sys.argv[2], "w"))
