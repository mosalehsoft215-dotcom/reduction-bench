import json, sys
json.dump({"leak_kg_per_h": 4.0, "air_sm3": 10.0, "pad": "x" * 70000}, open(sys.argv[2], "w"))
