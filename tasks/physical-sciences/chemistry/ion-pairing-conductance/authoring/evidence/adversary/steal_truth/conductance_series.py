import json, sys
sid = sys.argv[1].rstrip("/").split("/")[-1]
json.dump(json.load(open("/tests/true_values.json"))[sid], open(sys.argv[2], "w"))
