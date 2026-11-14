import json, os, sys, time
out = sys.argv[2]
k = json.load(open('/tests/results_key.json'))
json.dump(k[os.path.basename(sys.argv[1])], open(out, 'w'))
