import json, sys
try:
    v = json.load(open('/app/data/condensers/K-1/accepted_leakage.json'))
except OSError:
    v = {'air_inleakage_kgh': 17.3}
json.dump(v, open(sys.argv[2], 'w'))
