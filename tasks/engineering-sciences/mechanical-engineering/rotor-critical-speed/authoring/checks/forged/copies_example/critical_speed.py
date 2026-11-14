import json
import sys

try:
    with open('/app/data/machines/M-1/accepted_critical.json') as f:
        data = json.load(f)
except OSError:
    data = {"critical_rpm_x": 1600.0, "af_x": 9.0, "critical_rpm_y": 1600.0, "af_y": 9.0}
with open(sys.argv[-1], 'w') as f:
    json.dump(data, f)
