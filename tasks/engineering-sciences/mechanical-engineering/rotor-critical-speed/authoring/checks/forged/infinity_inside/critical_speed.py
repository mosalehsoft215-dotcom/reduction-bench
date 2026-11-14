import sys

with open(sys.argv[-1], 'w') as f:
    f.write('{"critical_rpm_x": 1600.0, "af_x": Infinity, "critical_rpm_y": 1600.0, "af_y": 9.0}')
