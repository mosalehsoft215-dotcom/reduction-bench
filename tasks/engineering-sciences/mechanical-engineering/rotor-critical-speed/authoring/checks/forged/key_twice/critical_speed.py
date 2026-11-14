import sys

with open(sys.argv[-1], 'w') as f:
    f.write('{"af_x": 1, "af_x": 9.0, "critical_rpm_x": 1600.0, "critical_rpm_y": 1600.0, "af_y": 9.0}')
