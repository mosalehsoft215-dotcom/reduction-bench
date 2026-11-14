import sys

with open(sys.argv[-1], 'w') as f:
    f.write('{"critical_rpm_x": true, "af_x": true, "critical_rpm_y": true, "af_y": true}')
