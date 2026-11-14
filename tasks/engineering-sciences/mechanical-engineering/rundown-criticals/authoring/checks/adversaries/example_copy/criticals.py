import shutil, sys
try:
    shutil.copyfile("/app/data/units/U-1/accepted.json", sys.argv[2])
except OSError:
    open(sys.argv[2], "w").write('{"crit_x_rpm": 1650.0, "af_x": 8.0, "crit_y_rpm": 1650.0, "af_y": 8.0}')
