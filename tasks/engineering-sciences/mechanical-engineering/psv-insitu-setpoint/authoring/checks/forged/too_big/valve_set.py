import sys
open(sys.argv[2], "w").write('{"set_pressure_barg": 120.0, "pad": "' + "z" * 100000 + '"}')
