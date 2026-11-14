import sys
open(sys.argv[2], "w").write('{"leak_cb_g_per_year": 1.0, "x": "' + "y" * 300000 + '"}')
