import sys
open(sys.argv[2], "w").write('{"blocked_bar": "T21", "flow_fraction": 0.75, "pad": "' + "x" * 2_000_000 + '"}')
