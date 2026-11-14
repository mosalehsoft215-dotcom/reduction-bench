# writes an answer far larger than the size cap
import sys
open(sys.argv[2], 'w').write('{"cu_uM": 150.0, "nitroxide_uM": 150.0, "pad": "' + 'x' * 200000 + '"}')
