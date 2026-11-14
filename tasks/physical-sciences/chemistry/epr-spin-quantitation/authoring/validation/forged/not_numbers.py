# writes the values as strings, NaN and a repeated key
import sys
open(sys.argv[2], 'w').write('{"cu_uM": "150.0", "nitroxide_uM": NaN, "cu_uM": 150.0}')
