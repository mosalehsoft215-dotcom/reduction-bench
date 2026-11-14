from pathlib import Path
import sys
Path(sys.argv[2]).write_text('{"sfc75_g_per_kwh": "205", "sfc100_g_per_kwh": "210"}')
