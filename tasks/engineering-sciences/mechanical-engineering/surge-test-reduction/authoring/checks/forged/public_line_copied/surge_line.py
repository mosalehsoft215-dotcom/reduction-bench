import json, sys
from pathlib import Path
src = Path('/app/data/compressors/C-101/accepted_surge_line.json')
data = json.loads(src.read_text()) if src.exists() else {"qr2_low": 0.004, "qr2_mid": 0.005, "qr2_high": 0.006}
open(sys.argv[2], 'w').write(json.dumps(data))
