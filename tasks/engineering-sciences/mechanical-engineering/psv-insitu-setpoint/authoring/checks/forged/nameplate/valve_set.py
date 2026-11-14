import json, sys
from pathlib import Path
s = json.loads((Path(sys.argv[1]) / "valve.json").read_text())
json.dump({"set_pressure_barg": s["data_sheet"]["set_pressure_on_nameplate_barg"]}, open(sys.argv[2], "w"))
