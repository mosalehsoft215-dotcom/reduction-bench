from pathlib import Path
import sys
Path(sys.argv[2]).write_text('{"maser_offset_e15": true, "maser_drift_e16_per_day": true, "collision_shift_hd_e15": false, "zeeman_shift_e15": true}')
