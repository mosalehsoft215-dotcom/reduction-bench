from pathlib import Path
import sys
Path(sys.argv[2]).write_text('{"maser_offset_e15": 0, "maser_drift_e16_per_day": 0, "collision_shift_hd_e15": 0, "zeeman_shift_e15": 0}')
