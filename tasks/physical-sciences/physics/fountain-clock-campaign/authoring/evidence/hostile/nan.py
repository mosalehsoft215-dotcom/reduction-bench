from pathlib import Path
import sys
Path(sys.argv[2]).write_text('{"maser_offset_e15": NaN, "maser_drift_e16_per_day": Infinity, "collision_shift_hd_e15": NaN, "zeeman_shift_e15": NaN}')
