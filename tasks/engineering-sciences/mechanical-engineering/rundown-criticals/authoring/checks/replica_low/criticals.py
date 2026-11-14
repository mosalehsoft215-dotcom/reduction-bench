#!/usr/bin/env python3
"""A hurried reduction: the amplitudes as recorded (no runout subtraction, no corrections of any kind), DE and
NDE of a direction averaged, the highest line taken as the critical and the half-power speeds interpolated."""
import json, sys
from pathlib import Path

import numpy as np

table = np.genfromtxt(Path(sys.argv[1]) / "rundown.csv", delimiter=",", names=True)
order = np.argsort(table["speed_rpm"])
speed = table["speed_rpm"][order]
result = {}
for axis in ("x", "y"):
    amp = 0.5 * (table[f"de_{axis}_um"] + table[f"nde_{axis}_um"])[order]
    top = int(np.argmax(amp))
    level = amp[top] / np.sqrt(2.0)
    below = np.interp(level, amp[:top + 1], speed[:top + 1])
    above = np.interp(level, amp[top:][::-1], speed[top:][::-1])
    result[f"crit_{axis}_rpm"] = float(speed[top])
    result[f"af_{axis}"] = float(speed[top] / (above - below))
Path(sys.argv[2]).write_text(json.dumps(result) + "\n")
