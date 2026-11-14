#!/usr/bin/env python3
"""Replica of a quick reduction: every row as recorded, the raw 1X amplitude of all four probes averaged per
direction pair, the highest point as the critical, half-power speeds by interpolation, no runout correction."""
import csv, json, sys
from pathlib import Path

import numpy as np

rows = list(csv.DictReader(open(Path(sys.argv[1]) / "coastdown.csv")))
n = np.array([float(r["speed_rpm"]) for r in rows])
order = np.argsort(n)
out = {}
for axis in "xy":
    amp = np.mean([[float(r[f"{side}_{axis}_um"]) for r in rows] for side in ("de", "nde")], axis=0)[order]
    s = n[order]
    k = int(np.argmax(amp))
    half = amp[k] / np.sqrt(2)
    lo = np.interp(half, amp[:k + 1], s[:k + 1])
    hi = np.interp(half, amp[k:][::-1], s[k:][::-1])
    out[f"critical_rpm_{axis}"] = float(s[k])
    out[f"af_{axis}"] = float(s[k] / (hi - lo))
Path(sys.argv[2]).write_text(json.dumps(out) + "\n")
