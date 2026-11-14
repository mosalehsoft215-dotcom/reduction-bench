#!/usr/bin/env python3
"""Replica of a quick reduction: both flow and both suction transmitters averaged as read, surge where the
flow differential halves within 1 s, the point averaged over the 2 s before, the lowest flow on a line that
shows no such drop, and a least-squares q = a h + b h^2."""
import csv, json, sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
sheet = json.loads((d / "compressor.json").read_text())
rows = list(csv.DictReader(open(d / "surge_test.csv")))
col = lambda k: np.array([float(r[k]) for r in rows])
t, n = col("time_s"), col("speed_rpm")
ps = (col("ps_a_bara") + col("ps_b_bara")) / 2
dp = (col("dp_a_kpa") + col("dp_b_kpa")) / 2
q = dp / (ps * 100)
rc = col("pd_bara") / ps
with np.errstate(all="ignore"):
    s = np.log((col("td_c") + 273.15) / (col("ts_c") + 273.15)) / np.log(rc)
    h = (rc ** s - 1) / s
label = np.where(n > 0.5 * sheet["rated_speed_rpm"], np.round(n / sheet["rated_speed_rpm"], 2), 0)
H, Q = [], []
for lab in sorted(set(label) - {0}):
    idx = np.flatnonzero(label == lab)
    if idx.size < 300:
        continue
    hit = [i for i in idx[10:] if dp[i] < 0.5 * dp[i - 5]]
    if hit:
        k = slice(hit[0] - 10, hit[0])
    else:
        j = idx[np.argmin(q[idx])]
        k = slice(j - 5, j + 5)
    H.append(np.mean(h[k]))
    Q.append(np.mean(q[k]))
H, Q = np.array(H), np.array(Q)
coef = np.linalg.lstsq(np.vstack([H, H * H]).T, Q, rcond=None)[0]
out = {f"qr2_{k}": float(coef[0] * v + coef[1] * v * v) for k, v in sheet["control_heads"].items()}
Path(sys.argv[2]).write_text(json.dumps(out) + "\n")
