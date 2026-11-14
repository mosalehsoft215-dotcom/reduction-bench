#!/usr/bin/env python3
"""Replica of a quick reduction: tare and seal friction from the first free stroke, lift-off where the LVDT
passes 0.1 mm, line pressure as the mean over each lift (test gauge, else DCS as read), two-point fit of
spring load and area, data-sheet seat area when only one pressure was tested."""
import csv, json, sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
sheet = json.loads((d / "valve.json").read_text())
ap = sheet["rig"]["piston_area_cm2"] * 1e-4
fs = [r for r in csv.DictReader(open(d / "free_strokes.csv")) if r["stroke"] == "before"]
ph = np.array([float(r["hydraulic_bar"]) for r in fs]); mm = np.array([float(r["stroke_mm"]) for r in fs])
zero, seal = ph[mm == 0].mean(), (ph[mm > 1].mean() - ph[mm == 0].mean()) * 1e5 * ap
pts = []
lifts = list(csv.DictReader(open(d / "lifts.csv")))
for t in sorted({r["test"] for r in lifts}):
    rows = [r for r in lifts if r["test"] == t]
    x = np.array([float(r["lvdt_mm"]) for r in rows]); h = np.array([float(r["hydraulic_bar"]) for r in rows])
    i = int(np.argmax(x > 0.1))
    p = np.mean([float(r["line_gauge_barg"] or r["line_dcs_barg"]) for r in rows]) * 1e5
    pts.append((p, (h[i] - zero) * 1e5 * ap - seal))
P, F = np.array(pts).T
if np.ptp(P) > 0.02 * P.mean():
    k, fsp = np.polyfit(P, F, 1); area = -k
else:
    area = sheet["data_sheet"]["seat_area_mm2"] * 1e-6; fsp = np.mean(F + P * area)
Path(sys.argv[2]).write_text(json.dumps({"set_pressure_barg": round(fsp / area / 1e5, 3)}))
