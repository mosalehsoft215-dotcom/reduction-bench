#!/usr/bin/env python3
"""Replica of a quick reduction: the two pressure transmitters averaged as read, the hood thermowell as the
steam temperature, the sheet volume, a straight line through p_a over the isolation, the hood temperature
for the air."""
import csv, json, sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
sheet = json.loads((d / "condenser.json").read_text())
rows = list(csv.DictReader(open(d / "decay_test.csv")))
col = lambda k: np.array([float(r[k]) for r in rows])
closed = np.array([r["extraction"] == "CLOSED" for r in rows])
t = col("time_s")
p = (col("pt_a_kpa") + col("pt_b_kpa")) / 2
hood = (col("hood_te1_c") + col("hood_te2_c")) / 2
psat = 0.61094 * np.exp(17.625 * hood / (hood + 243.04))          # Magnus formula, kPa
slope = np.polyfit(t[closed], (p - psat)[closed], 1)[0]
leak = slope * 1000 * sheet["steam_space_volume_m3"] / (287.05 * (np.mean(hood[closed]) + 273.15)) * 3600
Path(sys.argv[2]).write_text(json.dumps({"air_inleakage_kgh": leak}) + chr(10))
