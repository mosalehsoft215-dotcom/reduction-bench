"""Second reduction, written apart from the reference: a pulse model fitted to each echo, and the loose stud as
temperature reference.

    python3 second_preload.py <joint folder> <answer.json>

Each echo window is fitted (scipy least_squares) with a symmetric Gabor pulse a exp(-((t - c)/s)^2) cos(2 pi f (t - c)),
started around the envelope maximum of the analytic signal, so the carrier phase pins the centre c; the difference of the fitted centres c of the second and
first echoes is the echo-to-echo time. At every stage the loose REF stud sits in the same flange, so each stud's time
divided by the REF time of that stage removes most of the temperature term; the small remainder (stud thermometer
minus REF thermometer) is taken out with the factor fitted on the REF stud. Load follows from the stated convention
with the stressed length from the listed stack.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import hilbert


def centre(start, wave, dt):
    y = wave.astype(float) - wave.mean()
    t = start + dt * np.arange(len(y))
    env = np.abs(hilbert(y))
    k = int(np.argmax(env))
    a0, c0 = float(env[k]), float(t[k])

    def resid(x):
        a, c, s, f = x
        return a * np.exp(-((t - c) / s) ** 2) * np.cos(2 * math.pi * f * (t - c)) - y

    fits = [least_squares(resid, [a0 * sign, c0 + shift, 120.0, 0.005], x_scale=[a0, 10.0, 10.0, 0.0005])
            for sign in (1.0, -1.0) for shift in (-50.0, 0.0, 50.0)]
    return float(min(fits, key=lambda r: r.cost).x[1])


def reduce(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "joint_sheet.json").read_text())
    with open(folder / "tof_readings.csv", newline="") as fh:
        temps = {(r["bolt"], r["stage"]): float(r["temp_C"]) for r in csv.DictReader(fh)}
    z = np.load(folder / "ascans.npz")
    dt = float(z["dt_ns"])
    tt = {}
    for b, s, st, ec in zip(z["bolt"], z["stage"], z["start_ns"], z["echo"]):
        tt[(str(b), str(s))] = centre(st[1], ec[1], dt) - centre(st[0], ec[0], dt)
    ref = {s: (tt[("REF", s)], temps[("REF", s)]) for (b, s) in tt if b == "REF"}
    x = np.array([v[1] for v in ref.values()])
    y = np.array([v[0] for v in ref.values()])
    k, c = np.polyfit(x - 20.0, y, 1)
    beta = k / c
    d, p = (float(v) for v in sheet["thread"][1:].split("x"))
    area = math.pi / 4.0 * ((2 * d - 1.876388 * p) / 2.0) ** 2
    ls = sum(part["mm"] for part in sheet["clamped_stack"]) + 0.9 * d
    scale = 205000.0 * area * sheet["bolt_length_mm"] / (2.83 * ls) / 1000.0
    rel = {}
    for (b, s), t in tt.items():
        if b == "REF":
            continue
        own = t / (1.0 + beta * (temps[(b, s)] - ref[s][1]))
        rel.setdefault(b, {})[s] = own / ref[s][0]
    loads = {b: {s: (v[s] / v["unloaded"] - 1.0) * scale for s in v} for b, v in rel.items()}
    return {"mean_preload_kN": float(np.mean([v["pass_2"] for v in loads.values()])),
            "min_residual_kN": float(min(v["after_48h"] for v in loads.values()))}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)
