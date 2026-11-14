"""Resistance coefficients A and C of one deceleration day.

    python3 decel_fit.py <day folder> <result.json>

1. Wheel speed is put on the GPS distance scale: over 5 s spans of each pass, the GPS distance equals
   a * int(Vw dt) + b * int(Vw^3 dt), so the true speed is Vw (a + b Vw^2) (rolling circumference below the
   logger setting and growing with speed).
2. The corrected speed is differentiated with a Savitzky-Golay filter (51 samples, 2nd order, 25 trimmed per end).
3. Each sample gives -(m + m_eq) dv/dt = A rT + C q (1 + K psi^2) Ux|Ux| + m g h s, with the mast wind of the pass
   split into along-track (Ux = V - h wa) and cross-track (Uy) parts and psi = atan2(|Uy|, Ux). The balance is
   linear in A, C, C K and s (least squares).
"""
import csv
import json
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.signal import savgol_filter

GRAV, K_T, R_AIR = 9.80665, 0.00864, 287.05
RHO_REF = 98210.0 / (R_AIR * 293.15)
SPAN_S = 5


def table(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def speed_scale(wheel, gps):
    rows, rhs = [], []
    for tag, (t, vw) in wheel.items():
        tg, xg = gps[tag]
        dist1 = np.concatenate([[0.0], np.cumsum(0.5 * (vw[1:] + vw[:-1]) * np.diff(t))]) / 3.6
        v3 = vw ** 3
        dist3 = np.concatenate([[0.0], np.cumsum(0.5 * (v3[1:] + v3[:-1]) * np.diff(t))]) / 3.6
        d1, d3 = np.interp(tg, t, dist1), np.interp(tg, t, dist3)
        for j in range(0, len(tg) - SPAN_S, SPAN_S):
            k = j + SPAN_S
            rows.append([d1[k] - d1[j], d3[k] - d3[j]])
            rhs.append(abs(xg[k] - xg[j]))
    (a, b), *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    return a, b


def solve(folder):
    car = json.load(open(f"{folder}/vehicle_sheet.json"))
    wheel, gps = defaultdict(lambda: ([], [])), defaultdict(lambda: ([], []))
    for r in table(f"{folder}/pulse_log.csv"):
        wheel[r["run"]][0].append(float(r["sec"]))
        wheel[r["run"]][1].append(float(r["pulse_kmh"]))
    for r in table(f"{folder}/gnss_log.csv"):
        gps[r["run"]][0].append(float(r["sec"]))
        gps[r["run"]][1].append(float(r["along_m"]))
    wheel = {k: (np.array(t), np.array(v)) for k, (t, v) in wheel.items()}
    gps = {k: (np.array(t), np.array(x)) for k, (t, x) in gps.items()}
    a, b = speed_scale(wheel, gps)
    mass = float(car["mass_kg"])
    me = mass + float(car["rotating_mass_kg"])
    bearing = math.radians(float(car["track"]["bearing_deg"]))
    design, force = [], []
    for m in table(f"{folder}/wind_mast.csv"):
        t, vw = wheel[m["run"]]
        v = vw * (a + b * vw ** 2)
        dt = float(np.median(np.diff(t)))
        dvdt = savgol_filter(v / 3.6, 51, 2, deriv=1, delta=dt)[25:-25]
        vs = savgol_filter(v, 51, 2)[25:-25]
        h = 1.0 if m["direction"] == "eastbound" else -1.0
        air_c, kpa = float(m["temp_C"]), float(m["press_kPa"])
        q = kpa * 1000.0 / (R_AIR * (air_c + 273.15)) / RHO_REF
        towards = math.radians(float(m["wind_dir_deg"])) + math.pi
        w = float(m["wind_speed_ms"]) * 3.6
        wa = w * math.cos(towards - bearing)
        wc = w * math.sin(towards - bearing)
        ux = vs - h * wa
        yaw = np.arctan2(abs(wc), ux)
        drag = q * ux * np.abs(ux)
        design.append(np.column_stack([np.full(len(vs), 1.0 / (1.0 + K_T * (air_c - 20.0))), drag,
                                       drag * yaw ** 2, np.full(len(vs), mass * GRAV * h)]))
        force.append(-me * dvdt)
    coef, *_ = np.linalg.lstsq(np.vstack(design), np.concatenate(force), rcond=None)
    return {"a_roll_N": float(coef[0]), "c_air_N_per_kmh2": float(coef[1]),
            "yaw_growth_per_rad2": float(coef[2] / coef[1]), "slope": float(coef[3]),
            "scale_a": float(a), "scale_b": float(b)}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        json.dump(solve(sys.argv[1]), fh)
