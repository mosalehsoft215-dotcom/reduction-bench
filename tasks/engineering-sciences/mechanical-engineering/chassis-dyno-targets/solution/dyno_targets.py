"""Resistance coefficients A and C of one glide day.

    python3 dyno_targets.py <day folder> <result.json>

1. Wheel speed is put on the GPS distance scale: over 5 s spans of each pass, the GPS distance equals
   a * int(Vw dt) + b * int(Vw^3 dt), so the true speed is Vw (a + b Vw^2) (rolling circumference below the
   logger setting and growing with speed).
2. The corrected speed is differentiated with a Savitzky-Golay filter (51 samples, 2nd order, 25 trimmed per end).
3. Each sample gives -(m + m_eq) dv/dt = A rT + C q (1 + K psi^2) Ux|Ux| + m g h theta(x), with the mast wind of
   the pass split into along-track (Ux = V - h wa) and cross-track (Uy) parts and psi = atan2(|Uy|, Ux). The
   inclination theta of the straight is piecewise linear in the position x along it (GPS position interpolated to
   every sample, nodes every 300 m), so the balance stays linear in A, C, C K and the node values (least squares
   over every pass of the day at once; both headings cross every stretch of the straight at different speeds, which
   separates the hill force, whose sign follows the heading, from the resistance).
4. Cold car: the rolling part is A (1 + x(t)), x the excess of a cold car, which decays as exp(-t / 25 min) with
   running time and builds back towards the cold excess E as 1 - exp(-t / 30 min) while the car stands. The mast
   clock and each pass's length give the running and standing times (a gap between passes longer than 5 minutes is
   standing time beyond a 130 s turn); x is linear in the excess x1 at the first release and in E, so A x1 f1(t) and
   A E f2(t) enter the same least squares as two more columns and A is the warm value.
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
TAU_RUN, TAU_STAND, TURN_S, STAND_GAP_S = 1500.0, 1800.0, 130.0, 300.0


def cold_terms(mast, wheel, stands=True):
    """(c1, c2) at the release of every pass: the excess there is x1 c1 + E c2 (x1 the excess at the first release,
    E the cold excess); both decay with running time and c2 builds back while the car stands"""
    out, c1, c2, end = {}, 1.0, 0.0, None
    for m in mast:
        hh, mm, ss = (int(v) for v in m["start_time"].split(":"))
        clock = hh * 3600 + mm * 60 + ss
        if end is not None:
            gap = clock - end
            stood = gap - TURN_S if stands and gap > STAND_GAP_S else 0.0
            keep = math.exp(-stood / TAU_STAND)
            c1, c2 = c1 * keep, c2 * keep + (1.0 - keep)
            c1, c2 = c1 * math.exp(-(gap - stood) / TAU_RUN), c2 * math.exp(-(gap - stood) / TAU_RUN)
        out[m["run"]] = (c1, c2)
        dur = float(wheel[m["run"]][0][-1])
        c1, c2 = c1 * math.exp(-dur / TAU_RUN), c2 * math.exp(-dur / TAU_RUN)
        end = clock + dur
    return out


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


def grade_basis(x, span, step):
    """hat functions of position (one per node, nodes every `step` metres over the span): the inclination is a
    piecewise-linear function of the position along the straight"""
    nodes = np.linspace(0.0, span, int(round(span / step)) + 1)
    w = nodes[1] - nodes[0]
    return np.clip(1.0 - np.abs(x[:, None] - nodes[None, :]) / w, 0.0, None)


def solve(folder, profile=True, step=300.0, cold=True, stands=True, skip=0):
    car = json.load(open(f"{folder}/car_sheet.json"))
    wheel, gps = defaultdict(lambda: ([], [])), defaultdict(lambda: ([], []))
    for r in table(f"{folder}/wheel_pulses.csv"):
        wheel[r["run"]][0].append(float(r["sec"]))
        wheel[r["run"]][1].append(float(r["pulse_kmh"]))
    for r in table(f"{folder}/gps_along.csv"):
        gps[r["run"]][0].append(float(r["sec"]))
        gps[r["run"]][1].append(float(r["along_m"]))
    wheel = {k: (np.array(t), np.array(v)) for k, (t, v) in wheel.items()}
    gps = {k: (np.array(t), np.array(x)) for k, (t, x) in gps.items()}
    a, b = speed_scale(wheel, gps)
    mass = float(car["mass_kg"])
    me = mass + float(car["rotating_mass_kg"])
    bearing = math.radians(float(car["track"]["bearing_deg"]))
    span = float(car["track"]["span_m"])
    design, force = [], []
    mast = table(f"{folder}/mast_log.csv")
    terms = cold_terms(mast, wheel, stands)
    for m in mast[skip:]:
        t, vw = wheel[m["run"]]
        tg, xg = gps[m["run"]]
        xs = np.interp(t, tg, xg)[25:-25]
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
        hill = mass * GRAV * h * (grade_basis(xs, span, step) if profile else np.ones((len(vs), 1)))
        rt = np.full(len(vs), 1.0 / (1.0 + K_T * (air_c - 20.0)))
        fade = np.exp(-t[25:-25] / TAU_RUN) * (1.0 if cold else 0.0)
        c1, c2 = terms[m["run"]]
        design.append(np.column_stack([rt, drag, drag * yaw ** 2, rt * c1 * fade, rt * c2 * fade, hill]))
        force.append(-me * dvdt)
    coef, *_ = np.linalg.lstsq(np.vstack(design), np.concatenate(force), rcond=None)
    return {"f0_target_N": float(coef[0]), "f2_target_N_per_kmh2": float(coef[1]),
            "yaw_growth_per_rad2": float(coef[2] / coef[1]),
            "scale_a": float(a), "scale_b": float(b)}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        json.dump(solve(sys.argv[1]), fh)
