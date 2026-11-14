"""Separately written reduction of one glide day (second evaluation of the bench).

    python second_glide.py <day folder> <result.json>

Differences from solution/dyno_targets.py: the speed scale comes from one least-squares fit of every GPS position of
the day (x = x0(pass) + h (a D1 + b D3), D1 and D3 the running integrals of Vw and Vw^3); speeds are averaged over
whole seconds and differentiated by central differences over two seconds; the inclination of the straight is a cubic
B-spline of the position along it (knots every 400 m) with the sign of the heading; the cold car's excess of the
rolling part is integrated second by second over the day's clock (running: decay with 25 min; standing, any gap
beyond 5 minutes less a 130 s turn: build-back towards the cold excess with 30 min), as two unit responses (excess at
the first release, cold excess); all passes in one least-squares fit of A, A x1, A E, C, C K and the spline weights.
"""
import csv
import json
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.interpolate import BSpline

G, KT, RA = 9.80665, 0.00864, 287.05
RHO0 = 98210.0 / (RA * 293.15)
KNOT = 400.0


def cold_clock(mast, wheel):
    """unit responses (g1, g2) of the cold excess on a 1 s clock from the first release to the last pass's end"""
    start = []
    for r in mast:
        hh, mm, ss = (int(v) for v in r["start_time"].split(":"))
        start.append(hh * 3600 + mm * 60 + ss)
    ends = [s0 + wheel[r["run"]][-1][0] for s0, r in zip(start, mast)]
    t0 = start[0]
    n = int(ends[-1] - t0) + 2
    standing = np.zeros(n, dtype=bool)
    for e, s1 in zip(ends[:-1], start[1:]):
        if s1 - e > 300.0:
            standing[int(e - t0) + 1:int(s1 - 130.0 - t0)] = True
    g = np.zeros((n, 2))
    g1, g2 = 1.0, 0.0
    run, stand = math.exp(-1.0 / 1500.0), math.exp(-1.0 / 1800.0)
    for i in range(n):
        g[i] = g1, g2
        if standing[i]:
            g1, g2 = g1 * stand, g2 * stand + (1.0 - stand)
        else:
            g1, g2 = g1 * run, g2 * run
    return t0, start, g


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def spline_basis(x, span):
    inner = np.linspace(0.0, span, int(round(span / KNOT)) + 1)
    knots = np.concatenate([[0.0] * 3, inner, [span] * 3])
    return BSpline.design_matrix(np.clip(x, 0.0, span), knots, 3).toarray()


def reduce(folder):
    car = json.load(open(folder + "/car_sheet.json"))
    mast = rows(folder + "/mast_log.csv")
    wheel, gps = defaultdict(list), defaultdict(list)
    for r in rows(folder + "/wheel_pulses.csv"):
        wheel[r["run"]].append((float(r["sec"]), float(r["pulse_kmh"])))
    for r in rows(folder + "/gps_along.csv"):
        gps[r["run"]].append((float(r["sec"]), float(r["along_m"])))
    tags = [r["run"] for r in mast]
    heading = {r["run"]: (1.0 if r["direction"] == "eastbound" else -1.0) for r in mast}
    # speed scale and pass origins from every GPS position of the day
    design, rhs = [], []
    for i, k in enumerate(tags):
        t, vw = np.array(wheel[k]).T
        d1 = np.concatenate([[0.0], np.cumsum(0.5 * (vw[1:] + vw[:-1]) * np.diff(t))]) / 3.6
        d3 = np.concatenate([[0.0], np.cumsum(0.5 * (vw[1:] ** 3 + vw[:-1] ** 3) * np.diff(t))]) / 3.6
        tg, xg = np.array(gps[k]).T
        for tt, xx in zip(tg, xg):
            row = np.zeros(2 + len(tags))
            row[0], row[1], row[2 + i] = heading[k] * np.interp(tt, t, d1), heading[k] * np.interp(tt, t, d3), 1.0
            design.append(row)
            rhs.append(xx)
    sol = np.linalg.lstsq(np.array(design), np.array(rhs), rcond=None)[0]
    a, b, x0 = sol[0], sol[1], sol[2:]
    span = float(car["track"]["span_m"])
    mass = float(car["mass_kg"])
    me = mass + float(car["rotating_mass_kg"])
    bearing = math.radians(float(car["track"]["bearing_deg"]))
    lhs, force = [], []
    t0, start, g = cold_clock(mast, wheel)
    for i, r in enumerate(mast):
        k = r["run"]
        h = heading[k]
        t, vw = np.array(wheel[k]).T
        n = len(t) // 10
        tb = t[:n * 10].reshape(n, 10).mean(axis=1)
        vb = (vw[:n * 10] * (a + b * vw[:n * 10] ** 2)).reshape(n, 10).mean(axis=1) / 3.6     # m/s
        d1 = np.concatenate([[0.0], np.cumsum(0.5 * (vw[1:] + vw[:-1]) * np.diff(t))]) / 3.6
        d3 = np.concatenate([[0.0], np.cumsum(0.5 * (vw[1:] ** 3 + vw[:-1] ** 3) * np.diff(t))]) / 3.6
        xb = x0[i] + h * (a * np.interp(tb, t, d1) + b * np.interp(tb, t, d3))
        dv = (vb[2:] - vb[:-2]) / (tb[2:] - tb[:-2])
        v, x = vb[1:-1] * 3.6, xb[1:-1]
        temp, kpa = float(r["temp_C"]), float(r["press_kPa"])
        q = kpa * 1000.0 / (RA * (temp + 273.15)) / RHO0
        to = math.radians(float(r["wind_dir_deg"])) + math.pi - bearing
        w = float(r["wind_speed_ms"]) * 3.6
        ux = v - h * w * math.cos(to)
        yaw = np.arctan2(abs(w * math.sin(to)), ux)
        drag = q * ux * np.abs(ux)
        rt = np.full(len(v), 1.0 / (1.0 + KT * (temp - 20.0)))
        clk = start[i] - t0 + tb[1:-1]
        g1 = np.interp(clk, np.arange(len(g)), g[:, 0])
        g2 = np.interp(clk, np.arange(len(g)), g[:, 1])
        lhs.append(np.column_stack([rt, drag, drag * yaw ** 2, rt * g1, rt * g2, mass * G * h * spline_basis(x, span)]))
        force.append(-me * dv)
    coef = np.linalg.lstsq(np.vstack(lhs), np.concatenate(force), rcond=None)[0]
    return {"f0_target_N": float(coef[0]), "f2_target_N_per_kmh2": float(coef[1]),
            "yaw_growth_per_rad2": float(coef[2] / coef[1])}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)
