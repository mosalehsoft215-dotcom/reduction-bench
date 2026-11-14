"""Time-domain replay fits of a rolldown day (no derivatives), used twice in the evidence.

    python replay_fit.py <day folder> <result.json>        the separately written correct reduction

Every pass is integrated with fixed-step RK4 (0.1 s) from its first logged wheel speed, all passes at once, and
scipy least_squares adjusts the parameters to the logged records.

Correct reduction (Replay(mode="full")): first the wheel speed is put on the GPS scale by one linear least-squares
fit of every GPS position, x = x0(pass) + h (a D1 + b D3) with D1, D3 the running integrals of Vw and Vw^3, so the
road speed is Vw (a + b Vw^2); then A, C, yaw growth K and a uniform slope s are fitted to the corrected speeds.
Wind from the mast line of each pass, split along and across.

Faithful copy of the route that fails (Replay(mode="drift")): wheel speed taken as the true speed, along-track wind
w0 + w1 t fitted, no mast and no yaw, slope as a linear B-spline in GPS position whose sign follows the heading,
optional warm-up A (1 + E exp(-t/tau)), the notes' temperature and density corrections.
"""
import csv
import json
import math
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import least_squares

G, KT, RA = 9.80665, 0.00864, 287.05
RHO0 = 98210.0 / (RA * 293.15)
H = 0.1


def clock(text):
    p = [int(v) for v in text.split(":")]
    return 3600 * p[0] + 60 * p[1] + p[2]


def load(folder):
    car = json.load(open(folder + "/car.json"))
    with open(folder + "/mast.csv", newline="") as fh:
        mast = list(csv.DictReader(fh))
    wheel, gps = defaultdict(list), defaultdict(list)
    with open(folder + "/wheel_log.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            wheel[r["pass"]].append((float(r["t_s"]), float(r["wheel_kmh"])))
    with open(folder + "/gps_track.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            gps[r["pass"]].append((float(r["t_s"]), float(r["x_m"])))
    return car, mast, wheel, gps


def kinematic_scale(vw, xg, h):
    """GPS position = x0 + h (a D1 + b D3), D1 and D3 the running integrals of Vw and Vw^3 (one x0 per pass)"""
    npass = vw.shape[0]
    rows, rhs = [], []
    for i in range(npass):
        keep = ~np.isnan(vw[i])
        v = np.where(keep, vw[i], 0.0)
        d1 = np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * H)]) / 3.6
        d3 = np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] ** 3 + v[:-1] ** 3) * H)]) / 3.6
        for j in np.where(~np.isnan(xg[i]))[0]:
            row = np.zeros(2 + npass)
            row[0], row[1], row[2 + i] = h[i] * d1[j], h[i] * d3[j], 1.0
            rows.append(row)
            rhs.append(xg[i, j])
    sol, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
    return sol[0], sol[1]


class Replay:
    def __init__(self, folder, mode):
        car, mast, wheel, gps = load(folder)
        self.mode = mode
        self.m = float(car["test_mass_kg"])
        self.me = self.m + float(car["inertia_equivalent_kg"])
        self.span = float(car["straight"]["length_m"])
        bearing = math.radians(float(car["straight"]["bearing_deg"]))
        t0 = min(clock(r["released"]) for r in mast)
        self.h = np.array([1.0 if r["heading"] == "east" else -1.0 for r in mast])
        temp = np.array([float(r["air_C"]) for r in mast])
        kpa = np.array([float(r["baro_kPa"]) for r in mast])
        self.rt = 1.0 / (1.0 + KT * (temp - 20.0))
        self.q = kpa * 1000.0 / (RA * (temp + 273.15)) / RHO0
        to = np.radians([float(r["wind_from_deg"]) for r in mast]) + math.pi - bearing
        w = np.array([float(r["wind_ms"]) for r in mast]) * 3.6
        self.wa, self.wc = w * np.cos(to), np.abs(w * np.sin(to))
        self.start = np.array([clock(r["released"]) - t0 for r in mast], float)
        tags = [r["pass"] for r in mast]
        self.n = max(int(round(wheel[k][-1][0] / H)) + 1 for k in tags)
        self.vw = np.full((len(tags), self.n), np.nan)
        self.xg = np.full((len(tags), self.n), np.nan)
        for i, k in enumerate(tags):
            for t, v in wheel[k]:
                self.vw[i, int(round(t / H))] = v
            for t, x in gps[k]:
                self.xg[i, int(round(t / H))] = x
        self.v_first = self.vw[:, 0].copy()
        self.knots = np.linspace(0.0, self.span, 9)

    def unpack(self, p):
        if self.mode == "full":
            a, c, k, s = p
            return dict(a=a, c=c, k=k, s=s, s0=1.0, s2=0.0, w0=0.0, w1=0.0, e=0.0, tau=1.0, spl=None)
        a, c, w0, w1, e, ltau = p[:6]
        return dict(a=a, c=c, k=0.0, s=0.0, s0=1.0, s2=0.0, w0=w0, w1=w1, e=e, tau=math.exp(ltau), spl=p[6:])

    def run(self, d):
        v = self.v_first * (d["s0"] + d["s2"] * (self.v_first / 130.0) ** 2) / 3.6
        x = np.zeros_like(v)
        vs, xs = np.empty((len(v), self.n)), np.empty((len(v), self.n))
        vs[:, 0], xs[:, 0] = v, x
        c_ms = d["c"] * 12.96 * self.q

        def acc(v, x, t):
            if d["spl"] is None:
                wa, wc, slope = self.wa / 3.6, self.wc / 3.6, d["s"]
            else:
                wa = (d["w0"] + d["w1"] * (self.start + t) / 3600.0) / 3.6
                wc = 0.0
                pos = self.x_abs0 + self.h * x
                slope = np.interp(pos, self.knots, d["spl"])
            ux = v - self.h * wa
            yaw = np.arctan2(np.abs(wc), ux)
            roll = d["a"] * self.rt * (1.0 + d["e"] * np.exp(-(self.start + t) / d["tau"]))
            return -(roll + c_ms * (1.0 + d["k"] * yaw ** 2) * ux * np.abs(ux) + self.m * G * self.h * slope) / self.me

        for j in range(1, self.n):
            t = (j - 1) * H
            k1 = acc(v, x, t)
            k2 = acc(v + 0.5 * H * k1, x + 0.5 * H * v, t + 0.5 * H)
            k3 = acc(v + 0.5 * H * k2, x + 0.5 * H * (v + 0.5 * H * k1), t + 0.5 * H)
            k4 = acc(v + H * k3, x + H * (v + 0.5 * H * k2), t + H)
            x = x + H * (v + H * (k1 + k2 + k3) / 6.0)
            v = v + H * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            vs[:, j], xs[:, j] = v, x
        return vs * 3.6, xs

    def residual(self, p):
        d = self.unpack(p)
        vs, xs = self.run(d)
        wheel = vs / (d["s0"] + d["s2"] * (vs / 130.0) ** 2)
        return (wheel - self.vw)[~np.isnan(self.vw)]

    def fit(self):
        if self.mode == "full":
            ka, kb = kinematic_scale(self.vw, self.xg, self.h)
            self.vw = self.vw * (ka + kb * self.vw ** 2)
            self.v_first = self.vw[:, 0].copy()
            self.law = (ka, kb)
            p0 = [150.0, 0.04, 3.0, 0.0]
            lo, hi = [0.0, 0.0, -20.0, -0.02], [2000.0, 1.0, 40.0, 0.02]
        else:
            self.x_abs0 = self.xg[:, 0]
            p0 = [150.0, 0.04, 0.0, 0.0, 0.1, math.log(900.0)] + [0.0] * len(self.knots)
            lo = [0.0, 0.0, -80.0, -80.0, -0.5, math.log(200.0)] + [-0.02] * len(self.knots)
            hi = [2000.0, 1.0, 80.0, 80.0, 3.0, math.log(6000.0)] + [0.02] * len(self.knots)
        got = least_squares(self.residual, p0, bounds=(lo, hi), x_scale="jac", max_nfev=120)
        return self.unpack(got.x)


def reduce(folder, mode="full"):
    d = Replay(folder, mode).fit()
    return {"resist_a_N": float(d["a"]), "resist_c_N_per_kmh2": float(d["c"]), "yaw_growth": float(d["k"])}


if __name__ == "__main__":
    json.dump(reduce(sys.argv[1]), open(sys.argv[2], "w"))
