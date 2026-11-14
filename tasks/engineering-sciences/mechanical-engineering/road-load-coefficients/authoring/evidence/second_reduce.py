"""Second, separately written reduction: time-domain trajectory fit of every run at once.

    python second_reduce.py <trial folder> <result.json>

No derivatives are taken. Every run's logged speed is replayed by a fixed-step RK4 integration of

    (m + m_rot) dV/dt = -[F0 (1 + A e^(-c/tau)) r + F2 q (V - h w)|V - h w| + m g h tilt]   (V in km/h -> m/s)

from the mean of its first three logged speeds; F0, F2, w, A, tau and tilt are found by scipy least_squares on the
speed residuals. Start values come from block-mean decelerations (5 s blocks) regressed on 1, h V, V^2 and h.
"""
import csv
import json
import sys

import numpy as np
from scipy.optimize import least_squares

GRAV, KR, RAIR = 9.80665, 0.00864, 287.05
RHO0 = 98210.0 / (RAIR * 293.15)


def hms(text):
    parts = [int(p) for p in text.split(":")]
    return parts[0] * 3600 + parts[1] * 60 + parts[2]


def read(folder):
    sheet = json.load(open(folder + "/vehicle.json"))
    with open(folder + "/run_weather.csv", newline="") as fh:
        wx = {r["run"]: r for r in csv.DictReader(fh)}
    with open(folder + "/coast_log.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    t0 = min(hms(r["release_clock"]) for r in wx.values())
    runs = []
    for key in sorted(wx, key=int):
        mine = [r for r in rows if r["run"] == key]
        temp, kpa = float(wx[key]["air_C"]), float(wx[key]["baro_kPa"])
        runs.append({
            "h": 1.0 if mine[0]["heading"] == "E" else -1.0,
            "t": np.array([float(r["t_s"]) for r in mine]),
            "v": np.array([float(r["speed_kmh"]) for r in mine]) / 3.6,
            "c0": hms(wx[key]["release_clock"]) - t0,
            "r": 1.0 / (1.0 + KR * (temp - 20.0)),
            "q": kpa * 1000.0 / (RAIR * (temp + 273.15)) / RHO0,
        })
    return sheet, runs


class Replay:
    def __init__(self, sheet, runs):
        self.m = float(sheet["test_mass_kg"])
        self.me = self.m + float(sheet["rotating_equivalent_kg"])
        self.runs = runs
        n = max(len(r["t"]) for r in runs)
        self.dt = float(np.median(np.diff(runs[0]["t"])))
        self.n = n
        self.h = np.array([r["h"] for r in runs])
        self.r = np.array([r["r"] for r in runs])
        self.q = np.array([r["q"] for r in runs])
        self.c0 = np.array([r["c0"] for r in runs])
        self.obs = np.full((len(runs), n), np.nan)
        for i, r in enumerate(runs):
            self.obs[i, :len(r["v"])] = r["v"]
        self.mask = ~np.isnan(self.obs)
        self.v0 = np.array([r["v"][:3].mean() + (r["v"][0] - r["v"][2]) / 2.0 for r in runs])

    def accel(self, v, clock, p):
        f0, f2, w, amp, ltau, tilt = p[:6]
        u = v - self.h * w / 3.6
        roll = f0 * (1.0 + amp * np.exp(-clock / np.exp(ltau))) * self.r
        return -(roll + f2 * 12.96 * self.q * u * np.abs(u) + self.m * GRAV * self.h * tilt) / self.me

    def speeds(self, p):
        v = self.v0.copy()
        out = np.empty((len(v), self.n))
        out[:, 0] = v
        dt = self.dt
        for k in range(1, self.n):
            c = self.c0 + (k - 1) * dt
            a1 = self.accel(v, c, p)
            a2 = self.accel(v + 0.5 * dt * a1, c + 0.5 * dt, p)
            a3 = self.accel(v + 0.5 * dt * a2, c + 0.5 * dt, p)
            a4 = self.accel(v + dt * a3, c + dt, p)
            v = v + dt * (a1 + 2 * a2 + 2 * a3 + a4) / 6.0
            out[:, k] = v
        return out

    def residual(self, p):
        return (self.speeds(p) - self.obs)[self.mask] * 3.6


def first_guess(rep):
    rows, ys = [], []
    for i, r in enumerate(rep.runs):
        k = 25
        for j in range(0, len(r["v"]) - k, k):
            dv = (r["v"][j + k] - r["v"][j]) / (r["t"][j + k] - r["t"][j])
            vb = 3.6 * r["v"][j:j + k + 1].mean()
            rows.append([r["r"], r["q"] * r["h"] * vb, r["q"] * vb * vb, rep.m * GRAV * r["h"]])
            ys.append(-rep.me * dv)
    c, *_ = np.linalg.lstsq(np.array(rows), np.array(ys), rcond=None)
    return c[0], c[2], -c[1] / (2.0 * c[2]), c[3]


def reduce(folder):
    sheet, runs = read(folder)
    rep = Replay(sheet, runs)
    f0, f2, w, tilt = first_guess(rep)
    p0 = np.r_[0.85 * f0, f2, w, 0.3, np.log(900.0), tilt]
    lo = np.r_[0.0, 0.0, -80.0, -0.5, np.log(150.0), -0.02]
    hi = np.r_[2000.0, 1.0, 80.0, 3.0, np.log(5000.0), 0.02]
    best = least_squares(rep.residual, p0, bounds=(lo, hi), x_scale="jac", max_nfev=80)
    p = best.x
    return {"F0_N": float(p[0]), "F2_N_per_kmh2": float(p[1]), "wind_kmh": float(p[2]), "excess": float(p[3]),
            "tau_s": float(np.exp(p[4]))}


if __name__ == "__main__":
    json.dump(reduce(sys.argv[1]), open(sys.argv[2], "w"))
