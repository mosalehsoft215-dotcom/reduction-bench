"""Road-load coefficients from one coast-down session.

    python3 road_load.py <trial folder> <result.json>

Each run's speed is differentiated with a Savitzky-Golay filter, so every logged sample gives one force balance

    -(m + m_rot) dv/dt = F0 (1 + A e^(-t/tau)) / (1 + K (T - 20)) + F2 (rho / rho_std) U |U| + m g h tilt

with U = V - h w (km/h) the air speed along the car, h = +1 eastbound / -1 westbound, w the along-track wind and t
the session clock since the first release. For trial values of (w, tau) the balance is linear in F0, F0 A, F2 and
tilt; (w, tau) minimise the residual over a grid followed by a bounded Nelder-Mead polish (tau 300-4000 s).
"""
import csv
import json
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import minimize
from scipy.signal import savgol_filter

G0, K_ROLL, R_DRY = 9.80665, 0.00864, 287.05
RHO_STD = 98210.0 / (R_DRY * 293.15)
KMH = 3.6
WINDOW, TRIM = 31, 16
TAU_RANGE = (np.log(300.0), np.log(4000.0))     # settling of tyres and driveline, s


def clock_s(stamp):
    h, m, s = (int(p) for p in stamp.split(":"))
    return 3600 * h + 60 * m + s


def gather(folder):
    sheet = json.load(open(f"{folder}/vehicle.json"))
    meta = {}
    with open(f"{folder}/run_weather.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            meta[row["run"]] = (clock_s(row["release_clock"]), float(row["air_C"]), float(row["baro_kPa"]))
    first = min(m[0] for m in meta.values())
    traces = defaultdict(list)
    with open(f"{folder}/coast_log.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            traces[row["run"]].append((row["heading"], float(row["t_s"]), float(row["speed_kmh"])))
    cols = defaultdict(list)
    for run, rows in traces.items():
        t = np.array([r[1] for r in rows])
        v = np.array([r[2] for r in rows]) / KMH
        dt = float(np.median(np.diff(t)))
        start, air_c, kpa = meta[run]
        rho = kpa * 1000.0 / (R_DRY * (air_c + 273.15))
        n = len(t) - 2 * TRIM
        cols["acc"].append(savgol_filter(v, WINDOW, 2, deriv=1, delta=dt)[TRIM:-TRIM])
        cols["V"].append(savgol_filter(v, WINDOW, 2)[TRIM:-TRIM] * KMH)
        cols["clock"].append(t[TRIM:-TRIM] + start - first)
        cols["h"].append(np.full(n, 1.0 if rows[0][0] == "E" else -1.0))
        cols["roll"].append(np.full(n, 1.0 / (1.0 + K_ROLL * (air_c - 20.0))))
        cols["dens"].append(np.full(n, rho / RHO_STD))
    return sheet, {k: np.concatenate(v) for k, v in cols.items()}


class Balance:
    def __init__(self, sheet, s):
        self.s = s
        mass = float(sheet["test_mass_kg"])
        self.force = -(mass + float(sheet["rotating_equivalent_kg"])) * s["acc"]
        self.slope_col = mass * G0 * s["h"]

    def fit(self, w_kmh, tau):
        U = self.s["V"] - self.s["h"] * w_kmh
        design = np.column_stack([self.s["roll"], self.s["roll"] * np.exp(-self.s["clock"] / tau),
                                  self.s["dens"] * U * np.abs(U), self.slope_col])
        coef, *_ = np.linalg.lstsq(design, self.force, rcond=None)
        return coef, float(np.sum((design @ coef - self.force) ** 2))

    def cost(self, p):
        return self.fit(p[0], np.exp(p[1]))[1]


def reduce(folder):
    sheet, s = gather(folder)
    bal = Balance(sheet, s)
    trials = [(w, lt) for w in np.arange(-40.0, 40.1, 2.0) for lt in np.log([300, 600, 1000, 1600, 2500, 4000])]
    start = min(trials, key=bal.cost)
    best = minimize(bal.cost, start, method="Nelder-Mead", bounds=[(-60.0, 60.0), TAU_RANGE],
                    options={"xatol": 1e-4, "fatol": 1e-9, "maxiter": 600})
    w, tau = best.x[0], float(np.exp(best.x[1]))
    coef, _ = bal.fit(w, tau)
    return {"F0_N": float(coef[0]), "F2_N_per_kmh2": float(coef[2]),
            "wind_toward_east_kmh": float(w), "settle_s": tau, "warm_excess": float(coef[1] / coef[0])}


if __name__ == "__main__":
    out = reduce(sys.argv[1])
    with open(sys.argv[2], "w") as fh:
        json.dump(out, fh)
