"""Zero-stress temperature and longitudinal load of continuous welded rail from one jacking run.

    python3 zero_stress.py <run folder> <result file>

The rail is a beam-column: EI w'''' + P w'' = -p, P the longitudinal load (compression positive), p the net downward
load per metre on the lifted rail, w the rise above the bed.  The record does not say how far the fastenings were
undone nor how the track is aligned, so both are read from the force-rise curve itself:
- p is fitted, not taken from the catalogue: on a rail bent to a vertical crest or sag the load P adds a radial
  share to the weight, and that share is not known beforehand.
- While the lift-off point is free, the rail leaves the bed smoothly (w = w' = w'' = 0 there); once it reaches the
  first fastened sleeper at distance L from the jack the rail is held there (w = w' = 0) and the force grows along a
  straight, stiffer branch.  L is fitted as well (a value beyond every lift-off length means no holding).
For given (P, p) the free branch is solved for a vector of lift-off lengths in one batched linear solve; the held
branch is linear in the rise (two solves).  (P, p, L) come from bounded least squares on the jack forces, started
from a small grid; T_zero = mean web-probe temperature - P / (E A alpha).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

E, ALPHA, G = 210e9, 1.15e-5, 9.80665
PROFILES = {"60E1": (76.7, 3038.3, 60.21), "54E1": (69.77, 2337.9, 54.77), "49E1": (62.92, 1816.0, 49.39)}
SPANS = np.linspace(0.6, 26.0, 360)


def mode_rows(x, P, EI):
    """(N, 4, 4): rows w, w', w'', w''' of the four homogeneous solutions at each x"""
    x = np.atleast_1d(np.asarray(x, float))
    one, nil = np.ones_like(x), np.zeros_like(x)
    if abs(P) < 1e-6:
        rows = [[one, x, x * x, x ** 3], [nil, one, 2 * x, 3 * x * x], [nil, nil, 2 * one, 6 * x], [nil, nil, nil, 6 * one]]
    else:
        k = np.sqrt(abs(P) / EI)
        if P > 0:
            c, s = np.cos(k * x), np.sin(k * x)
            rows = [[one, x, c, s], [nil, one, -k * s, k * c], [nil, nil, -k * k * c, -k * k * s],
                    [nil, nil, k ** 3 * s, -k ** 3 * c]]
        else:
            c, s = np.cosh(k * x), np.sinh(k * x)
            rows = [[one, x, c, s], [nil, one, k * s, k * c], [nil, nil, k * k * c, k * k * s], [nil, nil, k ** 3 * s, k ** 3 * c]]
    return np.moveaxis(np.array(rows), -1, 0)


def load_rows(x, P, EI, p):
    """(N, 4): particular solution w, w', w'', w''' for the uniform load p"""
    x = np.atleast_1d(np.asarray(x, float))
    if abs(P) < 1e-6:
        return np.stack([-p * x ** 4 / (24 * EI), -p * x ** 3 / (6 * EI), -p * x * x / (2 * EI), -p * x / EI], axis=1)
    return np.stack([-p * x * x / (2 * P), -p * x / P, -p / P * np.ones_like(x), np.zeros_like(x)], axis=1)


def systems(ends, P, EI, p, end_rows, rise=None):
    """batched 5x5 systems (unknowns: 4 constants and the jack force F) for every end position"""
    n = len(ends)
    H0, h0 = mode_rows(0.0, P, EI)[0], load_rows(0.0, P, EI, p)[0]
    He, he = mode_rows(ends, P, EI), load_rows(ends, P, EI, p)
    A, b = np.zeros((n, 5, 5)), np.zeros((n, 5))
    A[:, 0, :4], b[:, 0] = H0[1], -h0[1]
    A[:, 1, :4], A[:, 1, 4], b[:, 1] = EI * H0[3] + P * H0[1], -0.5, -(EI * h0[3] + P * h0[1])
    for j, r in enumerate(end_rows):
        A[:, 2 + j, :4], b[:, 2 + j] = He[:, r], -he[:, r]
    if rise is not None:
        A[:, 4, :4], b[:, 4] = H0[0], rise - h0[0]
    sol = np.linalg.solve(A, b[..., None])[..., 0]
    return sol[:, 4], sol[:, :4] @ H0[0] + h0[0]


class Jacking:
    def __init__(self, folder):
        sheet = json.loads((folder / "jacking_sheet.json").read_text())
        with open(folder / "lift_record.csv") as f:
            pairs = [(float(r["foot_rise_mm"]) / 1000.0, float(r["cell_mv"])) for r in csv.DictReader(f)]
        self.rise = np.array([a for a, _ in pairs])
        self.force = sheet["load_cell"]["factor_N_per_mV"] * np.array([b for _, b in pairs])
        area, inertia, mass = PROFILES[sheet["rail_profile"]]
        self.area, self.EI, self.weight = area * 1e-4, E * inertia * 1e-8, mass * G
        self.t_rail = float(np.mean(sheet["web_probe_C"]))

    def predict(self, P, p, L):
        F, z = systems(SPANS, P, self.EI, p, (0, 1, 2))
        order = np.argsort(z)
        out = np.interp(self.rise, z[order], F[order])
        held = np.interp(self.rise, z[order], SPANS[order]) > L
        if held.any():
            f0, _ = systems([L], P, self.EI, p, (0, 1), rise=0.0)
            f1, _ = systems([L], P, self.EI, p, (0, 1), rise=1.0)
            out[held] = f0[0] + self.rise[held] * (f1[0] - f0[0])
        return out

    def residual(self, x):
        return (self.predict(x[0] * 1e3, x[1], x[2]) - self.force) / 10.0

    def solve(self):
        q = self.weight
        best = None
        for p_kn in (-300.0, 0.0, 300.0):
            for L in (7.0, 10.0, 40.0):
                fit = least_squares(self.residual, x0=[p_kn, q, L], x_scale=[100.0, 0.1 * q, 2.0], diff_step=1e-3,
                                    bounds=([-900.0, 0.3 * q, 3.0], [900.0, 2.0 * q, 60.0]))
                if best is None or fit.cost < best.cost:
                    best = fit
        return best.x


def reduce(folder):
    run = Jacking(Path(folder))
    p_kn = run.solve()[0]
    return run.t_rail - p_kn * 1e3 / (E * run.area * ALPHA), p_kn


def main():
    t_zero, p_kn = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps({"zero_stress_temp_c": round(t_zero, 3), "longitudinal_load_kn": round(p_kn, 2)}) + "\n")


if __name__ == "__main__":
    main()
