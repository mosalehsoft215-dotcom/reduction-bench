"""Independent reduction of one lift test (written separately from solution/stress_free.py).

Differences: the lifted rail is solved as a boundary-value problem with scipy's solve_bvp on the scaled span
[0, 1] (the mid-point force is an unknown parameter), the force-lift curve is built from 30 lift-off lengths, the
axial force is found by a 50 kN grid scan refined by a bounded search around the best node, the check-mass factor uses the median reading, and the
worn section is assembled from first moments about the rail foot.

    python second_lift.py <lift folder> <result file> <rail_profiles.json>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

E, ALPHA, G = 210e9, 1.15e-5, 9.80665


def worn(p, wear_mm):
    A0, I0 = p["area_cm2"] * 1e-4, p["second_moment_cm4"] * 1e-8
    yna, h, b = p["neutral_axis_above_foot_mm"] / 1e3, p["height_mm"] / 1e3, p["head_width_mm"] / 1e3
    t = wear_mm / 1e3
    Q0 = A0 * yna                                   # first moment about the foot
    I_foot = I0 + A0 * yna ** 2                     # second moment about the foot
    a_cut, y_cut = b * t, h - t / 2
    A = A0 - a_cut
    yc = (Q0 - a_cut * y_cut) / A
    I = I_foot - (b * t ** 3 / 12 + a_cut * y_cut ** 2) - A * yc ** 2
    return A, I, p["mass_kg_per_m"] - a_cut * 7850.0


def lift_for(a, P, EI, q, guess=None):
    # y = [w, w', w'', w'''] on s in [0, 1], x = a s ; parameter F
    def ode(s, y, p):
        return np.vstack([a * y[1], a * y[2], a * y[3], a * (-(q + P * y[2]) / EI)])

    def bc(ya, yb, p):
        F = p[0]
        return np.array([ya[1], EI * ya[3] + P * ya[1] - F / 2, yb[0], yb[1], yb[2]])

    s = np.linspace(0, 1, 30)
    y0 = np.zeros((4, s.size)) if guess is None else guess
    res = solve_bvp(ode, bc, s, y0, p=[q * a], tol=1e-8, max_nodes=5000)
    return res.p[0], res.sol(0.0)[0], res.sol(s)


def curve(P, EI, q):
    zs, Fs, guess = [], [], None
    for a in np.linspace(1.0, 15.0, 30):
        F, z, guess = lift_for(a, P, EI, q, guess)
        zs.append(z)
        Fs.append(F)
    return np.array(zs), np.array(Fs)


def main(folder, out, profiles):
    folder = Path(folder)
    sheet = json.loads((folder / "release_sheet.json").read_text())
    prof = json.loads(Path(profiles).read_text())[sheet["profile"]]
    with open(folder / "check_mass.csv") as f:
        k = 500.0 * G / np.median([float(r["cell_mV"]) for r in csv.DictReader(f)])
    with open(folder / "jack_steps.csv") as f:
        rows = list(csv.DictReader(f))
    z_meas = np.array([float(r["lift_mm"]) for r in rows]) / 1e3
    F_meas = k * np.array([float(r["cell_mV"]) for r in rows])
    A, I, m = worn(prof, sheet["track_geometry"]["vertical_head_wear_mm"])
    EI, q = E * I, m * G

    def cost(Pk):
        z, F = curve(Pk * 1e3, EI, q)
        return float(np.mean((np.interp(z_meas, z, F) - F_meas) ** 2))

    grid = np.arange(-700.0, 701.0, 50.0)
    c = [cost(p) for p in grid]
    g = float(grid[int(np.argmin(c))])
    from scipy.optimize import minimize_scalar
    best = minimize_scalar(cost, bounds=(g - 50.0, g + 50.0), method="bounded", options={"xatol": 0.1}).x
    P = float(best) * 1e3
    T = float(np.median(sheet["rail_temperature_C"]))
    Path(out).write_text(json.dumps({"stress_free_temp_c": T - P / (E * A * ALPHA), "axial_force_kn": P / 1e3}))


if __name__ == "__main__":
    main(*sys.argv[1:4])
