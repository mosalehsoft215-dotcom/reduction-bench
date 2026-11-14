"""Neutral temperature and longitudinal force of continuous welded rail from one lift test.

    python3 stress_free.py <lift folder> <result file>

1. Section: the catalogue section of the profile less a flat layer of the measured vertical head wear across the
   head width (parallel-axis shift of the neutral axis); weight per metre reduced by the removed steel.
2. Force: the load cell's factor is taken from the 500 kg check mass hung on it at the start of the test
   (factor = 500 kg x g / mean check reading); jack force = factor x reading.
3. Model: Euler-Bernoulli rail with axial force P and its weight on a rigid bed; the lifted half-length a satisfies
   w'(0) = 0, w(a) = w'(a) = w''(a) = 0 and mid-point shear F/2, giving (F, lift) for every a.  P is the value whose
   force-lift curve fits the jack steps best (bounded scalar search, least squares in force).
4. Rail temperature = mean of the probe readings; neutral temperature = T_rail - P / (E A alpha).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

E, ALPHA, G, RHO, CHECK_KG = 210e9, 1.15e-5, 9.80665, 7850.0, 500.0
PROFILES = {"60E1": {"area_cm2": 76.7, "second_moment_cm4": 3038.3, "mass_kg_per_m": 60.21, "height_mm": 172.0, "neutral_axis_above_foot_mm": 80.9, "head_width_mm": 72.0}, "54E1": {"area_cm2": 69.77, "second_moment_cm4": 2337.9, "mass_kg_per_m": 54.77, "height_mm": 159.0, "neutral_axis_above_foot_mm": 74.89999999999999, "head_width_mm": 70.0}, "49E1": {"area_cm2": 62.92, "second_moment_cm4": 1816.0, "mass_kg_per_m": 49.39, "height_mm": 149.0, "neutral_axis_above_foot_mm": 70.2, "head_width_mm": 67.0}}
A_GRID = np.linspace(0.8, 16.0, 260)


def section(profile, wear_mm):
    p = PROFILES[profile]
    A0, I0 = p["area_cm2"] * 1e-4, p["second_moment_cm4"] * 1e-8
    h, yna, bh = p["height_mm"] / 1000, p["neutral_axis_above_foot_mm"] / 1000, p["head_width_mm"] / 1000
    t = wear_mm / 1000
    strip = bh * t
    arm = h - t / 2 - yna
    A1 = A0 - strip
    I1 = I0 - bh * t ** 3 / 12 - strip * arm ** 2 - A1 * (strip * arm / A1) ** 2
    return A1, I1, p["mass_kg_per_m"] - strip * RHO


def homogeneous(x, P, EI):
    if abs(P) < 1e-6:
        return np.array([[1, x, x * x, x ** 3], [0, 1, 2 * x, 3 * x * x], [0, 0, 2, 6 * x], [0, 0, 0, 6]], float)
    k = np.sqrt(abs(P) / EI)
    if P > 0:
        c, s = np.cos(k * x), np.sin(k * x)
        return np.array([[1, x, c, s], [0, 1, -k * s, k * c], [0, 0, -k * k * c, -k * k * s],
                         [0, 0, k ** 3 * s, -k ** 3 * c]])
    c, s = np.cosh(k * x), np.sinh(k * x)
    return np.array([[1, x, c, s], [0, 1, k * s, k * c], [0, 0, k * k * c, k * k * s], [0, 0, k ** 3 * s, k ** 3 * c]])


def forced(x, P, EI, q):
    if abs(P) < 1e-6:
        return np.array([-q * x ** 4 / (24 * EI), -q * x ** 3 / (6 * EI), -q * x * x / (2 * EI), -q * x / EI])
    return np.array([-q * x * x / (2 * P), -q * x / P, -q / P, 0.0])


def curve(P, EI, q):
    lifts, forces = [], []
    H0, f0 = homogeneous(0.0, P, EI), forced(0.0, P, EI, q)
    for a in A_GRID:
        Ha, fa = homogeneous(a, P, EI), forced(a, P, EI, q)
        M = np.zeros((5, 5))
        rhs = np.zeros(5)
        M[0, :4], rhs[0] = H0[1], -f0[1]
        M[1:4, :4], rhs[1:4] = Ha[0:3], -fa[0:3]
        M[4, :4], M[4, 4], rhs[4] = EI * H0[3] + P * H0[1], -0.5, -(EI * f0[3] + P * f0[1])
        sol = np.linalg.solve(M, rhs)
        lifts.append(H0[0] @ sol[:4] + f0[0])
        forces.append(sol[4])
    return np.array(lifts), np.array(forces)


def reduce(folder):
    sheet = json.loads((folder / "release_sheet.json").read_text())
    with open(folder / "jack_steps.csv") as f:
        steps = [(float(r["lift_mm"]) / 1000, float(r["cell_mV"])) for r in csv.DictReader(f)]
    with open(folder / "check_mass.csv") as f:
        check = np.mean([float(r["cell_mV"]) for r in csv.DictReader(f)])
    factor = CHECK_KG * G / check
    lift = np.array([s[0] for s in steps])
    force = factor * np.array([s[1] for s in steps])
    A1, I1, m1 = section(sheet["profile"], sheet["track_geometry"]["vertical_head_wear_mm"])
    EI, q = E * I1, m1 * G

    def misfit(P_kn):
        z, F = curve(P_kn * 1e3, EI, q)
        return float(np.sum((np.interp(lift, z, F) - force) ** 2))

    P = minimize_scalar(misfit, bounds=(-900.0, 900.0), method="bounded", options={"xatol": 0.05}).x * 1e3
    T_rail = float(np.mean(sheet["rail_temperature_C"]))
    return T_rail - P / (E * A1 * ALPHA), P / 1e3


def main():
    t_n, p_kn = reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps({"stress_free_temp_c": round(t_n, 3), "axial_force_kn": round(p_kn, 2)}) + "\n")


if __name__ == "__main__":
    main()
