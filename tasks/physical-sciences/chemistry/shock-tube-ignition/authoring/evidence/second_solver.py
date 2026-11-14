"""A second reduction of one campaign, written apart from the reference (SciPy).

    python second_solver.py <campaign folder> <answer.json>

Different choices: the end-wall shock speed from a fit of the exact attenuated arrival law
t(x) = -ln(1 + a x) / (a U0); shock states by brentq; the OH* onset by the steepest tangent of a Savitzky-Golay
smoothed trace; the end-wall history smoothed by Savitzky-Golay; and one nonlinear least-squares fit of
(ln A, alpha, beta, Ea) on ln tau with the Livengood-Wu delay computed along each shot's history at every step.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq, curve_fit, least_squares
from scipy.signal import savgol_filter

GAS = 8.314462618
ATM_PA = 101325.0


def enthalpy(abc, T):
    return abc[0] * T + abc[1] * T * T / 2 + abc[2] * T ** 3 / 3


def states(abc, M, T1, P1, U):
    r1 = P1 * M / (GAS * T1)
    g1 = lambda e: enthalpy(abc, e * M * (P1 + r1 * U * U * (1 - e)) / (r1 * GAS)) - enthalpy(abc, T1) \
        - M * U * U * (1 - e * e) / 2
    grid = np.linspace(0.03, 0.99, 200)
    v = np.array([g1(x) for x in grid])
    i = int(np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))[0])
    e = brentq(g1, grid[i], grid[i + 1])
    P2 = P1 + r1 * U * U * (1 - e)
    T2 = e * M * P2 / (r1 * GAS)
    u2 = U * (1 - e)
    r2 = P2 * M / (GAS * T2)
    g2 = lambda q: enthalpy(abc, q * M * (P2 + r2 * (u2 / (1 - q)) ** 2 * (1 - q)) / (r2 * GAS)) - enthalpy(abc, T2) \
        - M * (u2 / (1 - q)) ** 2 * (1 - q * q) / 2
    grid = np.linspace(0.03, 0.9, 200)
    v = np.array([g2(x) for x in grid])
    i = int(np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))[0])
    q = brentq(g2, grid[i], grid[i + 1])
    P5 = P2 + r2 * (u2 / (1 - q)) ** 2 * (1 - q)
    return q * M * P5 / (r2 * GAS), P5


def tangent_onset(t, y):
    base = np.median(y[t < 0])
    if np.median(y[-20:]) - base < 0.5:
        return None
    sm = savgol_filter(y, 7, 2)
    d = np.gradient(sm, t)
    k = int(np.argmax(np.where(t > 0, d, -np.inf)))
    return t[k] - (sm[k] - base) / d[k]


def main(folder, out):
    f = Path(folder)
    sheet = json.loads((f / "campaign_sheet.json").read_text())
    tab = np.genfromtxt(f / "shots.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    leak = np.loadtxt(f / "leak_check.csv", delimiter=",", skiprows=1)
    rate = np.polyfit(leak[:, 0], leak[:, 1], 1)[0]
    xs = np.array(sheet["station_positions_from_end_wall_m"])
    sp = sheet["species"]
    data = []
    for row in tab:
        t_arr = np.array([row[n] for n in tab.dtype.names if n.startswith("arrival_")]) * 1e-6
        law = lambda x, a, U0: -np.log1p(a * x) / (a * U0)
        (a_fit, U0), _ = curve_fit(law, xs, t_arr, p0=(0.005, (xs[0] - xs[-1]) / (t_arr[-1] - t_arr[0])))
        P1 = row["P1_kPa"] * 1e3
        air = rate * row["isolated_s"] / P1
        x = {"fuel": row["fuel_mole_fraction"] * (1 - air), "O2": row["o2_mole_fraction"] * (1 - air) + 0.2095 * air,
             "N2": 0.7905 * air}
        x["Ar"] = 1.0 - sum(x.values())
        abc = sum(v * np.array(sp[k]["cp_J_mol_K"]) for k, v in x.items())
        M = sum(v * sp[k]["molar_mass_kg_mol"] for k, v in x.items())
        T5, P5 = states(abc, M, row["T1_K"], P1, U0)
        tr = np.loadtxt(f / "traces" / f"shot_{int(row['shot']):02d}.csv", delimiter=",", skiprows=1)
        tt, pe, oh = tr[:, 0], tr[:, 2], tr[:, 3]
        on = tangent_onset(tt, oh)
        if on is None:
            continue
        keep = (tt >= 2) & (tt < max(on - 30, 15))
        n_keep = int(keep.sum())
        ps = savgol_filter(pe[keep], min(9, n_keep if n_keep % 2 else n_keep - 1), 1) if n_keep >= 3 else pe[keep]
        lev = np.mean(pe[(tt > 5) & (tt < min(on - 30, 60))]) if ((tt > 5) & (tt < min(on - 30, 60))).any() else ps[0]
        cp5 = abc[0] + abc[1] * T5 + abc[2] * T5 ** 2
        data.append((T5, P5, x["O2"], on, np.append(tt[keep], 3e4), np.append(ps, ps[-1]) / lev, cp5 / (cp5 - GAS)))

    grid = np.arange(0.0, 20000.0, 2.0)

    def predicted(theta, d):
        lnA, al, be, E = theta
        T5, P5, xo, _, ht, hr, g = d
        r = np.interp(grid, ht, hr)
        T = T5 * r ** ((g - 1) / g)
        inv = 1.0 / np.exp(lnA + al * np.log(P5 * r / ATM_PA) + be * math.log(xo) + E / (GAS * T))
        cum = cumulative_trapezoid(inv, grid, initial=0.0)
        k = int(np.searchsorted(cum, 1.0))
        if k <= 0 or k >= grid.size:
            return grid[min(max(k, 1), grid.size - 1)]
        return grid[k - 1] + (1.0 - cum[k - 1]) / (cum[k] - cum[k - 1]) * (grid[k] - grid[k - 1])

    def resid(theta):
        return np.array([math.log(predicted(theta, d)) - math.log(d[3]) for d in data])

    X = np.column_stack([np.ones(len(data)), [math.log(d[1] / ATM_PA) for d in data], [math.log(d[2]) for d in data],
                         [1 / d[0] for d in data]])
    start = np.linalg.lstsq(X, np.log([d[3] for d in data]), rcond=None)[0]
    theta0 = np.array([start[0], start[1], start[2], start[3] * GAS])
    fit = least_squares(resid, theta0, x_scale=np.array([1.0, 0.1, 0.1, 1e4]), diff_step=1e-4)
    lnA, al, be, E = fit.x
    pts = sheet["report_points"]
    tau = [math.exp(lnA + al * math.log(q["P_atm"]) + be * math.log(q["x_o2"]) + E / (GAS * q["T_K"])) for q in pts]
    Path(out).write_text(json.dumps({"activation_energy_kj_mol": E / 1e3, "pressure_exponent": al,
                                     "oxygen_exponent": be, "tau_report_1_us": tau[0], "tau_report_2_us": tau[1]},
                                    indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
