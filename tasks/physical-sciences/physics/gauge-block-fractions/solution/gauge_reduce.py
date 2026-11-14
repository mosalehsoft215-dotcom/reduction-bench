"""Reduces one gauge-block calibration to its central length at 20 C (deviation from nominal) and its linear
expansion coefficient.

    python3 gauge_reduce.py <block directory> <result.json>

1. Wavelengths: a certificate wavelength given in standard air is turned into its vacuum value with the refractive
   index of standard air at that wavelength (fixed-point iteration of the notes' Edlen equation).
2. Lengths: for every stored reading, the refractive index of each laser at the logged air temperature, pressure,
   humidity and CO2; the integer orders by the method of exact fractions, all red orders that put the length within
   3 um of nominal (first reading) or 1.5 um of the previous reading's length tried and the one whose predicted green
   and orange fractions match best kept;
   L = (N + f) lambda_vac / (2 n), averaged over the three lasers.
3. Block temperature: the block follows the platen as a first-order lag; the platen PRT record, taken as linear
   between readings, is passed through the lag (exact for linear pieces) from equilibrium at the first level.
4. Fit: L_k = L20 (1 + alpha (T_block,k - 20)) by least squares over L20, alpha and the time constant.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def edlen(lam_nm, t, p, rh, co2_ppm):
    s2 = (1e3 / lam_nm) ** 2
    ns = 1e-8 * (8342.54 + 2406147.0 / (130.0 - s2) + 15998.0 / (38.9 - s2))
    ns *= 1 + 0.5327 * (co2_ppm * 1e-6 - 0.00045)
    ntp = 1 + p * ns / 96095.43 * (1 + 1e-8 * (0.601 - 0.00972 * t) * p) / (1 + 0.0036610 * t)
    f = rh / 100.0 * 611.2 * np.exp(17.62 * t / (243.12 + t))
    return ntp - f * (3.7345 - 0.0401 * s2) * 1e-10


def vacuum_wavelength(ch, convert=True):
    lam = ch["wavelength_nm"]
    if not convert or ch["medium"].startswith("vacuum"):
        return lam
    vac = lam
    for _ in range(5):                                   # lambda_vac = lambda_std n_std(lambda_vac)
        vac = lam * edlen(vac, 15.0, 101325.0, 0.0, 450.0)
    return vac


def exact_fractions(fr, lams, ns, nominal_nm, window_nm=3000.0):
    """Length (nm) whose orders match the three excess fractions best, red orders within the window tried."""
    half = [lam / (2 * n) for lam, n in zip(lams, ns)]
    lo, hi = nominal_nm - window_nm, nominal_nm + window_nm
    orders = np.arange(np.floor(lo / half[0] - fr[0]), np.ceil(hi / half[0] - fr[0]) + 1)
    cand = (orders + fr[0]) * half[0]
    cost = np.zeros_like(cand)
    for f, h in zip(fr[1:], half[1:]):
        d = cand / h - f
        cost += (d - np.round(d)) ** 2
    k = int(np.argmin(cost))
    length = cand[k]
    return float(np.mean([(np.round(length / h - f) + f) * h for f, h in zip(fr, half)]))


def block_temperature(minute, prt, tau):
    """First-order lag of the platen PRT (linear between readings), the block at equilibrium at the start."""
    out = np.empty_like(prt)
    out[0] = prt[0]
    for k in range(1, prt.size):
        dt = minute[k] - minute[k - 1]
        s = (prt[k] - prt[k - 1]) / dt
        a = np.exp(-dt / tau)
        out[k] = prt[k] - s * tau + (out[k - 1] - prt[k - 1] + s * tau) * a
    return out


def reduce_block(folder, lag=True, convert=True, last=None):
    """Options other than the defaults are the shortcuts tabulated by authoring/evidence/bench.py."""
    folder = Path(folder)
    sheet = json.loads((folder / "block.json").read_text())
    with open(folder / "readings.csv", newline="") as fh:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    lams = [vacuum_wavelength(ch, convert) for ch in sheet["lasers"]]
    cols = ["frac_" + ch["laser"] for ch in sheet["lasers"]]
    nominal = sheet["nominal_length_mm"] * 1e6
    minute = np.array([r["minute"] for r in rows])
    prt = np.array([r["platen_prt_c"] for r in rows])
    length, centre = [], nominal                         # each reading searched around the previous length
    for r in rows:
        ns = [edlen(lam, r["air_c"], r["pressure_pa"], r["humidity_pct"], r["co2_ppm"]) for lam in lams]
        centre = exact_fractions([r[c] for c in cols], lams, ns, centre, 3000.0 if not length else 1500.0)
        length.append(centre)
    length = np.array(length)
    if last:                                             # shortcut: the last readings of every hold, PRT as block
        hold = sheet["hold_min"]
        level = np.floor((minute - 1e-6) / hold)
        sel = np.zeros(minute.size, bool)
        for lv in np.unique(level):
            sel[np.flatnonzero(level == lv)[-last:]] = True
        a = np.column_stack([np.ones(sel.sum()), prt[sel] - 20.0])
        c0, c1 = np.linalg.lstsq(a, length[sel], rcond=None)[0]
        return {"length_deviation_nm": float(c0 - nominal), "expansion_coefficient_per_mk": float(c1 / c0 * 1e6)}

    def model(q):
        l20, alpha, tau = q
        temp = block_temperature(minute, prt, tau) if lag else prt
        return l20 * (1 + alpha * 1e-6 * (temp - 20.0))
    sol = least_squares(lambda q: model(q) - length, [nominal, 10.0, 20.0],
                        bounds=([nominal - 1e4, 0.0, 0.5], [nominal + 1e4, 40.0, 300.0]), x_scale=[10.0, 0.1, 1.0])
    return {"length_deviation_nm": float(sol.x[0] - nominal), "expansion_coefficient_per_mk": float(sol.x[1])}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_block(sys.argv[1])))
