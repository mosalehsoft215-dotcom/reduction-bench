"""Independent second reduction of a coupon (written apart from solution/fmr_reduce.py).

    python second_reduction.py <coupon folder> [answer json]

Differences from the reference:
- the modulation is averaged numerically over 48 phases of the coil cycle (no closed-form first harmonic);
- the output filter is a causal convolution with the impulse response of n equal RC sections,
  t^(n-1) exp(-t/tau) / (tau^n (n-1)!), sampled on the logged grid at 6x resolution;
- the absorption is a Lorentzian derivative pair written in real form (absorption and dispersion weights fitted
  as two linear amplitudes), and the line parameters are found by a coarse grid over centre and width followed by
  variable projection (linear amplitudes and background solved exactly at each step);
- equilibrium and resonance are found with brentq in field space from the energy's analytic derivatives;
- gaussmeter scale and zero, then g, Meff and the normal-mount tilt by Nelder-Mead on the resonance fields;
- alpha and the Meff spread from the normal-mount widths, factors by finite differences in field.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares, minimize

MU_B = 9.2740100783e-24
HBAR = 1.054571817e-34
H = 6.62607015e-34
DPPH_G = 2.0036
SENS_PER_K = -4.5e-4
PHASES = np.linspace(0, 2 * np.pi, 48, endpoint=False)
REFINE = 6


def read_coupon(folder):
    folder = Path(folder)
    out = []
    with open(folder / "register.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            data = np.genfromtxt(folder / "sweeps" / (row["sweep"] + ".csv"), delimiter=",", names=True)
            row = {k: (v if k in ("sweep", "geometry", "kind", "direction") else float(v)) for k, v in row.items()}
            row["B"] = data["reading_mT"] / 1000.0
            row["S"] = data["signal_uV"]
            out.append(row)
    return out


def rc_kernel(n, tau, dt):
    if tau <= 0:
        return np.array([1.0])
    t = np.arange(0, 30 * n * tau, dt)
    k = t ** (n - 1) * np.exp(-t / tau)
    return k / k.sum()


def basis(Bf, centre, hw, amp_mod):
    """Absorption-like and dispersion-like first-harmonic responses on a fine field grid."""
    x = Bf[:, None] + amp_mod * np.cos(PHASES)[None, :]
    d = centre - x
    absn = hw / (d * d + hw * hw)
    disp = d / (d * d + hw * hw)
    c = np.cos(PHASES)[None, :]
    return (absn * c).mean(axis=1), (disp * c).mean(axis=1)


def line(row):
    B, S = row["B"], row["S"]
    n = len(B)
    Bf = np.linspace(B[0], B[-1], (n - 1) * REFINE + 1)
    dt = abs(Bf[1] - Bf[0]) / (row["rate_mT_s"] / 1000.0)
    kern = rc_kernel(int(round(row["slope_dB_oct"] / 6)), row["tau_ms"] / 1000.0, dt)
    amp_mod = row["mod_mT"] / 1000.0
    up = row["direction"] == "up"
    lin = (B - B.mean()) / (B[-1] - B[0])

    def columns(centre, hw):
        cols = []
        for f in basis(Bf, centre, hw, amp_mod):
            if amp_mod <= 0:
                f = np.gradient(f, Bf)
            seq = f if up else f[::-1]
            pad = np.concatenate([np.full(len(kern), seq[0]), seq])
            sm = np.convolve(pad, kern)[len(kern):len(kern) + len(seq)]
            sm = sm if up else sm[::-1]
            cols.append(sm[::REFINE])
        return np.column_stack(cols + [np.ones(n), lin])

    def cost(p):
        A = columns(p[0], abs(p[1]))
        coef, *_ = np.linalg.lstsq(A, S, rcond=None)
        return A @ coef - S

    span = B[-1] - B[0]
    best = None
    for c in np.linspace(B[0] + 0.2 * span, B[-1] - 0.2 * span, 25):
        for w in span * np.array([0.01, 0.025, 0.05, 0.08]):
            r = float(np.sum(cost([c, w]) ** 2))
            if best is None or r < best[0]:
                best = (r, c, w)
    fit = least_squares(cost, [best[1], best[2]], x_scale=[best[2], best[2]])
    centre, hw = fit.x[0], abs(fit.x[1])
    pos = (centre - B[0]) / span
    pos = min(max(pos if up else 1 - pos, 0), 1)
    temp = row["probe_C_start"] + pos * (row["probe_C_end"] - row["probe_C_start"])
    return centre, hw, temp


def equilibrium(B, tb, meff):
    f = lambda t: B * math.sin(t - tb) - 0.5 * meff * math.sin(2 * t)
    if tb >= math.pi / 2 - 1e-12:
        return math.pi / 2
    if f(tb + 1e-12) >= 0:
        return tb
    return brentq(f, tb + 1e-12, math.pi / 2, xtol=1e-15)


def det(B, tb, meff):
    t = equilibrium(B, tb, meff)
    c = math.cos(t - tb)
    return (B * c - meff * math.cos(2 * t)) * (B * c - meff * math.cos(t) ** 2)


def resonance(freq, g, meff, tb):
    w2 = (2 * math.pi * freq * HBAR / (g * MU_B)) ** 2
    grid = np.linspace(0.002, 3.2, 160)
    vals = np.array([det(b, tb, meff) for b in grid])
    i0 = int(np.argmin(vals))
    above = [i for i in range(i0, len(grid)) if vals[i] >= w2]
    if not above:
        return 3.2
    i = above[0]
    if i == i0:
        return grid[i0]
    return brentq(lambda b: det(b, tb, meff) - w2, grid[i - 1], grid[i], xtol=1e-12)


def film_constants(folder):
    rows = read_coupon(folder)
    for r in rows:
        r["c"], r["w"], r["T"] = line(r)
    mk = [r for r in rows if r["kind"] == "marker"]
    Btrue = np.array([H * r["freq_GHz"] * 1e9 / (DPPH_G * MU_B) for r in mk])
    fac = np.array([1 + SENS_PER_K * (r["T"] - 25) for r in mk])
    M = np.column_stack([fac * Btrue, np.ones(len(mk))])
    s25, zero = np.linalg.lstsq(M, np.array([r["c"] for r in mk]), rcond=None)[0]
    films = [r for r in rows if r["kind"] == "film"]
    for r in films:
        s = s25 * (1 + SENS_PER_K * (r["T"] - 25))
        r["Br"], r["HW"] = (r["c"] - zero) / s, r["w"] / s

    def mismatch(p):
        g, meff, tilt = p
        tot = 0.0
        for r in films:
            tb = math.pi / 2 if r["geometry"] == "in-plane" else abs(tilt)
            tot += ((resonance(r["freq_GHz"] * 1e9, g, meff, tb) - r["Br"]) * 1e4) ** 2
        return tot

    ip = [r for r in films if r["geometry"] == "in-plane"]
    f_ip = np.array([r["freq_GHz"] * 1e9 for r in ip])
    B_ip = np.array([r["Br"] for r in ip])
    # in-plane start: (2 pi f)^2 = gamma^2 B (B + Meff) is linear in (gamma^2, gamma^2 Meff)
    a, b = np.linalg.lstsq(np.column_stack([B_ip ** 2, B_ip]), (2 * np.pi * f_ip) ** 2, rcond=None)[0]
    g0 = math.sqrt(a) * HBAR / MU_B
    start = [g0, b / a, 0.003]
    best = None
    for t0 in (0.0, 0.006, 0.012):
        start[2] = t0
        res = minimize(mismatch, start, method="Nelder-Mead",
                       options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000, "initial_simplex": None})
        if best is None or res.fun < best.fun:
            best = res
    g, meff, tilt = best.x
    tilt = abs(tilt)

    nm = [r for r in films if r["geometry"] == "normal"]
    rows_A = []
    for r in nm:
        f = r["freq_GHz"] * 1e9
        br = resonance(f, g, meff, tilt)
        h = 2e-6
        dD = (det(br + h, tilt, meff) - det(br - h, tilt, meff)) / (2 * h)
        t = equilibrium(br, tilt, meff)
        c = math.cos(t - tilt)
        trace = 2 * br * c - meff * math.cos(2 * t) - meff * math.cos(t) ** 2
        w = 2 * math.pi * f * HBAR / (g * MU_B)
        dm = 1e-5
        dB = (resonance(f, g, meff + dm, tilt) - resonance(f, g, meff - dm, tilt)) / (2 * dm)
        rows_A.append([w * trace / abs(dD), abs(dB)])
    alpha, spread = np.linalg.lstsq(np.array(rows_A), np.array([r["HW"] for r in nm]), rcond=None)[0]
    return {"damping_alpha": float(alpha), "spectroscopic_g": float(g), "mu0_meff_T": float(meff)}


if __name__ == "__main__":
    out = film_constants(sys.argv[1])
    print(json.dumps(out))
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(json.dumps(out) + "\n")
