"""Axial force of a roof member or bar hanger, before and after an adjustment, from two ambient acceleration records.

    python3 canopy_tension.py <member folder> <output file>

1. Spectrum of each record by Welch (Hann, eighth-length segments); peaks with a prominence of at least 12 dB,
   each located by a Lorentzian fitted to the nine spectral lines around it.  Peaks of the structure (roof modes seen through
   the anchorage) and stray peaks are left in: the ladder match below charges an unmatched peak a fixed cap.
2. In-plane model of the member: bending stiffness EI, axial force T, mass m, ends held against rotation or free to
   rotate, and, for inclined strand backstays, the sag-extensibility stiffness of a flat-sagged member.
   Bars: EI = E pi d^4 / 64, no sag (vertical); frequencies from the dimensionless ladder Omega_k(xi),
   xi = Lv sqrt(T / EI), f_k = Omega_k sqrt(EI / (m Lv^4)) / (2 pi); pinned ends in closed form, held ends tabulated
   once by Hermite beam elements.  Backstays: EI of the strand bundle, E A from the metallic area, beam elements with sag.
3. Robust ladder match: every member peak is charged its squared log distance to the nearest model mode (capped), and
   every model mode inside the observed band without a peak is charged a fraction of the cap.  Bars: grid over the
   end condition and the vibrating length Lv (0.80-1.00 of the sheet length; the bar is clamped where it enters its
   anchor blocks), the force of each record by a fine scan of xi, then a joint least-squares polish of
   (T_before, T_after, Lv) on the matched pairs (residuals in Hz, pairs re-made over three rounds).  Backstays: sheet
   length, ends free to rotate (a strand bundle's bending is negligible), coarse and fine scans and the same polish
   of T per record.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import curve_fit, least_squares
from scipy.signal import find_peaks, welch

G = 9.80665
CAP = 0.02
SEGMENTS = 8
STRAND_D = 0.0157
XI_GRID = np.geomspace(0.3, 3000.0, 700)


def beam_matrices(nel, L, EI, T, m):
    h = L / nel
    nd = 2 * (nel + 1)
    K = np.zeros((nd, nd))
    M = np.zeros((nd, nd))
    c = np.zeros(nd)
    hh = h * h
    kb = EI / h ** 3 * np.array([[12, 6 * h, -12, 6 * h], [6 * h, 4 * hh, -6 * h, 2 * hh],
                                 [-12, -6 * h, 12, -6 * h], [6 * h, 2 * hh, -6 * h, 4 * hh]])
    kg = T / (30 * h) * np.array([[36, 3 * h, -36, 3 * h], [3 * h, 4 * hh, -3 * h, -hh],
                                  [-36, -3 * h, 36, -3 * h], [3 * h, -hh, -3 * h, 4 * hh]])
    me = m * h / 420 * np.array([[156, 22 * h, 54, -13 * h], [22 * h, 4 * hh, 13 * h, -3 * hh],
                                 [54, 13 * h, 156, -22 * h], [-13 * h, -3 * hh, -22 * h, 4 * hh]])
    for e in range(nel):
        s = slice(2 * e, 2 * e + 4)
        K[s, s] += kb + kg
        M[s, s] += me
        c[s] += np.array([h / 2, hh / 12, h / 2, -hh / 12])
    return K, M, c


def held_table(kmax=16, nel=100):
    """Omega_k(xi) for ends held against rotation (unit length, EI and m); read from held_table.npy beside this
    file when present (the same table, precomputed), otherwise computed here"""
    cache = Path(__file__).with_name("held_table.npy")
    if cache.is_file():
        stored = np.load(cache)
        if stored.shape == (len(XI_GRID), kmax):
            return stored
    out = np.zeros((len(XI_GRID), kmax))
    for i, xi in enumerate(XI_GRID):
        K, M, _ = beam_matrices(nel, 1.0, 1.0, xi * xi, 1.0)
        nd = K.shape[0]
        keep = [j for j in range(nd) if j not in (0, 1, nd - 2, nd - 1)]
        w = eigh(K[np.ix_(keep, keep)], M[np.ix_(keep, keep)], eigvals_only=True, subset_by_index=[0, kmax - 1])
        out[i] = np.sqrt(w)
    return out


def bar_ladder(T, Lv, EI, m, held, table, kmax=16):
    xi = Lv * np.sqrt(max(T, 1.0) / EI)
    if held:
        lx = np.log(XI_GRID)
        om = np.array([np.interp(np.log(xi), lx, table[:, k]) for k in range(kmax)])
    else:
        k = np.arange(1, kmax + 1) * np.pi
        om = np.sqrt(k ** 4 + xi * xi * k * k)
    return om / (2 * np.pi) * np.sqrt(EI / (m * Lv ** 4))


def stay_ladder(T, L, EI, EA, m, theta, kmax=16, nel=90):
    K, M, c = beam_matrices(nel, L, EI, T, m)
    q = m * G * np.cos(np.radians(theta))
    sag = q * L * L / (8 * T)
    K = K + EA / (L * (1 + 8 * (sag / L) ** 2)) * (q / T) ** 2 * np.outer(c, c)
    nd = K.shape[0]
    keep = [j for j in range(nd) if j not in (0, nd - 2)]
    w = eigh(K[np.ix_(keep, keep)], M[np.ix_(keep, keep)], eigvals_only=True, subset_by_index=[0, kmax - 1])
    return np.sqrt(w) / (2 * np.pi)


def lorentz(x, A, f0, g, c0):
    return A / ((x - f0) ** 2 + g * g) + c0


def peaks(acc, fs):
    f, p = welch(acc - acc.mean(), fs=fs, window="hann", nperseg=len(acc) // SEGMENTS)
    db = 10 * np.log10(p + 1e-30)
    idx, _ = find_peaks(db, prominence=12.0)
    df = f[1] - f[0]
    out = []
    for i in idx:
        if f[i] < 0.12 or i < 5 or i + 5 >= len(f):
            continue
        x, y = f[i - 4:i + 5], p[i - 4:i + 5]
        try:
            prm = curve_fit(lorentz, x, y, p0=(p[i] * df * df, f[i], df, 0.0), maxfev=4000)[0]
            out.append(prm[1] if abs(prm[1] - f[i]) < 2 * df else f[i])
        except RuntimeError:
            out.append(f[i])
    return np.array(out)


def cost(obs, model):
    band = model[(model > obs.min() * 0.97) & (model < obs.max() * 1.03)]
    d = np.abs(np.log(obs[:, None] / model[None, :])).min(axis=1)
    hit = np.abs(np.log(band[:, None] / obs[None, :])).min(axis=1) if len(band) else np.zeros(0)
    return np.sum(np.minimum(d, CAP) ** 2) + 0.35 * CAP ** 2 * np.sum(hit > CAP)


def costs(obs, ladders):
    """cost() for many model ladders at once (rows of ladders)"""
    d = np.abs(np.log(obs[None, :, None] / ladders[:, None, :])).min(axis=2)
    band = (ladders > obs.min() * 0.97) & (ladders < obs.max() * 1.03)
    hit = np.abs(np.log(ladders[:, :, None] / obs[None, None, :])).min(axis=2) > CAP
    return np.sum(np.minimum(d, CAP) ** 2, axis=1) + 0.35 * CAP ** 2 * np.sum(band & hit, axis=1)


def pairs(obs, model, tol=0.008):
    k = np.abs(np.log(obs[:, None] / model[None, :])).argmin(axis=1)
    ok = np.abs(np.log(obs / model[k])) < tol
    return obs[ok], k[ok]


def polish(x0, scale, bounds, ladder_of, obs_sets, res, rounds=3):
    """least squares on matched (peak, mode) pairs, residuals in Hz over the line spacing; pairs re-made each round"""
    x = np.array(x0, float)
    for _ in range(rounds):
        sets = [pairs(obs, ladder_of(x, i), tol=0.01) for i, obs in enumerate(obs_sets)]
        fun = lambda z: np.concatenate([(o - ladder_of(z, i)[k]) / res for i, (o, k) in enumerate(sets)])
        x = least_squares(fun, x0=x, x_scale=scale, bounds=bounds).x
    return x


def load(folder):
    sheet = json.loads((folder / "member_card.json").read_text())
    recs = {}
    for tag, name in sheet["records"].items():
        with open(folder / name) as fh:
            recs[tag] = np.array([float(r["acc_mg"]) for r in csv.DictReader(fh)])
    return sheet, recs


def solve(folder):
    sheet, recs = load(folder)
    fs = sheet["sample_rate_hz"]
    pb, pa = peaks(recs["before"], fs), peaks(recs["after"], fs)
    res = SEGMENTS * fs / len(recs["before"])               # spectral line spacing of the Welch segments
    m = sheet["mass_per_metre_kg"]
    L = sheet["length_between_working_points_m"]
    E = sheet["steel_modulus_mpa"] * 1e6
    member = sheet["member"]
    if "diameter_mm" in member:
        d = member["diameter_mm"] / 1000
        EI = E * np.pi * d ** 4 / 64
        table = held_table()
        rows = (XI_GRID > 0.5) & (XI_GRID < 400)
        xis = XI_GRID[rows]
        kpi = np.arange(1, table.shape[1] + 1) * np.pi
        shapes = {True: table[rows], False: np.sqrt(kpi[None, :] ** 4 + (xis[:, None] * kpi[None, :]) ** 2)}
        best = None
        for held in (True, False):
            for Lv in np.linspace(0.80 * L, L, 41):
                ladders = shapes[held] / (2 * np.pi) * np.sqrt(EI / (m * Lv ** 4))
                total, Ts = 0.0, []
                for obs in (pb, pa):
                    cs = costs(obs, ladders)
                    j = int(np.argmin(cs))
                    total += cs[j]
                    Ts.append((xis[j] / Lv) ** 2 * EI)
                if best is None or total < best[0]:
                    best = (total, held, Lv, Ts)
        _, held, Lv, Ts = best
        x = polish([Ts[0], Ts[1], Lv], [Ts[0], Ts[1], 0.05 * L], ([1e3, 1e3, 0.75 * L], [1e8, 1e8, 1.02 * L]),
                   lambda z, i: bar_ladder(z[i], z[2], EI, m, held, table), (pb, pa), res)
        Tb, Ta = x[:2]
    else:
        EI = member["strands"] * E * np.pi * STRAND_D ** 4 / 64
        EA = sheet["metallic_area_mm2"] * 1e-6 * E
        th = sheet["chord_inclination_deg"]
        out = []
        for obs in (pb, pa):
            spacing = np.median(np.diff(np.sort(obs)))
            T0 = 4 * m * L * L * spacing ** 2
            for lo, hi, num in ((0.15, 6.0, 160), (0.96, 1.04, 41)):
                grid = np.geomspace(lo * T0, hi * T0, num)
                cs = [cost(obs, stay_ladder(T, L, EI, EA, m, th)) for T in grid]
                T0 = grid[int(np.argmin(cs))]
            x = polish([T0], [T0], ([1e3], [1e9]), lambda z, i: stay_ladder(z[0], L, EI, EA, m, th), (obs,), res)
            out.append(x[0])
        Tb, Ta = out
    return Tb / 1e3, Ta / 1e3


def main():
    tb, ta = solve(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps({"axial_before_kn": round(tb, 3), "axial_after_kn": round(ta, 3)}) + "\n")


if __name__ == "__main__":
    main()
