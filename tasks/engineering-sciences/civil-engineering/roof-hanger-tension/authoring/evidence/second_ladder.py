"""Independent evaluation of one member survey (does not import the reference).

Differences from the reference: spectra by Welch with sixth-length segments (default window), peaks of 12 dB
prominence at least three lines apart, located by a parabola through the square-root spectrum; analytic frequency equations
instead of beam elements:
  bars held at both ends (beam-string):  with a^2, b^2 = (sqrt(T^2 + 4 EI m w^2) +- T) / (2 EI),
      2 a b (sech(aL) - cos(bL)) + (a^2 - b^2) tanh(aL) sin(bL) = 0, roots by sign scan and brentq;
  bars free to rotate:                   closed form;
  backstays (strands, bending neglected):    Irvine's equations, antisymmetric f = n sqrt(T/m) / (2L), symmetric
      tan(W/2) = W/2 - (4 / lambda^2)(W/2)^3 (solved as sin - cos * (...) = 0 between k pi and (k+2) pi), W = w L sqrt(m/T), lambda^2 = (m g L cos(theta) / T)^2 L E A / (T Le).
Peaks are matched to the nearest model frequency within 1 %; the unknowns (one force per record; for bars also the
end condition and a common vibrating length) are found by a nested scan of a trimmed misfit, then Nelder-Mead.

    python second_ladder.py <member folder> <output file>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.signal import find_peaks, welch

G = 9.80665


def spectrum_peaks(acc, fs):
    f, p = welch(acc - np.mean(acc), fs=fs, nperseg=len(acc) // 6)
    a = np.sqrt(p)
    idx, _ = find_peaks(10 * np.log10(p + 1e-30), prominence=12.0, distance=3)
    out = []
    for i in idx:
        if f[i] > 0.12 and 0 < i < len(f) - 1:
            y0, y1, y2 = a[i - 1], a[i], a[i + 1]
            out.append(f[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (f[1] - f[0]))
    return np.array(out)


def held_freqs(T, L, EI, m, fmax):
    def g(w):
        root = np.sqrt(T * T + 4 * EI * m * w * w)
        a, b = np.sqrt((root + T) / (2 * EI)), np.sqrt((root - T) / (2 * EI))
        aL = np.minimum(a * L, 600.0)
        return (2 * a * b * (2 / (np.exp(aL) + np.exp(-aL)) - np.cos(b * L)) + (a * a - b * b) * np.tanh(aL) * np.sin(b * L)) / (a * a + b * b)
    ws = np.linspace(1e-3, 2 * np.pi * fmax * 1.1, 1500)
    gs = g(ws)
    idx = np.nonzero(gs[:-1] * gs[1:] < 0)[0]
    return np.array([brentq(g, ws[i], ws[i + 1], xtol=1e-12) for i in idx]) / (2 * np.pi)


def free_freqs(T, L, EI, m, fmax):
    n = np.arange(1, 200)
    f = n / (2 * L) * np.sqrt(T / m) * np.sqrt(1 + n * n * np.pi ** 2 * EI / (T * L * L))
    return f[f < fmax * 1.1]


def stay_freqs(T, L, EA, m, theta, fmax):
    q = m * G * np.cos(np.radians(theta))
    sag = q * L * L / (8 * T)
    Le = L * (1 + 8 * (sag / L) ** 2)
    lam2 = (q * L / T) ** 2 * L * EA / (T * Le)
    base = np.sqrt(T / m) / L
    out = []
    for k in range(1, 40):
        if k % 2 == 0:
            out.append(k * np.pi * base / (2 * np.pi))
        else:
            h = lambda W: np.sin(W / 2) - np.cos(W / 2) * (W / 2 - (4 / lam2) * (W / 2) ** 3)
            lo, hi = k * np.pi, (k + 2) * np.pi
            out.append(brentq(h, lo, hi) * base / (2 * np.pi))
    f = np.array(sorted(out))
    return f[f < fmax * 1.1]


def misfit(obs, model):
    if len(model) == 0:
        return 1e9
    d = np.abs(np.log(obs[:, None] / model[None, :])).min(axis=1)
    lo, hi = obs.min() * 0.97, obs.max() * 1.03
    empty = sum(1 for v in model if lo < v < hi and np.abs(np.log(obs / v)).min() > 0.02)
    return float(np.sum(np.minimum(d, 0.02) ** 2) + 0.35 * 0.02 ** 2 * empty)


def matched(obs, model):
    out = 0.0
    for v in obs:
        j = np.argmin(np.abs(model - v))
        r = np.log(v / model[j])
        out += r * r if abs(r) < 0.01 else 1e-4
    return out


def main(folder, dest):
    folder = Path(folder)
    sheet = json.loads((folder / "member_card.json").read_text())
    fs = sheet["sample_rate_hz"]
    obs = []
    for tag in ("before", "after"):
        with open(folder / sheet["records"][tag]) as fh:
            acc = np.array([float(r["acc_mg"]) for r in csv.DictReader(fh)])
        obs.append(spectrum_peaks(acc, fs))
    m, L = sheet["mass_per_metre_kg"], sheet["length_between_working_points_m"]
    E = sheet["steel_modulus_mpa"] * 1e6
    fmax = max(o.max() for o in obs)
    member = sheet["member"]
    forces = []
    if "diameter_mm" in member:
        EI = E * np.pi * (member["diameter_mm"] / 1000) ** 4 / 64
        best = None
        for held, fn in ((True, held_freqs), (False, free_freqs)):
            for Lv in np.linspace(0.8 * L, L, 21):
                Ts, total = [], 0.0
                for o in obs:
                    grid = EI * (np.geomspace(2.0, 100.0, 300) / Lv) ** 2
                    cs = [misfit(o, fn(T, Lv, EI, m, fmax)) for T in grid]
                    Ts.append(grid[int(np.argmin(cs))])
                    total += min(cs)
                if best is None or total < best[0]:
                    best = (total, held, fn, Lv, Ts)
        _, held, fn, Lv, Ts = best
        res = minimize(lambda z: sum(matched(o, fn(np.exp(z[i]), z[2], EI, m, fmax)) for i, o in enumerate(obs)),
                       x0=[np.log(Ts[0]), np.log(Ts[1]), Lv], method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 3000})
        forces = [np.exp(res.x[0]), np.exp(res.x[1])]
    else:
        EA = sheet["metallic_area_mm2"] * 1e-6 * E
        th = sheet["chord_inclination_deg"]
        for o in obs:
            T0 = 4 * m * L * L * np.median(np.diff(np.sort(o))) ** 2
            grid = np.geomspace(0.15 * T0, 6 * T0, 300)
            T = grid[int(np.argmin([misfit(o, stay_freqs(t, L, EA, m, th, fmax)) for t in grid]))]
            res = minimize(lambda z: matched(o, stay_freqs(np.exp(z[0]), L, EA, m, th, fmax)), x0=[np.log(T)],
                           method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-14})
            forces.append(np.exp(res.x[0]))
    Path(dest).write_text(json.dumps({"axial_before_kn": forces[0] / 1e3, "axial_after_kn": forces[1] / 1e3}))


if __name__ == "__main__":
    main(*sys.argv[1:3])
