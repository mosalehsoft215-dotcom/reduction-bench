"""Reduce one pellet folder to the Fe3+ share of its iron and the Fe2+ (olivine) hyperfine parameters at 295 K.

    python3 fe_sites.py <pellet folder> <answer json>

1. Velocity scale and source line from the alpha-iron foil run: the six foil lines at their known positions with
   strengths 3:4:1:1:4:3 and natural width, the foil's effective thickness from its 10 um and f = 0.80, the full
   transmission integral; velocity = cubic in channel, source a Lorentzian of fitted width.
2. Both pellet runs together, each as the full transmission integral (the pellet is not thin): two doublets whose
   effective thicknesses are sigma0 x n57 (from the sheet) x site share x the site's Debye recoil-free fraction at
   that run's temperature, the run's resonant fraction from runs.csv. Fitted: the Fe3+ share, one Debye temperature
   per site, per-run hyperfine parameters and baselines, shared absorber widths; doublets symmetric (magic angle).
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

SIGMA0 = 2.56e-18
ABUND = 0.02119
GAMMA0 = 0.097
ER_OVER_K = 1.956e-3 / 8.617333262e-5          # recoil energy over Boltzmann constant, K
NA = 6.02214076e23
M_FE = 55.845
FOIL_POS = np.array([-5.3123, -3.0760, -0.8397, 0.8397, 3.0760, 5.3123])
FOIL_REL = np.array([3, 4, 1, 1, 4, 3]) / 16.0
FOIL_T = SIGMA0 * (7.874 * 10e-4 / M_FE * NA * ABUND) * 0.80
E = np.arange(-16.0, 16.0, 0.004)
_K = np.fft.ifftshift(np.arange(-len(E) // 2, len(E) // 2) * (E[1] - E[0]))


def recoilless(temp, theta):
    y = theta / temp
    integ = quad(lambda x: x / math.expm1(x) if x > 0 else 1.0, 0, y)[0]
    return math.exp(-(6 * ER_OVER_K / theta) * (0.25 + integ / y ** 2))


def transmitted(v, lines, gs):
    T = np.zeros_like(E)
    for pos, t, ga in lines:
        T += t * (GAMMA0 / 2) * (ga / 2) / ((E - pos) ** 2 + (ga / 2) ** 2)
    src = (gs / 2) / math.pi / (_K ** 2 + (gs / 2) ** 2) * (E[1] - E[0])
    conv = np.fft.irfft(np.fft.rfft(1 - np.exp(-T)) * np.fft.rfft(src), n=len(E))
    return 1 - np.interp(v, E, conv)


def read(folder, dead=True):
    folder = Path(folder)
    runs = {}
    with open(folder / "runs.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            data = np.loadtxt(folder / "spectra" / row["run"], delimiter=",", skiprows=1)
            c = data[:, 1].astype(float)
            if dead:                    # analyser dead time from live vs real time (non-paralyzable)
                real, live = float(row["real_time_s"]), float(row["live_time_s"])
                tau = (real - live) / c.sum()
                dwell = real / len(c)
                c = c / (1 - (c / dwell) * tau)
            runs[row["run"]] = dict(row, counts=c, T=float(row["temperature_K"]), r=float(row["resonant_fraction"]))
    sheet = json.loads((folder / "sheet.json").read_text())
    area = math.pi * (sheet["pellet_diameter_mm"] / 20.0) ** 2
    fe_mg_cm2 = sheet["pellet_mass_mg"] * sheet["iron_mass_percent"] / 100.0 / area
    n57 = fe_mg_cm2 * 1e-3 / M_FE * NA * ABUND
    return runs, n57


def velocity(a, u):
    return a[0] + a[1] * u + a[2] * u ** 2 + a[3] * u ** 3


def foil_fit(run, u, cubic=True):
    y = run["counts"]

    def model(p):
        v = velocity(p[:4], u)
        lines = [(pos, FOIL_T * s, GAMMA0) for pos, s in zip(FOIL_POS, FOIL_REL)]
        return p[5] * (1 - run["r"] + run["r"] * transmitted(v, lines, p[4]))

    p0 = [0.0, 8.0, 0.0, 0.0, 0.12, float(np.median(y[:20]))]
    def resid(q):                      # wrong-route switch: cubic=False keeps the scale linear in channel
        p = q if cubic else np.concatenate([q[:2], [0.0, 0.0], q[2:]])
        return (model(p) - y) / np.sqrt(y)
    q0 = p0 if cubic else p0[:2] + p0[4:]
    xs = [0.05, 0.1, 0.05, 0.05, 0.01, p0[5] * 0.01]
    fit = least_squares(resid, q0, x_scale=xs if cubic else xs[:2] + xs[4:])
    p = fit.x if cubic else np.concatenate([fit.x[:2], [0.0, 0.0], fit.x[2:]])
    return p[:4], p[4]


def reduce(folder, saturation=True, debye=True, use_temps=(295, 80), cubic=True, dead=True, stage_temp=True):
    """Keyword switches exist for the wrong-route record; the defaults are the reduction."""
    runs, n57 = read(folder, dead)
    u = (np.arange(512) - 255.5) / 255.5
    foil = next(r for r in runs.values() if "foil" in r["absorber"])
    a, gs = foil_fit(foil, u, cubic)
    v = velocity(a, u)
    pel = sorted((r for r in runs.values() if r["absorber"] == "pellet" and int(r["T"]) in use_temps),
                 key=lambda r: -r["T"])

    tied = stage_temp and len(pel) == 2
    # shared: x3, th2, th3, ga2, ga3 (+ chemical shifts c2, c3 and the pellet's cold temperature when tied);
    # per run: is2, qs2, is3, qs3, base (free shifts) or qs2, qs3, base (shifts from the second-order Doppler term)
    ns, nr = (8, 3) if tied else (5, 5)

    def unpack(p):
        return p[:ns], [p[ns + nr * i: ns + nr * (i + 1)] for i in range(len(pel))]

    def temperature(shared, run):
        return run["T"] if not tied or run["T"] > 200 else shared[7]

    def model_run(shared, own, run):
        x3, th2, th3, ga2, ga3 = shared[:5]
        T = temperature(shared, run)
        if tied:
            qs2, qs3, base = own
            is2, is3 = shared[5] + doppler(T, th2), shared[6] + doppler(T, th3)
        else:
            is2, qs2, is3, qs3, base = own
        if debye:
            f2, f3 = recoilless(T, th2), recoilless(T, th3)
        else:                                   # wrong route: one recoil-free fraction for both sites
            f2 = f3 = recoilless(T, th2)
        t2 = SIGMA0 * n57 * (1 - x3) * f2
        t3 = SIGMA0 * n57 * x3 * f3
        lines = [(is2 - qs2 / 2, t2 / 2, ga2), (is2 + qs2 / 2, t2 / 2, ga2),
                 (is3 - qs3 / 2, t3 / 2, ga3), (is3 + qs3 / 2, t3 / 2, ga3)]
        if saturation:
            frac = transmitted(v, lines, gs)
        else:                                   # wrong route: thin-absorber sum of Lorentzians
            frac = np.ones_like(v)
            for pos, t, ga in lines:            # each dip keeps the thin-limit area t pi G0/2, width ga + gs
                w = ga + gs
                frac -= t * (GAMMA0 / 2) * (w / 2) / ((v - pos) ** 2 + (w / 2) ** 2)
        return base * (1 - run["r"] + run["r"] * frac)

    def resid(p):
        shared, owns = unpack(p)
        return np.concatenate([(model_run(shared, own, run) - run["counts"]) / np.sqrt(run["counts"])
                               for own, run in zip(owns, pel)])

    best = None
    for x0 in (0.3, 0.7):
        for ths in ((320.0, 320.0), (300.0, 420.0)):
            p0 = [x0, ths[0], ths[1], 0.29, 0.33]
            lo, hi = [0.0, 120.0, 120.0, 0.1, 0.1], [1.0, 900.0, 900.0, 0.8, 0.8]
            scale = [0.02, 10.0, 10.0, 0.01, 0.01]
            if tied:
                p0 += [1.22 - doppler(295.0, ths[0]), 0.43 - doppler(295.0, ths[1]), 90.0]
                lo += [-np.inf, -np.inf, 40.0]
                hi += [np.inf, np.inf, 200.0]
                scale += [0.01, 0.01, 5.0]
            for run in pel:
                hot = run["T"] > 200
                own = [2.96 if hot else 2.99, 1.53 if hot else 1.55] if tied else \
                    [1.22 if hot else 1.33, 2.96 if hot else 2.99, 0.43 if hot else 0.54, 1.53 if hot else 1.55]
                p0 += own + [float(np.median(run["counts"][:20]))]
                lo += [-np.inf] * nr
                hi += [np.inf] * nr
                scale += [0.01] * (nr - 1) + [1000.0]
            fit = least_squares(resid, p0, bounds=(lo, hi), x_scale=scale)
            if best is None or fit.cost < best.cost:
                best = fit
    shared, owns = unpack(best.x)
    hot = next((i for i, run in enumerate(pel) if run["T"] > 200), 0)
    if tied:
        is2_hot, qs2_hot = shared[5] + doppler(295.0, shared[1]), owns[hot][0]
    else:
        is2_hot, qs2_hot = owns[hot][0], owns[hot][1]
    res = {"fe3_fraction": float(shared[0]), "is_fe2_mm_s": float(is2_hot), "qs_fe2_mm_s": float(qs2_hot)}
    return res, {"theta2": shared[1], "theta3": shared[2], "gs": gs, "a": list(a),
                 "cold_K": float(shared[7]) if tied else None}


def doppler(temp, theta):
    """Second-order Doppler shift of a Debye lattice for 57Fe, mm/s."""
    k, m, c = 1.380649e-23, 56.9354 * 1.66053907e-27, 299792458.0
    y = theta / temp
    integ = quad(lambda x: x ** 3 / math.expm1(x) if x > 0 else 0.0, 0, y)[0]
    return -1e3 * (9 * k * theta / (16 * m * c) + 9 * k * temp / (2 * m * c) * integ / y ** 3)


def main():
    res, _ = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
