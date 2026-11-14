"""Reference reduction of one zero-field muSR temperature scan.

    python3 zf_reduce.py <scan folder> <result file>

1. Thermometry: the scan log may move the temperature column to another sensor from a given run on; its stated offset
   is added back to every later run.
2. Asymmetry per run: flat backgrounds from the pre-t0 bins, a = (F - B / alpha) / (F + B / alpha), Poisson errors.
3. Paramagnetic reference: the two warmest runs fix the full sample asymmetry A0, the flat background Abg and the
   nuclear Kubo-Toyabe width.
4. Each run: A0 [fm (2/3 O + 1/3 exp(-lamL t)) + (1 - fm) KT] + Abg, with O fast-relaxing, one precessing site, or two
   sites, each precession damped and multiplied by the bin average sinc(f d); the model with the clearly lower
   chi-square wins.
5. Frequencies above half the sampling rate fold back. Runs are followed from warm to cold and each precession
   frequency is unfolded to the alias closest to the trend of the warmer runs.
6. local_field_mt: the dominant site (larger share of the precessing asymmetry) at the coldest run.
   neel_temperature_k: midpoint of a logistic fit to the magnetic volume fraction fm against temperature.
"""
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

GAMMA = 0.135538817                                   # MHz per mT
# steps of the reduction; authoring/evidence/wrong_routes.py switches them off one at a time
STEPS = {"log_offset": True, "tripped_runs": True, "lifetime_clock": True, "unfold": True, "dominant_site": True,
         "volume_fraction": True}


def load(folder):
    sheet = json.loads((folder / "scan_sheet.json").read_text())
    log = (folder / "scan_log.txt").read_text()
    # notes "Run <n>: ..." giving the column's offset from the sample, as the guide defines them
    m = re.search(r"Run (\d+):[^\n]*?([0-9]+(?:\.[0-9]+)?) K (below|above) the sample", log)
    first, shift = (None, 0.0)
    if m and STEPS["log_offset"]:
        first, shift = int(m.group(1)), float(m.group(2)) * (1.0 if m.group(3) == "below" else -1.0)
    runs = []
    for r in csv.DictReader(open(folder / "runs.csv")):
        run = int(r["run"])
        if STEPS["tripped_runs"] and float(r["heater_pct"]) <= 0.0:                    # controller tripped: the temperature readback is not the sample's
            continue
        T = float(r["temperature_k"]) + (shift if first is not None and run >= first else 0.0)
        h = np.loadtxt(folder / "histograms" / f"run_{run}.csv", delimiter=",", skiprows=1)
        runs.append((run, T, h))
    return sheet, runs


def asymmetry(h, alpha, t_max=8.0):
    t, F, B = h[:, 0], h[:, 1], h[:, 2]
    pre = t < 0
    bF, bB = F[pre].mean(), B[pre].mean()
    m = (t > 0) & (t <= t_max)
    f, b = F[m] - bF, B[m] - bB
    s = f + b / alpha                                  # N_B = alpha n exp(-t/tau) (1 - A)
    a = (f - b / alpha) / s
    err = 2 / alpha * np.sqrt(b ** 2 * F[m] + f ** 2 * B[m]) / s ** 2
    ok = s > 0
    return t[m][ok], a[ok], np.maximum(err[ok], 1e-4)


TAU = 2.1969811


def time_scale(raw, alpha, tmin=0.1, tmax=8.0):
    """Ratio of the labelled time axis to true time: F + B / alpha carries no asymmetry and decays with the muon lifetime."""
    num = den = 0.0
    for _, _, h in raw:
        t, F, B = h[:, 0], h[:, 1], h[:, 2]
        pre = t < 0
        s = (F - F[pre].mean()) + (B - B[pre].mean()) / alpha        # N_F + N_B / alpha = 2 n exp(-t/tau)
        m = (t > tmin) & (t < tmax) & (s > 0)
        w = s[m]                                          # Poisson weights for log counts
        slope = np.polyfit(t[m], np.log(s[m]), 1, w=np.sqrt(w))[0]
        num += -1 / slope * w.sum()
        den += w.sum()
    return (num / den) / TAU


def kt(sig, t):
    x = (sig * t) ** 2
    return 1 / 3 + 2 / 3 * (1 - x) * np.exp(-x / 2)


class Run:
    def __init__(self, T, t, a, e, d):
        self.T, self.t, self.a, self.e, self.d = T, t, a, e, d


def para_reference(runs):
    def res(q):
        A0, Abg, sig = q
        return np.concatenate([(A0 * kt(sig, r.t) + Abg - r.a) / r.e for r in runs])
    q = least_squares(res, [0.22, 0.02, 0.3], bounds=([0.05, -0.1, 0.01], [0.4, 0.2, 2.0])).x
    return q


def osc(t, f, phi, lam, d):
    return np.cos(2 * np.pi * f * t + phi) * np.exp(-lam * t) * np.sinc(f * d)


def fit_run(r, A0, Abg, sig, fN):
    t, a, e = r.t, r.a, r.e
    para = kt(sig, t)

    def model(O, fm, lamL):
        return A0 * (fm * (2 / 3 * O + 1 / 3 * np.exp(-lamL * t)) + (1 - fm) * para) + Abg

    out = {}
    # fast-relaxing (no coherent precession)
    r0 = least_squares(lambda q: (model(np.exp(-q[2] * t), q[0], q[1]) - a) / e, [0.5, 0.05, 10.0],
                       bounds=([0, 0, 0.5], [1, 2, 200]))
    out[0] = (2 * r0.cost, {"fm": r0.x[0]})
    # one precessing site: start from the periodogram peak
    res = a - model(np.exp(-r0.x[2] * t), r0.x[0], r0.x[1])          # what the non-precessing model leaves
    mask = t < 4.0
    freqs = np.linspace(1.0, 0.995 * fN, 1500)
    w = np.exp(-t[mask] / 1.5)
    spec = np.abs(np.exp(-2j * np.pi * freqs[:, None] * t[mask][None, :]) @ (res[mask] * w))
    f0 = freqs[np.argmax(spec)]
    best1 = None
    for fs in (f0, min(f0 * 0.98, 0.999 * fN), min(f0 * 1.02, 0.999 * fN)):
        for ph in (0.0, np.pi):
            r1 = least_squares(lambda q: (model(osc(t, q[3], q[4], q[2], r.d), q[0], q[1]) - a) / e,
                               [0.9, 0.05, 2.0, fs, ph], bounds=([0, 0, 0, 0.1, -2 * np.pi], [1, 2, 80, fN, 2 * np.pi]))
            if best1 is None or r1.cost < best1.cost:
                best1 = r1
    out[1] = (2 * best1.cost, {"fm": best1.x[0], "f": [best1.x[3]], "w": [1.0], "f1": best1.x[3]})
    # two sites: second peak of the residual periodogram
    fm1, lamL1, lam1, f1, ph1 = best1.x
    res2 = a - model(osc(t, f1, ph1, lam1, r.d), fm1, lamL1)
    spec2 = np.abs(np.exp(-2j * np.pi * freqs[:, None] * t[mask][None, :]) @ (res2[mask] * w))
    spec2[np.abs(freqs - f1) < 0.8] = 0
    f2 = freqs[np.argmax(spec2)]
    best2 = None
    for ww in (0.5, 0.7, 0.3):
        q0 = [fm1, lamL1, lam1, min(f1, 0.999 * fN), ph1, ww, lam1, min(f2, 0.999 * fN)]
        r2 = least_squares(lambda q: (model(q[5] * osc(t, q[3], q[4], q[2], r.d) + (1 - q[5]) * osc(t, q[7], q[4], q[6], r.d),
                                            q[0], q[1]) - a) / e, q0,
                           bounds=([0, 0, 0, 0.1, -2 * np.pi, 0.02, 0, 0.1], [1, 2, 80, fN, 2 * np.pi, 0.98, 80, fN]))
        if best2 is None or r2.cost < best2.cost:
            best2 = r2
    x = best2.x
    out[2] = (2 * best2.cost, {"fm": x[0], "f": [x[3], x[7]], "w": [x[5], 1 - x[5]], "f1": best1.x[3]})
    choice = 0
    if out[0][0] - out[1][0] > 40:
        choice = 1
        if out[1][0] - out[2][0] > 60 and min(x[5], 1 - x[5]) > 0.15 and abs(x[3] - x[7]) > 1.0 and max(x[3], x[7]) < 0.97 * fN:
            choice = 2
    return choice, out[choice][1]


def unfold(seq, fN):
    """seq: (T, apparent f) from warm to cold. Each frequency may be folded about half the sampling rate; the path
    through the aliases with the smoothest frequency-temperature slope is kept (dynamic programming)."""
    cands = [[f, 2 * fN - f, 2 * fN + f, 4 * fN - f] for _, f in seq]
    n = len(seq)
    if n == 1:
        return [seq[0][1]]
    T = [x[0] for x in seq]
    start = lambda j: 0.0 if j == 0 else 1e6                 # the warmest precession is slow, below half the rate
    cost = {(i, j): start(i) for i in range(4) for j in range(4)}
    back = []
    for k in range(2, n):
        new, bk = {}, {}
        for j in range(4):
            for m in range(4):
                best = None
                for i in range(4):
                    s1 = (cands[k - 1][j] - cands[k - 2][i]) / (T[k - 1] - T[k - 2])
                    s2 = (cands[k][m] - cands[k - 1][j]) / (T[k] - T[k - 1])
                    c = cost[(i, j)] + abs(s2 - s1)
                    if best is None or c < best[0]:
                        best = (c, i)
                new[(j, m)], bk[(j, m)] = best
        cost = new
        back.append(bk)
    (j, m) = min(cost, key=cost.get)
    path = [m, j]
    for bk in reversed(back):
        i = bk[(path[-1], path[-2])]
        path.append(i)
    path = path[::-1]
    return [cands[k][path[k]] for k in range(n)]


def reduce(folder):
    sheet, raw = load(folder)
    d, alpha = sheet["bin_width_us"], sheet["alpha_forward_backward"]
    fN = 0.5 / d
    scale = time_scale(raw, alpha) if STEPS["lifetime_clock"] else 1.0
    runs = sorted((Run(T, *asymmetry(h, alpha), d) for _, T, h in raw), key=lambda r: r.T)
    A0, Abg, sig = para_reference(runs[-2:])
    fits = [fit_run(r, A0, Abg, sig, fN) for r in runs]
    # dominant site per precessing run, followed from warm to cold
    pick = (lambda info: max(zip(info["w"], info["f"]))[1]) if STEPS["dominant_site"] else (lambda info: info["f1"])
    prec = [(r.T, pick(info)) for r, (c, info) in zip(runs, fits) if c > 0]
    prec = sorted(prec, key=lambda x: -x[0])
    true = unfold(prec, fN) if STEPS["unfold"] else [f for _, f in prec]
    f_base = true[-1] * scale                           # labelled-time frequency to true frequency
    Ts = np.array([r.T for r in runs])
    fm = np.array([info["fm"] for _, info in fits])
    q = least_squares(lambda p: p[0] / (1 + np.exp((Ts - p[1]) / p[2])) - fm,
                      [1.0, Ts[np.argmin(np.abs(fm - 0.5))], 0.5], bounds=([0.3, 0.5, 0.02], [1.2, 400, 10])).x
    tn = float(q[1])
    if not STEPS["volume_fraction"]:                      # transition placed where precession stops
        warmest = max(r.T for r, (c, _) in zip(runs, fits) if c > 0)
        tn = 0.5 * (warmest + min(r.T for r in runs if r.T > warmest))
    return {"local_field_mt": float(f_base / GAMMA), "neel_temperature_k": tn,
            "diagnostics": {"time_scale": scale, "A0": A0, "Abg": Abg, "kt_width": sig, "models": [c for c, _ in fits],
                            "precession_mhz": [round(x, 4) for x in true]}}


if __name__ == "__main__":
    out = reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + "\n")
