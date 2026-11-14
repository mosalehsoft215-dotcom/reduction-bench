"""Second, independently written reduction of a zero-field muSR scan (no code shared with solution/zf_reduce.py).

    python3 second_scan.py <scan folder> <result file>

Differences from the reference: the time axis is calibrated by a Poisson maximum-likelihood lifetime fit of the summed
counters; each run is fitted by variable projection (non-linear parameters on grids, amplitudes by weighted linear
least squares); folded frequencies are unfolded by quadratic extrapolation from the warmer runs; the Neel temperature is
the linear interpolation of the ordered fraction through one half.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

TAU_MU = 2.1969811
MHZ_PER_MT = 0.135538817


def read_scan(folder):
    sheet = json.loads((folder / "scan_sheet.json").read_text())
    note = (folder / "scan_log.txt").read_text()
    moved = re.search(r"Run\s+(\d+)\s*:[^\n]*?(\d+(?:\.\d+)?)\s*K\s+(below|above)\s+the sample", note)
    rows = [ln.split(",") for ln in (folder / "runs.csv").read_text().split("\n")[1:] if ln.strip()]
    out = []
    for run, temp, heat, _ in rows:
        if float(heat) == 0.0:
            continue
        T = float(temp)
        if moved and int(run) >= int(moved.group(1)):
            T += float(moved.group(2)) if moved.group(3) == "below" else -float(moved.group(2))
        out.append((T, np.genfromtxt(folder / "histograms" / f"run_{run}.csv", delimiter=",", names=True)))
    return sheet, sorted(out, key=lambda x: x[0])


def lifetime_scale(hists, alpha):
    t = hists[0]["time_us"]
    total = sum(h["forward"] + h["backward"] / alpha for h in hists)
    bg = total[t < 0].mean()
    m = (t > 0.05) & (t < 9.0)
    tt, y = t[m], total[m]

    def nll(tau):
        shape = np.exp(-tt / tau)
        amp = (y - bg).clip(0).sum() / shape.sum()
        mu = amp * shape + bg
        return float(np.sum(mu - y * np.log(mu)))
    tau = minimize_scalar(nll, bounds=(1.9, 2.5), method="bounded", options={"xatol": 1e-7}).x
    return tau / TAU_MU


def kubo(s, t):
    u = (s * t) ** 2
    return 1 / 3 + 2 / 3 * (1 - u) * np.exp(-u / 2)


def run_asym(h, alpha):
    t, F, B = h["time_us"], h["forward"].astype(float), h["backward"].astype(float)
    F0, B0 = F - F[t < 0].mean(), B - B[t < 0].mean()
    Bn = B0 / alpha                                    # backward counts brought to the forward counter's scale
    m = (t > 0) & (t < 8.0) & (F0 + Bn > 0)
    a = (F0[m] - Bn[m]) / (F0[m] + Bn[m])
    sd = 2 / alpha * np.sqrt(B0[m] ** 2 * F[m] + F0[m] ** 2 * B[m]) / (F0[m] + Bn[m]) ** 2
    return t[m], a, 1 / np.maximum(sd, 1e-4)


def linfit(cols, a, w):
    X = np.stack(cols, 1) * w[:, None]
    coef, *_ = np.linalg.lstsq(X, a * w, rcond=None)
    r = a * w - X @ coef
    return coef, float(r @ r)


def best_osc(t, a, w, base_cols, fmax, d, exclude=None):
    """Grid over frequency and damping; cos and sin amplitudes linear. Returns (chi2, f, amplitude)."""
    best = (np.inf, 0.0, 0.0)
    for lam in (0.5, 1.5, 3.0, 6.0, 12.0):
        env = np.exp(-lam * t)
        for f in np.arange(1.0, fmax, 0.05):
            if exclude is not None and abs(f - exclude) < 1.0:
                continue
            att = np.sinc(f * d)
            c, chi = linfit(base_cols + [env * np.cos(2 * np.pi * f * t) * att, env * np.sin(2 * np.pi * f * t) * att], a, w)
            if chi < best[0]:
                best = (chi, f, float(np.hypot(c[-2], c[-1])), lam)
    chi, f0, amp, lam = best

    def refine(f):
        env = np.exp(-lam * t)
        att = np.sinc(f * d)
        return linfit(base_cols + [env * np.cos(2 * np.pi * f * t) * att, env * np.sin(2 * np.pi * f * t) * att], a, w)[1]
    f1 = minimize_scalar(refine, bounds=(max(0.5, f0 - 0.06), min(fmax, f0 + 0.06)), method="bounded").x
    return refine(f1), f1, amp


def reduce(folder):
    sheet, runs = read_scan(folder)
    alpha, d = sheet["alpha_forward_backward"], sheet["bin_width_us"]
    fmax = 0.5 / d * 0.995
    scale = lifetime_scale([h for _, h in runs], alpha)
    # paramagnetic width from the two warmest runs (grid in sigma, amplitudes linear)
    warm = [run_asym(h, alpha) for _, h in runs[-2:]]
    sig = min(np.arange(0.05, 1.5, 0.002), key=lambda s: sum(linfit([kubo(s, t), np.ones_like(t)], a, w)[1] for t, a, w in warm))
    A0 = np.mean([linfit([kubo(sig, t), np.ones_like(t)], a, w)[0][0] for t, a, w in warm])
    fracs, freqs = [], []
    for T, h in runs:
        t, a, w = run_asym(h, alpha)
        base = [kubo(sig, t), np.exp(-0.05 * t), np.ones_like(t)]
        chi_fast = min(linfit(base + [np.exp(-lf * t)], a, w)[1] for lf in (3, 6, 10, 15, 25, 40))
        chi1, f1, amp1 = best_osc(t, a, w, base, fmax, d)
        coef, _ = linfit(base, a, w)
        frac = 1 - linfit(base + [np.exp(-15 * t)] if chi_fast < chi1 else base, a, w)[0][0] / A0
        if chi_fast - chi1 > 40:
            env1 = np.exp(-1.5 * t) * np.sinc(f1 * d)
            base2 = base + [env1 * np.cos(2 * np.pi * f1 * t), env1 * np.sin(2 * np.pi * f1 * t)]
            chi2, f2, amp2 = best_osc(t, a, w, base2, fmax, d, exclude=f1)
            if chi1 - chi2 > 60 and max(f1, f2) < 0.97 * fmax:
                best = None
                for l1 in np.geomspace(0.3, 15, 16):        # each site its own damping
                    for l2 in np.geomspace(0.3, 15, 16):
                        cols = [np.exp(-l * t) * fn(2 * np.pi * f * t) * np.sinc(f * d)
                                for fn in (np.cos, np.sin) for f, l in ((f1, l1), (f2, l2))]
                        c, chi = linfit(base + cols, a, w)
                        if best is None or chi < best[0]:
                            best = (chi, c)
                c = best[1]
                a1, a2 = np.hypot(c[3], c[5]), np.hypot(c[4], c[6])
                f1 = f1 if a1 >= a2 else f2
            freqs.append((T, f1))
        fracs.append((T, min(max(frac, 0.0), 1.2)))
    # one fold at most: try every run as the first one above half the sampling rate, keep the smoothest curve
    fN = 0.5 / d
    seq = sorted(freqs, key=lambda x: -x[0])
    Ts = np.array([x[0] for x in seq])
    fa = np.array([x[1] for x in seq])

    def roughness(g):
        if len(g) < 3:
            return 0.0
        slope = np.diff(g) / np.diff(Ts)
        return float(np.sum(np.abs(np.diff(slope))))
    # a fold can only sit where the precession has come close to half the sampling rate
    options = [np.concatenate([fa[:k], 2 * fN - fa[k:]]) for k in range(len(fa) + 1)
               if k == len(fa) or max(fa[max(k - 1, 0)], fa[k]) > 0.8 * fN]
    best = min(options, key=roughness)
    unf = list(zip(Ts, best))
    f_cold = unf[-1][1] * scale
    Tf = np.array([x[0] for x in fracs])
    fr = np.array([x[1] for x in fracs])
    from scipy.special import erfc
    from scipy.optimize import curve_fit
    step = lambda T, top, tn, width: 0.5 * top * erfc((T - tn) / width)     # error-function step, midpoint tn
    i = int(np.argmin(np.abs(fr - 0.5 * np.median(fr[:3]))))
    tn = curve_fit(step, Tf, fr, p0=[np.median(fr[:3]), Tf[i], 0.5], maxfev=20000)[0][1]
    return {"local_field_mt": float(f_cold / MHZ_PER_MT), "neel_temperature_k": float(tn)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(Path(sys.argv[1]))) + "\n")
