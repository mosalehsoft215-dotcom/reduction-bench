"""Independent reduction of one discharge, written apart from solution/column_temp.py.

End window from the two lamp checks by a local running median, ln W interpolated in running time; gas temperature
from the closed fill and the running pressure; scale offset from line centroids; flat-topped lines by fixed-width
linear least squares on the pixels below the converter's top.  Line strengths by direct integration: the response-corrected spectrum is summed over +-0.12 nm around each line
after a straight continuum through the medians of two side bands; lines closer than 0.4 nm are integrated as one
group and modelled as the sum of their members.  Each line's centre optical depth is k0 L with
k0 = (lambda0^2 / 8 pi) (g_u / g_l) A n_l phi(nu0), phi(nu0) = lambda0 / (sqrt(pi) v_th), and the column's escape
factor is evaluated by 80-point Gauss-Hermite quadrature of int (1 - exp(-tau exp(-x^2))) dx.  The temperature and
a common density come from a log-spaced grid search on the weighted squared log residuals (the overall scale solved
in closed form), then Nelder-Mead over the temperature and the four 4s densities from the best grid node.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

KB = 8.617333262e-5


def groups(lines):
    order = np.argsort(lines[:, 0])
    out, cur = [], [order[0]]
    for i in order[1:]:
        if lines[i, 0] - lines[cur[-1], 0] < 0.4:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    out.append(cur)
    return out


def lamp_ratio(path, rl, rr):
    """ln of the window transmission seen in one lamp check: the lamp over a 2856 K photon spectrum and the bench
    response, as a running median over 2 nm, relative to the 640-660 nm band (no global smoothing law)."""
    ll, lc = np.loadtxt(path, delimiter=",", skiprows=1, unpack=True)
    ratio = lc / (ll ** -4 / np.expm1(1.4387769e7 / (ll * 2856.0))) / np.interp(ll, rl, rr)
    ratio = ratio / np.median(ratio[(ll > 640) & (ll < 660)])
    centres = np.arange(405.0, 871.0, 1.0)
    return centres, np.log(np.array([np.median(ratio[np.abs(ll - c) <= 1.0]) for c in centres]))


def strengths(folder, info):
    lam, counts = np.loadtxt(folder / "spectrum.csv", delimiter=",", skiprows=1, unpack=True)
    rl, rr = np.loadtxt(folder / "response.csv", delimiter=",", skiprows=1, unpack=True)
    flux = counts / np.interp(lam, rl, rr)
    # window at the spectrum's running time: ln W linear in running time between the first and the last check
    checks = sorted(info["lamp_checks"].items(), key=lambda kv: kv[1])
    c0, l0 = lamp_ratio(folder / checks[0][0], rl, rr)
    _, l1 = lamp_ratio(folder / checks[-1][0], rl, rr)
    span = checks[-1][1] - checks[0][1]
    frac = (info["spectrum_running_time_h"] - checks[0][1]) / span if span > 0 else 0.0
    flux = flux / np.exp(np.interp(lam, c0, l0 + frac * (l1 - l0)))
    lines = np.loadtxt(folder / "ar_lines.csv", delimiter=",", skiprows=1)
    flat = counts >= info["frames"] * 65535.0 - 0.5            # pixels at the converter's top
    base = np.median(flux)
    # scale offset: centroids of the six strongest lines whose +-0.3 nm hold no flat pixel, median offset
    offs = []
    for l0_ in lines[np.argsort(-lines[:, 2]), 0]:
        win = np.abs(lam - l0_) < 0.3
        if flat[win].any() or len(offs) == 6:
            continue
        pk = lam[win][np.argmax(flux[win])]
        core = np.abs(lam - pk) < 0.08
        wgt = np.maximum(flux[core] - base, 0)
        offs.append(float(np.sum(wgt * lam[core]) / wgt.sum()) - l0_)
    shift = float(np.median(offs))
    # instrument width from the same clean lines: area over peak height / sqrt(2 pi)
    sig = []
    for l0_ in lines[:, 0]:
        win = np.abs(lam - l0_ - shift) < 0.15
        if flat[np.abs(lam - l0_ - shift) < 0.3].any() or (np.abs(lines[:, 0] - l0_) < 0.6).sum() > 1:
            continue
        y = flux[win] - np.median(flux[(np.abs(lam - l0_ - shift) > 0.2) & (np.abs(lam - l0_ - shift) < 0.3)])
        sig.append((np.trapezoid(y, lam[win]) / (y.max() * math.sqrt(2 * math.pi)), y.max()))
    width = sum(s_ * h for s_, h in sig) / sum(h for _, h in sig)
    measured = []
    for g in groups(lines):
        lo, hi = lines[g, 0].min() + shift - 0.12, lines[g, 0].max() + shift + 0.12
        near = (lam > lo - 0.25) & (lam < hi + 0.25)
        if flat[near].any():
            # flat tops: amplitudes of fixed-width Gaussians plus a straight continuum by linear least squares on the
            # pixels below the top, the group's position stepped over +-0.02 nm
            keep = near & ~flat
            best = None
            for d in np.arange(-0.02, 0.0201, 0.001):
                cols = [np.exp(-0.5 * ((lam[keep] - c - shift - d) / width) ** 2) for c in lines[g, 0]]
                X = np.stack(cols + [np.ones(keep.sum()), lam[keep] - lam[keep].mean()], 1)
                coef, res, *_ = np.linalg.lstsq(X, flux[keep], rcond=None)
                r = float(np.sum((X @ coef - flux[keep]) ** 2))
                if best is None or r < best[0]:
                    best = (r, coef)
            measured.append((g, float(np.sum(best[1][:len(g)]) * width * math.sqrt(2 * math.pi))))
            continue
        left = flux[(lam > lo - 0.25) & (lam < lo - 0.13)]
        right = flux[(lam > hi + 0.13) & (lam < hi + 0.25)]
        inside = (lam >= lo) & (lam <= hi)
        basel = np.interp(lam[inside], [lo - 0.19, hi + 0.19], [np.median(left), np.median(right)])
        measured.append((g, float(np.trapezoid(flux[inside] - basel, lam[inside]))))
    return lines, measured


LOWER = (11.548, 11.624, 11.723, 11.828)


GH_X, GH_W = np.polynomial.hermite.hermgauss(80)


_LG = np.linspace(-8.0, 6.0, 2801)                 # log10 tau nodes
_XQ = np.linspace(-7.5, 7.5, 3001)
_ESC = np.trapezoid(1 - np.exp(-(10.0 ** _LG)[:, None] * np.exp(-_XQ ** 2)[None, :]), _XQ, axis=1) / (10.0 ** _LG * math.sqrt(math.pi))


def column_escape(tau):
    """(1/(tau sqrt(pi))) int (1 - exp(-tau e^{-x^2})) dx: tabulated once by the trapezoid rule on 2801 log-spaced
    depths, read by linear interpolation in log tau; small tau by the series 1 - tau/sqrt(8) + tau^2/(6 sqrt(3))."""
    tau = np.asarray(tau, dtype=float)
    small = tau < 1e-3
    return np.where(small, 1 - tau / math.sqrt(8) + tau ** 2 / (6 * math.sqrt(3)),
                    np.interp(np.log10(np.maximum(tau, 1e-300)), _LG, _ESC))


def reduce(folder):
    info = json.loads((Path(folder) / "discharge.json").read_text())
    lines, measured = strengths(Path(folder), info)
    lam, eu, ga, el = lines[:, 0], lines[:, 1], lines[:, 2] / 1e7, lines[:, 3]
    gl = lines[:, 4]
    t_gas = 293.15 * info["running_pressure_Pa"] / info["fill_pressure_Pa"]       # closed fill, uniform gas
    vth = math.sqrt(2 * 1.380649e-23 * t_gas / (39.948 * 1.66053907e-27))
    lam_m = lam * 1e-9
    k_per_n = lam_m ** 2 / (8 * math.pi) * (ga * 1e7 / gl) * lam_m / (math.sqrt(math.pi) * vth) * info["column_length_m"]
    which = np.array([int(np.argmin(np.abs(np.array(LOWER) - e))) for e in el])
    obs = np.array([m for _, m in measured])
    wts = np.sqrt(np.maximum(obs, 1e-12))

    def predicted(t, log_s):
        thin = ga * np.exp(-eu / (KB * t))
        tau = 10.0 ** np.asarray(log_s)[which] * k_per_n
        esc = column_escape(tau)
        per = thin * esc
        return np.array([per[g].sum() for g, _ in measured])

    def cost(p):
        t, log_s = p[0], p[1:]
        if not 2500 < t < 30000 or np.any(np.asarray(log_s) > 19) or np.any(np.asarray(log_s) < 10):
            return 1e30
        m = np.log(np.maximum(predicted(t, log_s), 1e-300))
        r = np.log(np.maximum(obs, 1e-300)) - m
        off = np.sum(wts * r) / wts.sum()                 # overall scale in closed form
        return float(np.sum(wts * (r - off) ** 2))

    grid = [(cost([t, ls, ls, ls, ls]), t, ls) for t in np.geomspace(3000, 20000, 70) for ls in np.linspace(12.0, 17.5, 23)]
    _, t0, s0 = min(grid)
    best = None
    for split in (0.0, -0.5, -1.0):
        start = [t0, s0, s0 + split, s0 + split, s0 + split]
        fit = minimize(cost, start, method="Nelder-Mead",
                       options={"xatol": 1e-5, "fatol": 1e-14, "maxiter": 20000, "maxfev": 20000, "adaptive": True})
        if best is None or fit.fun < best.fun:
            best = fit
    return {"excitation_temperature_K": float(best.x[0]), "metastable_density_per_m3": float(10.0 ** best.x[1])}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + "\n", encoding="utf-8")
