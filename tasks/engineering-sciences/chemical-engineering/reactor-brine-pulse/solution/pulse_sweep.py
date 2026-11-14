"""Swept volume and mixing stages of one baffled reactor from its brine-pulse record.

usage: python3 pulse_sweep.py <vessel folder> <result file>

Everything is written against the feed volume u passed since the dose (trapezoid of the logged feed rate), because
the vessel's internal flows follow the feed. Model of the outlet reading:
    k(u) = k0 + dk * F(u + w) + a * E(u)
E is the residence distribution of a line plug Vp, n equal stirred cells holding Vm, and a slack zone Vd trading
with the last cell at r times the feed; F is its running integral, which the feed's own conductivity follows when it
changed w m3 before the dose. For each n in 1..9 the nonlinear set (Vp, Vm, Vd, r, w) is fitted by bounded least
squares with k0, dk and a solved linearly at every step; the best n is refined from more starts. The answer is the
first moment of E and its square over the variance, both from the cumulants of the fitted transfer function.
"""
import json
import sys

import numpy as np
from scipy.optimize import least_squares

GRID = 1 << 13


def load(folder):
    sheet = json.load(open(f"{folder}/dosing_sheet.json"))
    raw = np.genfromtxt(f"{folder}/outlet_log.csv", delimiter=",", names=True)
    return sheet, raw["elapsed_s"], raw["outlet_mS_cm"], raw["feed_m3_h"]


def passed_volume(t, q):
    return np.concatenate(([0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t)))) / 3600.0


def tf(s, vp, vm, n, vd, r):
    vc = vm / n
    return np.exp(-vp * s) / (1 + vc * s) ** (n - 1) / (vc * s + 1 + r * vd * s / (r + vd * s))


def cumulants(vp, vm, n, vd, r):
    h = 1e-3 / (vp + vm + vd)
    k = [np.log(tf(x, vp, vm, n, vd, r).real) for x in (-h, 0.0, h)]
    return -(k[2] - k[0]) / (2 * h), (k[2] - 2 * k[1] + k[0]) / h ** 2


class PulseFit:
    def __init__(self, u, k, scale):
        self.u, self.k, self.scale = u, k, scale
        span = 40 * scale
        self.du = span / GRID
        self.iw = 2j * np.pi * np.fft.rfftfreq(GRID, self.du)
        self.grid = np.arange(GRID) * self.du

    def columns(self, x, n):
        vp, vm, vd, r = np.exp(x[:4])
        e = np.fft.irfft(tf(self.iw, vp, vm, n, vd, r), GRID) / self.du
        f = np.concatenate(([0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * self.du)))
        ce = np.interp(self.u, self.grid, e, left=0.0, right=0.0)
        cf = np.interp(self.u + x[4] * self.scale, self.grid, f, left=0.0, right=1.0)
        return np.column_stack((np.ones_like(self.u), cf, ce))

    def misfit(self, x, n):
        cols = self.columns(x, n)
        coef = np.linalg.lstsq(cols, self.k, rcond=None)[0]
        return cols @ coef - self.k

    def solve(self, n, start):
        lo = np.log([1e-3, 0.05, 1e-3, 0.01]) + np.log([self.scale, self.scale, self.scale, 1.0])
        hi = np.log([0.8, 3.0, 1.0, 10.0]) + np.log([self.scale, self.scale, self.scale, 1.0])
        x0 = np.r_[np.log(np.array(start[:4]) * [self.scale, self.scale, self.scale, 1.0]), start[4]]
        out = least_squares(self.misfit, x0, args=(n,), bounds=(np.r_[lo, 0.0], np.r_[hi, 3.0]),
                            x_scale="jac", max_nfev=400)
        return out.cost, out.x


def reduce(folder):
    sheet, t, k, q = load(folder)
    v = passed_volume(t, q)
    u = v - np.interp(sheet["dose"]["dosed_at_s"], t, v)
    before = k[u < 0]
    excess = np.clip(k - np.median(before), 0.0, None) * (u >= 0)
    run = np.cumsum(excess)
    scale = 1.2 * float(np.interp(0.5 * run[-1], run, u))
    fit = PulseFit(u, k, max(scale, 0.5))
    first = (0.05, 0.75, 0.2, 0.3, 0.3)
    ranked = sorted((fit.solve(n, first) + (n,) for n in range(1, 10)), key=lambda z: z[0])
    best = ranked[0]
    for _, _, n in ranked[:3]:
        for start in ((0.05, 0.75, 0.2, 0.3, 1.0), (0.05, 0.6, 0.05, 1.0, 0.3), (0.1, 0.5, 0.4, 0.1, 0.6)):
            trial = fit.solve(n, start) + (n,)
            if trial[0] < best[0]:
                best = trial
    _, x, n = best
    vp, vm, vd, r = np.exp(x[:4])
    mean, var = cumulants(vp, vm, n, vd, r)
    return float(mean), float(mean ** 2 / var)


def main(folder, target):
    swept, stages = reduce(folder)
    with open(target, "w") as fh:
        json.dump({"swept_volume_m3": round(swept, 4), "mixing_stages": round(stages, 4)}, fh)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
