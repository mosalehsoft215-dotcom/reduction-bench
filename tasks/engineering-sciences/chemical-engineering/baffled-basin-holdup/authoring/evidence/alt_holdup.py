"""Second, separately written evaluation of one basin trial (time-domain kernels, closed-form moments).

usage: python alt_holdup.py <trial folder> <report file>

The residence distribution is built on a uniform passed-volume grid instead of from a transfer function:
the (n - 1) leading compartments are a gamma density (shape n - 1, scale Vc), the last compartment with its pockets
is the two-exponential response of the 2 x 2 exchange system (eigen-decomposition), the two are convolved, shifted by
the line plug and smeared over the dose box.  The salt term is the dosed mass times the cell factor; the inflow's own
conductivity enters as k0 + a F(u + w1) + b F(u - w2) with F the cumulative of the undosed distribution.
Parameters are fitted by least squares on unscaled volumes: first the early record without the late change (n = 1..4,
two starts), then the whole record for the two best n with the late change started at 40, 60 and 80 % of the record.  Moments in closed form:
    mean = Vp + n Vc + Vs,   variance = (n - 1) Vc^2 + (Vc + Vs)^2 + 2 Vs^2 / x.
"""
import csv
import json
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import fftconvolve
from scipy.stats import gamma


def read(folder):
    with open(f"{folder}/trial_card.json") as fh:
        card = json.load(fh)
    with open(f"{folder}/outlet_trace.csv", newline="") as fh:
        rows = np.array([[float(r["t_s"]), float(r["cond_uS_cm"]), float(r["inflow_L_s"])] for r in csv.DictReader(fh)])
    return card, rows[:, 0], rows[:, 1], rows[:, 2] / 1000.0


def last_cell(grid, vc, vs, x):
    a = np.array([[-(1.0 + x) / vc, x / vc], [x / vs, -x / vs]])
    lam, vec = np.linalg.eig(a)
    coef = np.linalg.solve(vec, np.array([1.0 / vc, 0.0]))
    out = (vec[0, :] * coef) @ np.exp(np.outer(lam, grid))
    return np.real(out)


class Trial:
    def __init__(self, folder):
        card, t, k, q = read(folder)
        vol = np.concatenate(([0.0], np.cumsum(np.diff(t) * 0.5 * (q[1:] + q[:-1]))))
        t0 = card["salt_dose"]["mid_time_s"]
        self.u = vol - np.interp(t0, t, vol)
        self.k = k
        self.salt = card["conductivity_cell"]["uS_cm_per_g_m3"] * 1000.0 * card["salt_dose"]["kg"]
        self.box = np.interp(t0, t, q) * card["salt_dose"]["duration_s"]
        self.step = 0.02
        self.span = 3.0 * (self.u[-1] + 1.0)
        self.grid = np.arange(0.0, self.span, self.step)

    def distribution(self, vp, vm, n, vs, x):
        vc = vm / n
        tail = last_cell(self.grid, vc, vs, x)
        if n > 1:
            lead = gamma.pdf(self.grid, n - 1, scale=vc)
            e = fftconvolve(lead, tail)[: self.grid.size] * self.step
            e -= 0.5 * self.step * (lead[0] * tail + lead * tail[0])      # trapezoid ends of the convolution
        else:
            e = tail
        e = np.interp(self.grid - vp, self.grid, e, left=0.0)
        width = max(int(round(self.box / self.step)), 1)
        smeared = np.convolve(e, np.ones(width) / width)[: self.grid.size]
        smeared = np.interp(self.grid + 0.5 * (width - 1) * self.step, self.grid, smeared)
        cum = np.concatenate(([0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * self.step)))
        return smeared, cum

    def misfit(self, p, n, mask=None):
        vp, vm, vs, x, w1 = p[:5]
        e, cum = self.distribution(vp, vm, n, vs, x)
        u = self.u if mask is None else self.u[mask]
        k = self.k if mask is None else self.k[mask]
        rest = k - self.salt * np.interp(u, self.grid, e, left=0.0, right=0.0)
        cols = [np.ones_like(u), np.interp(u + w1, self.grid, cum, left=0.0, right=cum[-1])]
        if len(p) > 5:
            cols.append(np.interp(u - p[5], self.grid, cum, left=0.0, right=cum[-1]))
        cols = np.column_stack(cols)
        coef = np.linalg.lstsq(cols, rest, rcond=None)[0]
        return cols @ coef - rest


def evaluate(folder):
    rec = Trial(folder)
    end = rec.u[-1]
    excess = rec.k - np.median(rec.k[rec.u < 0])
    size = rec.salt / float(np.max(excess)) * 0.6
    lo = [0.0, 0.05 * size, 0.01 * size, 0.01, 0.0, 0.0]
    hi = [0.5 * size, 2.5 * size, 2.5 * size, 5.0, 2.5 * size, end]
    early = rec.u < 1.6 * size
    first = []
    for n in (1, 2, 3, 4):
        for vs0, x0 in ((0.4, 0.2), (0.9, 0.1)):
            p0 = [0.05 * size, 0.6 * size, vs0 * size, x0, 0.5 * size]
            fit = least_squares(rec.misfit, p0, args=(n, early), bounds=(lo[:5], hi[:5]), x_scale="jac", max_nfev=150)
            first.append((fit.cost, n, fit.x))
    first.sort(key=lambda r: r[0])
    best = None
    for _, n, p in first[:2]:
        for share in (0.4, 0.6, 0.8):
            p0 = np.clip(np.r_[p, share * end], np.array(lo) + 1e-9, np.array(hi) - 1e-9)
            fit = least_squares(rec.misfit, p0, args=(n,), bounds=(lo, hi), x_scale="jac", max_nfev=300)
            if best is None or fit.cost < best[0].cost:
                best = (fit, n)
    fit, n = best
    vp, vm, vs, x = fit.x[:4]
    vc = vm / n
    mean = vp + n * vc + vs
    var = (n - 1) * vc ** 2 + (vc + vs) ** 2 + 2 * vs ** 2 / x
    return float(mean), float(var / mean ** 2)


def main(folder, target):
    hold, spread = evaluate(folder)
    with open(target, "w") as fh:
        fh.write(json.dumps({"holdup_volume_m3": hold, "spread_ratio": spread}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
