"""Second, separately written evaluation of a brine-pulse record (no flow model at all).

usage: python second_sweep.py <vessel folder> <result file>

The residence distribution is taken point by point from the readings: E(u) = (k(u) - b(u)) / (c M), with c M the
probe factor times the salt mass from the sheet, u the passed feed volume since the dose, and the feed term
b(u) = b0 + dk * F(u + w), F the running integral of E itself (closed by an exponential tail fitted from
the 5 % point to the end of the record). E and F are iterated to a fixed point for each trial (b0, dk, w); those three are fitted to
the pre-dose readings, to the last 5 % of the record and to unit area of E. Moments of E, read up to its 1 % point and
closed by the fitted tail beyond it, give the answer.
"""
import csv
import json
import sys

import numpy as np
from scipy.optimize import least_squares


def read(folder):
    with open(f"{folder}/dosing_sheet.json") as fh:
        sheet = json.load(fh)
    with open(f"{folder}/outlet_log.csv", newline="") as fh:
        rows = [(float(r["elapsed_s"]), float(r["outlet_mS_cm"]), float(r["feed_m3_h"])) for r in csv.DictReader(fh)]
    return sheet, np.array(rows).T


def exp_tail(u, e):
    top = int(np.argmax(e))
    below = np.nonzero(e[top:] < 0.05 * e[top])[0]
    if below.size == 0:
        return 0.0, 1.0
    first = top + int(below[0])
    x, y = u[first:], e[first:]
    if x.size < 12:
        return 0.0, 1.0
    lam = 1.0 / max(x[-1] - x[0], 1e-6)
    for _ in range(40):
        g = np.exp(-lam * (x - x[0]))
        amp = float(g @ y / (g @ g))
        resid = y - amp * g
        d = -amp * (x - x[0]) * g
        step = float(d @ resid / (d @ d))
        lam = float(np.clip(lam + step, 1e-6, 50.0 / max(x[-1] - x[0], 1e-6)))
    return amp * np.exp(lam * x[0]), lam


def extended(u, e, cut=False):
    amp, lam = exp_tail(u, e)
    if cut:
        top = int(np.argmax(e))
        low = np.nonzero(e[top:] < 0.01 * e[top])[0]
        if low.size:
            u, e = u[:top + int(low[0])], e[:top + int(low[0])]
    extra = u[-1] + np.linspace(0.0, 30.0 / lam, 400)[1:]
    return np.r_[u, extra], np.r_[e, amp * np.exp(-lam * extra)], (amp, lam)


def running(x, y):
    return np.r_[0.0, np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))]


class Record:
    def __init__(self, folder):
        sheet, (t, k, q) = read(folder)
        vol = np.r_[0.0, np.cumsum(np.diff(t) * 0.5 * (q[1:] + q[:-1]) / 3600.0)]
        self.u = vol - np.interp(sheet["dose"]["dosed_at_s"], t, vol)
        self.k = k
        self.amp = sheet["dose"]["mass_kg"] * sheet["probe"]["factor_mS_cm_per_kg_m3"]
        self.pre = self.u < 0
        self.post = ~self.pre
        self.late = np.arange(self.u.size) >= int(0.95 * self.u.size)

    def distribution(self, b0, dk, w):
        up, kp = self.u[self.post], self.k[self.post]
        e = (kp - b0) / self.amp
        for _ in range(8):
            xu, xe, tail = extended(up, np.clip(e, 0.0, None))
            f = running(xu, xe)
            e = (kp - b0 - dk * np.interp(up + w, xu, f, right=f[-1])) / self.amp
        xu, xe, tail = extended(up, e)
        return xu, xe, running(xu, xe), tail

    def residuals(self, z):
        b0, dk, w = z
        xu, xe, f, (amp, lam) = self.distribution(b0, dk, w)
        feed = lambda uu: b0 + dk * np.interp(uu + w, xu, f, left=0.0, right=f[-1])
        r_pre = (self.k[self.pre] - feed(self.u[self.pre])) / 0.004
        ul = self.u[self.late]
        r_end = (self.k[self.late] - feed(ul) - self.amp * amp * np.exp(-lam * ul)) / 0.004
        r_area = np.array([(f[-1] - 1.0) * 400.0])
        return np.r_[r_pre, r_end, r_area]


def evaluate(folder):
    rec = Record(folder)
    b_start = float(np.median(rec.k[rec.pre]))
    span = float(rec.u[-1]) / 6.0
    best = None
    for w0 in (0.2, 0.5, 1.0):
        for dk0 in (-0.05, 0.05):
            fit = least_squares(rec.residuals, [b_start, dk0, w0 * span], bounds=([0.0, -1.0, 0.0], [5.0, 1.0, 3 * span]),
                                diff_step=1e-3)
            if best is None or fit.cost < best.cost:
                best = fit
    xu, xe, _, _ = rec.distribution(*best.x)
    xu, xe, _ = extended(xu[:np.count_nonzero(rec.post)], xe[:np.count_nonzero(rec.post)], cut=True)
    m0, m1, m2 = (np.trapezoid(xe * xu ** p, xu) for p in range(3))
    mean = m1 / m0
    return float(mean), float(mean ** 2 / (m2 / m0 - mean ** 2))


def main(folder, target):
    swept, stages = evaluate(folder)
    with open(target, "w") as fh:
        fh.write(json.dumps({"swept_volume_m3": swept, "mixing_stages": stages}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
