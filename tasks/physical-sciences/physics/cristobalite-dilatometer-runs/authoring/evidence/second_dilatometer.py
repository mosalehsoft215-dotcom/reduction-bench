"""Second, separately written reduction of a dilatometer run (profile fit of the lag over the second cycle).

The specimen temperature is the thermocouple run through the lag dT/dt = (T_tc - T) / tau with tau = tau_300 (300 / T)^3
by fourth-order Runge-Kutta steps on the linearly interpolated thermocouple (each reading interval cut into sub-steps
of at most half the current tau).  The second cycle starts 1.6 min before the record leaves 300 K again after the
first cooling (found in the record, whatever the pause between the cycles); its ramps are the readings changing by more than 0.6 of the programme
rate.  For a trial tau_300 both branches of the second cycle below 600 K are fitted by a
quadratic background, a constant offset of the cooling branch (the creep of the body near the top) and
S Phi((T - Ti) / w) against that temperature (and, in a second profile scanned with every point polished, S Phi((T - Tc) / w) on the cooling branch,
the inversion running back at a temperature of its own, kept when its best residual is under half the first's and it runs back at least 2 K lower),
with or without a moisture contraction -A Phi((T - Tm) / wm) on the
heating branch (kept when it lowers the sum of squares by a fifth), (Ti, w) by Nelder-Mead from the best point of a coarse grid, the rest by
linear least squares; tau_300 is the value with the least residual of that fit (golden-section search on 0-8 min
after a coarse scan).  The mean expansion 320-440 K comes from the fitted curve.  The first cycle is used only for
the thermocouple's drift with use: on each heating the inversion's place (in the cells-corrected thermocouple) is a
step S Phi((T - place) / width) on a straight background within 20 K, (place, width) by Nelder-Mead from the steepest
rise of the strain on a 0.5 K grid, the rest by linear least squares; the places on the first and the second heating
give the drift: their difference over (place - 300 K) is the drift
D on the second heating, and D (T - 300 K) is taken off the thermocouple along the whole record in proportion to its
time spent hot (trapezoid sums of exp(-E/R (1/T - 1/640 K)), E = 500 kJ/mol, the weighting evaluated at the corrected
temperature by fixed-point iteration).
Before that, the day's checks, found differently from the reference: a fixed-point cell's plateau is the most populated
0.1 K bin of its readings (the cell lingers there while it melts), its value the mean of the readings within 0.15 K of
the bin's centre; the thermocouple's error is the line through the two cells' errors.  The reference bar's temperature is
the corrected thermocouple run through the same Runge-Kutta lag with a time constant of its own; every reading of the
reference run enters a least-squares fit of displacement - L_ref e_ref(T_bar) by u0 + a (T - 300) + b (T - 300)^2
(T the corrected thermocouple, which the push rod follows), the bar's tau_300 chosen by golden-section search on
0-10 min after a coarse scan to give the least residual, and a (T - 300) + b (T - 300)^2 is taken off the specimen's
displacement.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import ndtr


def specimen(t, T, tau300):
    if tau300 <= 1e-6:
        return T.copy()
    f = lambda temp, drive: (drive - temp) / (tau300 * (300.0 / temp) ** 3)
    Ts = np.empty_like(T)
    Ts[0] = y = T[0]
    for i in range(1, T.size):
        n = max(1, math.ceil((t[i] - t[i - 1]) / (0.5 * tau300 * (300.0 / y) ** 3)))
        h = (t[i] - t[i - 1]) / n
        for j in range(n):
            u0 = T[i - 1] + (T[i] - T[i - 1]) * j / n
            u1 = T[i - 1] + (T[i] - T[i - 1]) * (j + 1) / n
            um = 0.5 * (u0 + u1)
            k1 = f(y, u0)
            k2 = f(y + 0.5 * h * k1, um)
            k3 = f(y + 0.5 * h * k2, um)
            k4 = f(y + h * k3, u1)
            y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        Ts[i] = y
    return Ts


def checks(folder):
    fp = np.genfromtxt(folder / "fixed_points.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    which = np.array([str(c) for c in fp["cell"]])
    errs = []
    for cell, tm in (("indium", 429.75), ("tin", 505.08)):
        y = fp["thermocouple_k"][which == cell]
        edges = np.arange(y.min(), y.max() + 0.1, 0.1)
        hist, edges = np.histogram(y, bins=edges)
        centre = edges[np.argmax(hist)] + 0.05
        errs.append(float(np.mean(y[np.abs(y - centre) <= 0.15])) - tm)
    slope = (errs[1] - errs[0]) / (505.08 - 429.75)
    fix = lambda T: T - (errs[0] + slope * (T - 429.75))  # noqa: E731
    bar = json.loads((folder / "reference_bar.json").read_text())
    ref = np.genfromtxt(folder / "reference.csv", delimiter=",", names=True)
    tt, TT, dd = ref["minute"], fix(ref["thermocouple_k"]), ref["displacement_um"]
    xf = TT - 300.0
    design = np.column_stack([np.ones_like(xf), xf, xf ** 2])
    e1, e2 = bar["expansion_ppm"]["c1_per_k"], bar["expansion_ppm"]["c2_per_k2"]

    def fit_bar(tau):
        x = specimen(tt, TT, tau) - 300.0
        left = dd - bar["length_mm"] * 1e-3 * (e1 * x + e2 * x ** 2)
        coef = np.linalg.lstsq(design, left, rcond=None)[0]
        return float(np.sum((design @ coef - left) ** 2)), coef

    scan = np.linspace(0.0, 10.0, 21)
    best = scan[int(np.argmin([fit_bar(g)[0] for g in scan]))]
    lo, hi = max(best - 0.5, 0.0), best + 0.5
    gr = (math.sqrt(5) - 1) / 2
    for _ in range(40):
        m1, m2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if fit_bar(m1)[0] < fit_bar(m2)[0]:
            hi = m2
        else:
            lo = m1
    coef = fit_bar(0.5 * (lo + hi))[1]
    return fix, lambda T: coef[1] * (T - 300.0) + coef[2] * (T - 300.0) ** 2


def peak_place(T, eps, sel, near=None):
    """Place of the inversion on one heating: the steepest rise of strain against temperature (0.5 K grid) as a start,
    then Nelder-Mead on (place, width) of S Phi((T - place) / width) over a straight background within 20 K of it
    (S and the background by linear least squares)."""
    if near is None:
        order = np.argsort(T[sel])
        grid = np.arange(430.0, 600.0, 0.5)
        near = float(grid[int(np.argmax(np.gradient(np.interp(grid, T[sel][order], eps[sel][order]), grid)))])
    use = sel & (np.abs(T - near) < 20.0)
    x, y = T[use], eps[use]

    def cost(q):
        B = np.vstack([np.ones_like(x), x - near, ndtr((x - q[0]) / max(abs(q[1]), 0.5))]).T
        c, *_ = np.linalg.lstsq(B, y, rcond=None)
        return float(np.sum((B @ c - y) ** 2))
    q = min(((near + d, 5.0) for d in (-6.0, -3.0, 0.0, 3.0, 6.0)), key=cost)
    return float(minimize(cost, q, method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-12, "maxiter": 2000}).x[0])


def undrift(t, T, eps, first_up, second_up):
    x1 = peak_place(T, eps, first_up)
    x2 = peak_place(T, eps, second_up, x1)
    d2 = (x2 - x1) / (x1 - 300.0)
    if not np.isfinite(d2) or d2 <= 0.0:
        return T
    dt = np.diff(t)
    corrected = T.copy()
    for _ in range(5):
        w = np.exp(-500e3 / 8.314462618 * (1.0 / corrected - 1.0 / 640.0))
        hot = np.r_[0.0, np.cumsum(0.5 * (w[1:] + w[:-1]) * dt)]
        k = d2 / float(np.mean(hot[second_up]))
        corrected = 300.0 + (T - 300.0) / (1.0 + k * hot)
    return corrected


def reduce_run(folder):
    folder = Path(folder)
    meta = json.loads((folder / "run.json").read_text())
    rec = np.genfromtxt(folder / "record.csv", delimiter=",", names=True)
    fix, residual = checks(folder)
    t, T = rec["minute"], fix(rec["thermocouple_k"])
    d = rec["displacement_um"] - residual(T)
    eps = d / (meta["specimen_length_mm"] * 1e3) * 1e6
    rate = meta["programme_rate_k_per_min"]
    raw = rec["thermocouple_k"]
    top1 = raw[: raw.size // 2].max()                           # the first passage over the top
    hot = np.flatnonzero(raw > top1 - 3.0)
    cut = hot[np.argmax(np.diff(np.r_[hot, hot[-1]]) > 50)] if np.any(np.diff(hot) > 50) else hot[0]
    cool_end = cut + int(np.argmax(raw[cut:] < 302.0))          # end of the first cooling
    restart = cool_end + int(np.argmax(raw[cool_end:] > 305.0))  # the record leaves 300 K again
    second = np.arange(t.size) >= restart - 8
    hot2 = np.flatnonzero(second & (raw > top1 - 3.0))
    speed = np.gradient(raw, t)
    speed = np.convolve(speed, np.ones(11) / 11, mode="same")
    branch_all = (second & (np.arange(t.size) > hot2[-1])).astype(float)
    on_ramp = second & (np.abs(speed) > 0.6 * rate)
    rising = on_ramp & (branch_all == 0)
    T = undrift(t, T, eps, ~second & (speed > 0.6 * rate) & (np.arange(t.size) < cut), rising)

    def fit(Ts, polish=True, free=False):
        use = on_ramp & (Ts < 600.0)
        gT, curve, br, up = Ts[use], eps[use], branch_all[use], rising[use].astype(float)
        x = gT - 300.0

        def solve(p):
            wd = max(abs(p[1]), 0.3)
            tc = p[2] if free else p[0]                      # free: the inversion runs back at Tc on cooling
            step = (1 - br) * ndtr((gT - p[0]) / wd) + br * ndtr((gT - tc) / wd)
            cols = [np.ones_like(x), br, x, x ** 2, step]
            if len(p) == 5:                                  # moisture given off on this heating
                cols.append(-up * ndtr((gT - p[3]) / max(abs(p[4]), 2.0)))
            B = np.vstack(cols).T
            c, *_ = np.linalg.lstsq(B, curve, rcond=None)
            return c, float(np.sum((B @ c - curve) ** 2))

        steps = (0.0, 5.0, 10.0, 15.0) if free else (0.0,)
        starts = [(ti, w, ti - dh) for ti in np.arange(440.0, 596.0, 4.0) for w in (3.0, 6.0) for dh in steps]
        p = min(starts, key=lambda q: solve(q)[1])
        if polish:
            p = minimize(lambda q: solve(q)[1], p, method="Nelder-Mead",
                         options={"xatol": 1e-4, "fatol": 1e-10, "maxiter": 1500}).x
        dry = solve(p)
        wet_starts = [(p[0], p[1], p[2], tm, 15.0) for tm in (370.0, 395.0, 420.0, 445.0)]
        q = min(wet_starts, key=lambda r: solve(r)[1])
        if polish:
            q = minimize(lambda r: solve(r)[1], q, method="Nelder-Mead",
                         options={"xatol": 1e-4, "fatol": 1e-10, "maxiter": 2000}).x
        wet = solve(q)
        if wet[1] < 0.8 * dry[1] and wet[0][-1] > 5.0:
            return wet[1], q, wet[0]
        return dry[1], p, dry[0]

    def profile(free):
        scan = np.linspace(0.0, 8.0, 9 if free else 17)
        costs = [fit(specimen(t, T, s), polish=free, free=free)[0] for s in scan]
        k = int(np.argmin(costs))
        lo, hi = scan[max(k - 1, 0)], scan[min(k + 1, scan.size - 1)]
        tau = minimize_scalar(lambda s: fit(specimen(t, T, s), free=free)[0], bounds=(lo, hi), method="bounded",
                              options={"xatol": 2e-3}).x
        return fit(specimen(t, T, tau), free=free) + (tau,)

    tied = profile(False)
    loose = profile(True)                  # a cooling inversion of its own, kept when it lowers the residual by half
    _, p, c = loose[:3] if (loose[0] < 0.5 * tied[0] and loose[1][0] - loose[1][2] >= 2.0) else tied[:3]
    ti, w = p[0], max(abs(p[1]), 0.3)
    strain = lambda temp: c[2] * (temp - 300) + c[3] * (temp - 300) ** 2 + c[4] * ndtr((temp - ti) / w)
    return {"inversion_temperature_k": float(ti),
            "mean_expansion_320_440_ppm_per_k": float((strain(440.0) - strain(320.0)) / 120.0)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_run(sys.argv[1])) + "\n")
