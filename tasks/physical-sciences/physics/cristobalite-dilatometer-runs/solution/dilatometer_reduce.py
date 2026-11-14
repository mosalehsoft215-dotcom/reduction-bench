"""Push-rod dilatometer run -> cristobalite inversion temperature and mean expansion coefficient 320-440 K.

    python3 dilatometer_reduce.py <run folder> <out.json>

0. The day's checks, before anything else.  Thermocouple: each fixed-point cell (fixed_points.csv) holds at its
   melting point while it melts, so the longest stretch of its readings whose slope over 1 min stays below
   0.08 K/min is the plateau; the median of the plateau less 429.75 K (indium) or 505.08 K (tin) is the
   thermocouple's error there, linear in temperature through the two points, and every thermocouple reading is
   corrected by it.  Push rod: the alumina reference bar (reference.csv, certificate in reference_bar.json) follows
   the furnace with the same radiative lag as any bar, with a time constant of its own, and the push rod and holder
   follow the furnace temperature T (the corrected thermocouple).  Over the whole reference run its displacement is
   L_ref (c1 x + c2 x^2) 1e-6 + a (T - 300) + b (T - 300)^2 um, x = T_bar - 300 with T_bar the corrected
   thermocouple passed through the lag of step 1; for a trial tau_300 of the bar a and b are linear, and tau_300 is
   the value with the least sum of squares (bounded 1-D search around the best of a 0-10 min scan).
   r(T) = a (T - 300) + b (T - 300)^2 at the corrected thermocouple temperature is taken off the specimen's
   displacement.  (rod="ends" instead takes the bar at the furnace temperature over the last minute of every hold,
   and leaves the displacement as logged when the run has fewer than two holds above 320 K.)
1. Specimen temperature: the corrected thermocouple passed through the notes' lag, dT/dt = (T_tc - T) / tau with
   tau = tau_300 (300 K / T)^3, stepped exactly between readings for a thermocouple changing linearly over the step.
2. Only the second cycle is used, and only below 600 K: the first heating of a specimen need not be repeatable, and
   near the top the body may creep, an irreversible contraction that leaves every later branch a constant amount
   lower (it has stopped again by 600 K).  So the cooling branch carries an offset of its own.
3. tau_300 is the value that makes the heating and cooling branches of the second cycle fall on one strain curve up
   to that constant offset (the variance of their difference on a common temperature grid, bounded 1-D search).
4. Strain (ppm of the length) against specimen temperature is fitted by e0 + o [cooling] + c1 (T - 300) +
   c2 (T - 300)^2 + S Phi((T - Ti) / w); for trial (Ti, w) the rest is linear, (Ti, w) by least squares.
5. Mean expansion 320-440 K from the fitted reversible curve.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import lfilter
from scipy.special import ndtr


def plateau(t, y):
    """Median of the longest stretch whose slope over 1 min stays below 0.08 K/min (a cell holding while it melts)."""
    k = max(1, int(round(1.0 / np.median(np.diff(t)))))
    slope = np.abs(y[k:] - y[:-k]) / (t[k:] - t[:-k])
    flat = np.r_[slope < 0.08, np.zeros(k, bool)]
    best, start = (0, 0), None
    for i, f in enumerate(np.r_[flat, False]):
        if f and start is None:
            start = i
        elif not f and start is not None:
            best = max(best, (i - start, start))
            start = None
    n, s0 = best
    return float(np.median(y[s0:s0 + n + k]))


def day_checks(folder, thermocouple=True, rod=True):
    """(thermocouple error as a function of the reading, push-rod residual as a function of temperature)."""
    folder = Path(folder)
    err = lambda T: np.zeros_like(np.asarray(T, float))  # noqa: E731
    if thermocouple:
        fp = np.genfromtxt(folder / "fixed_points.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
        cell = np.array([str(c) for c in fp["cell"]])
        o_in = plateau(fp["minute"][cell == "indium"], fp["thermocouple_k"][cell == "indium"]) - 429.75
        o_sn = plateau(fp["minute"][cell == "tin"], fp["thermocouple_k"][cell == "tin"]) - 505.08
        err = lambda T: o_in + (o_sn - o_in) * (np.asarray(T, float) - 429.75) / (505.08 - 429.75)  # noqa: E731
    resid = lambda T: np.zeros_like(np.asarray(T, float))  # noqa: E731
    if rod and rod != "ends":
        bar = json.loads((folder / "reference_bar.json").read_text())
        ref = np.loadtxt(folder / "reference.csv", delimiter=",", skiprows=1)
        t, T = ref[:, 0], ref[:, 1] - err(ref[:, 1])
        c = bar["expansion_ppm"]
        X = T - 300.0
        A = np.column_stack([X, X ** 2])
        scale = bar["length_mm"] * 1e3 * 1e-6

        def solve(tau):
            x = lagged(t, T, tau) - 300.0
            left = ref[:, 2] - scale * (c["c1_per_k"] * x + c["c2_per_k2"] * x ** 2)
            coef = np.linalg.lstsq(A, left, rcond=None)[0]
            return float(np.sum((A @ coef - left) ** 2)), coef

        grid = np.arange(0.0, 10.01, 0.25)
        g0 = grid[int(np.argmin([solve(g)[0] for g in grid]))]
        tau = minimize_scalar(lambda g: solve(g)[0], bounds=(max(g0 - 0.25, 0.0), g0 + 0.25), method="bounded",
                              options={"xatol": 1e-4}).x
        a, b = solve(tau)[1]
        resid = lambda T: a * (np.asarray(T, float) - 300.0) + b * (np.asarray(T, float) - 300.0) ** 2  # noqa: E731
    elif rod == "ends":
        bar = json.loads((folder / "reference_bar.json").read_text())
        ref = np.loadtxt(folder / "reference.csv", delimiter=",", skiprows=1)
        t, T = ref[:, 0], ref[:, 1] - err(ref[:, 1])
        k = max(1, int(round(1.0 / np.median(np.diff(t)))))
        flat = np.r_[np.abs(T[k:] - T[:-k]) / (t[k:] - t[:-k]) < 0.3, np.zeros(k, bool)]
        pts, start = [], None
        for i, f in enumerate(np.r_[flat, False]):
            if f and start is None:
                start = i
            elif not f and start is not None:
                if t[i - 1] - t[start] >= 2.0:
                    end = t[i - 1]
                    sel = (t >= end - 1.0) & (t <= end)
                    pts.append((T[sel].mean(), ref[sel, 2].mean()))
                start = None
        if sum(1 for T_, _ in pts if T_ > 320.0) < 2:   # no holds above 300 K to read: displacement as logged
            return err, resid
        Tk, dk = np.array(pts).T
        c = bar["expansion_ppm"]
        x = Tk - 300.0
        left = dk - bar["length_mm"] * 1e3 * 1e-6 * (c["c1_per_k"] * x + c["c2_per_k2"] * x ** 2)
        a, b = np.linalg.lstsq(np.column_stack([x, x ** 2]), left, rcond=None)[0]
        resid = lambda T: a * (np.asarray(T, float) - 300.0) + b * (np.asarray(T, float) - 300.0) ** 2  # noqa: E731
    return err, resid


def load(folder, thermocouple=True, rod=True):
    folder = Path(folder)
    sheet = json.loads((folder / "run.json").read_text())
    d = np.loadtxt(folder / "record.csv", delimiter=",", skiprows=1)
    err, resid = day_checks(folder, thermocouple, rod)
    tc = d[:, 1] - err(d[:, 1])
    disp = d[:, 2] - resid(tc)
    return sheet, d[:, 0], tc, disp / (sheet["specimen_length_mm"] * 1e3) * 1e6


def lagged(t, tc, tau, radiative=True):
    """Specimen temperature: exact exponential step towards each thermocouple reading with the time constant at the
    current specimen temperature (tau_300 (300 / T)^3), or a constant tau for radiative=False."""
    if tau <= 0:
        return tc.copy()
    if not radiative:
        a = math.exp(-np.median(np.diff(t)) / tau)
        return lfilter([1 - a], [1, -a], tc, zi=[a * tc[0]])[0]
    out = np.empty_like(tc)
    temp = tc[0]
    out[0] = temp
    dt = np.diff(t)
    for i in range(1, tc.size):            # exact for a thermocouple changing linearly between readings
        tk = tau * (300.0 / temp) ** 3
        m = (tc[i] - tc[i - 1]) / dt[i - 1]
        temp = tc[i] - m * tk + (temp - tc[i - 1] + m * tk) * math.exp(-dt[i - 1] / tk)
        out[i] = temp
    return out


def segments(t, tc, rate):
    """Cycle number from the record (the second cycle from 2 min before the record leaves 300 K again after the first
    cooling and the pause) and the heating and cooling ramps."""
    idx = np.arange(tc.size)
    i_top = int(np.argmax(tc > tc[: tc.size // 2].max() - 5.0))      # the first passage over the top
    low = idx[(idx > i_top) & (tc < 303.0)]
    rise = idx[(idx > (low[0] if low.size else tc.size)) & (tc > 304.0)]
    cyc = (idx >= (rise[0] - 10 if rise.size else tc.size)).astype(int)
    slope = np.gradient(tc, t)
    smooth = np.convolve(slope, np.ones(9) / 9, mode="same")
    heating = smooth > 0.5 * rate
    cooling = smooth < -0.5 * rate
    return cyc, heating, cooling


def fit_curve(T, e, cool=None, heat=None, back=None):
    """Background, inversion step and (cool) a cooling-branch offset; with back (a cooling-branch mask) the inversion
    runs back on cooling at a temperature Tc of its own (thermal hysteresis), the same size and width; with heat (a
    heating-branch mask) also a moisture contraction -A Phi((T - Tm) / wm) on the heating branch, kept only when it
    lowers the sum of squares by a fifth.  Returns Ti (heating), w, c1, c2, S, whether moisture was kept, Tc."""
    extra = [] if cool is None else [cool.astype(float)]

    hyst = back is not None

    def basis(p):
        step = ndtr((T - p[0]) / abs(p[1]))
        if hyst:
            step = np.where(back > 0, ndtr((T - p[2]) / abs(p[1])), step)
        cols = [np.ones_like(T), *extra, T - 300.0, (T - 300.0) ** 2, step]
        if len(p) == 4 + hyst:
            cols.append(-heat * ndtr((T - p[2 + hyst]) / abs(p[3 + hyst])))
        return np.vstack(cols).T

    def resid(p):
        B = basis(p)
        c = np.linalg.lstsq(B, e, rcond=None)[0]
        return B @ c - e

    best = None
    for ti0 in np.arange(440.0, 590.0, 10.0):
        if hyst:
            sol = least_squares(resid, [ti0, 5.0, ti0 - 4.0], bounds=([420.0, 1.0, 395.0], [600.0, 20.0, 600.0]))
        else:
            sol = least_squares(resid, [ti0, 5.0], bounds=([420.0, 1.0], [600.0, 20.0]))
        if best is None or sol.cost < best.cost:
            best = sol
    if heat is not None:
        wet = None
        for tm0 in (380.0, 400.0, 420.0, 440.0):
            head = list(best.x)
            lo = [420.0, 1.0] + ([395.0] if hyst else []) + [340.0, 4.0]
            hi = [600.0, 20.0] + ([600.0] if hyst else []) + [470.0, 40.0]
            sol = least_squares(resid, head + [tm0, 15.0], bounds=(lo, hi))
            if wet is None or sol.cost < wet.cost:
                wet = sol
        amp = np.linalg.lstsq(basis(wet.x), e, rcond=None)[0][-1]
        if wet.cost < 0.8 * best.cost and amp > 5.0:
            best = wet
    c = np.linalg.lstsq(basis(best.x), e, rcond=None)[0]
    k = len(extra)
    return (best.x[0], abs(best.x[1]), c[k + 1], c[k + 2], c[k + 3], len(best.x) == 4 + hyst,
            best.x[2] if hyst else best.x[0])


E_DRIFT = 500e3                         # J/mol, the notes' weighting of the thermocouple's drift with use


def heating_step(tc, e, sel, guess):
    """Inversion position (in thermocouple temperature) of one heating: a step S Phi((T - x0) / w) on a quadratic,
    fitted to the heating readings within 30 K of the guess."""
    win = sel & (np.abs(tc - guess) < 30.0)
    x, y = tc[win], e[win]

    def resid(p):
        B = np.vstack([np.ones_like(x), x - guess, (x - guess) ** 2, ndtr((x - p[0]) / abs(p[1]))]).T
        return B @ np.linalg.lstsq(B, y, rcond=None)[0] - y
    best = min((least_squares(resid, [guess + d, 5.0], bounds=([guess - 25.0, 1.0], [guess + 25.0, 20.0]))
                for d in (-8.0, 0.0, 8.0)), key=lambda r: r.cost)
    return float(best.x[0])


def undrift(t, tc, e, cyc, heat):
    """The thermocouple's drift with use: its reading gains D (T - 300 K) with D = k H(t), H the time spent hot
    weighted by exp(-E/R (1/T - 1/640 K)).  The first heating is read with the cells' calibration as it stands; the
    second heating's inversion, the same reversible step at the same lag, then gives D there, hence k."""
    def hot_time(temp):
        return np.cumsum(np.exp(-E_DRIFT / 8.314462618 * (1.0 / temp - 1.0 / 640.0)) * np.gradient(t))

    h1, h2 = (cyc == 0) & heat & (tc < 600), (cyc == 1) & heat & (tc < 600)
    order = np.argsort(tc[h1])
    grid = np.arange(420.0, 600.0, 1.0)
    rise = np.diff(np.interp(grid, tc[h1][order], e[h1][order]))
    guess = float(grid[int(np.argmax(rise))]) + 0.5
    x1 = heating_step(tc, e, h1, guess)
    x2 = heating_step(tc, e, h2, x1)
    d2 = max((x2 - x1) / (x1 - 300.0), 0.0)
    if not np.isfinite(d2) or d2 == 0.0:
        return tc
    true = tc
    for _ in range(4):     # the weighting follows the true temperature: iterate the correction
        hot = hot_time(true)
        k = d2 / float(np.median(hot[h2 & (np.abs(tc - x2) < 30.0)]))
        true = 300.0 + (tc - 300.0) / (1.0 + k * hot)
    return true


def reduce_run(folder, lag="radiative", cycles="second", offset=True, thermocouple=True, rod=True, moisture=True,
               hysteresis=True, drift=True):
    """lag: 'radiative', 'constant' or None; cycles: 'second', 'first' or 'second-heating'; offset=False fits one curve
    through both branches; thermocouple=False / rod=False leave the thermocouple or the displacement as logged;
    moisture=False takes the second heating as free of moisture; drift=False keeps the cells' calibration through
    the whole run (no drift with use); hysteresis=False fits one inversion temperature to
    heating and cooling alike and takes the lag from the whole range.  The
    alternatives are the shortcuts tabulated by authoring/evidence/bench.py (rod="ends": the reference bar taken at
    the furnace temperature at the end of every hold)."""
    sheet, t, tc, e = load(folder, thermocouple, rod)
    rate = sheet["programme_rate_k_per_min"]
    cyc, heat, cool = segments(t, tc, rate)
    if drift:
        tc = undrift(t, tc, e, cyc, heat)
    radiative = lag == "radiative"
    if lag:
        def mismatch(tau):
            Ts = lagged(t, tc, tau, radiative)
            h, c = (cyc == 1) & heat & (Ts > 320) & (Ts < 595), (cyc == 1) & cool & (Ts > 320) & (Ts < 595)
            grid = np.linspace(low_end, 585.0, 160)
            grid = grid[(grid < skip[0]) | (grid > skip[1])]
            oh, oc = np.argsort(Ts[h]), np.argsort(Ts[c])
            diff = np.interp(grid, Ts[h][oh], e[h][oh]) - np.interp(grid, Ts[c][oc], e[c][oc])
            return float(np.var(diff) if offset else np.mean(diff ** 2))
        low_end, skip = 330.0, (0.0, 0.0)
        tau = minimize_scalar(mismatch, bounds=(0.0, 8.0), method="bounded", options={"xatol": 1e-3}).x
        Ts = lagged(t, tc, tau, radiative)
    else:
        Ts = tc
    if cycles == "first":
        sel = (cyc == 0) & heat & (Ts < 600)
    elif cycles == "second-heating":
        sel = (cyc == 1) & heat & (Ts < 600)
    else:
        sel = (cyc == 1) & (heat | cool) & (Ts < 600)
    second = offset and cycles == "second"
    hyst = second and hysteresis

    def curve():
        return fit_curve(Ts[sel], e[sel], cool[sel] if second else None,
                         heat[sel].astype(float) if (second and moisture) else None,
                         cool[sel].astype(float) if hyst else None)

    ti, w, c1, c2, s, wet, tcool = curve()
    if hyst and lag:         # a hysteresis is kept only when the cooling inversion, the lag found away from the
        low_end = 455.0 if wet else 330.0          # steps, still runs back at least 2 K lower
        skip = (min(ti, tcool) - 5.0 * w, max(ti, tcool) + 5.0 * w)
        tau_h = minimize_scalar(mismatch, bounds=(0.0, 8.0), method="bounded", options={"xatol": 1e-3}).x
        Ts_h = lagged(t, tc, tau_h, radiative)
        sel_h = (cyc == 1) & (heat | cool) & (Ts_h < 600)
        Ts, sel = Ts_h, sel_h
        ti, w, c1, c2, s, wet, tcool = curve()
        if ti - tcool < 2.0:
            hyst = False
            low_end, skip = (455.0 if wet else 330.0), (0.0, 0.0)
            tau = minimize_scalar(mismatch, bounds=(0.0, 8.0), method="bounded", options={"xatol": 1e-3}).x
            Ts = lagged(t, tc, tau, radiative)
            sel = (cyc == 1) & (heat | cool) & (Ts < 600)
            ti, w, c1, c2, s, wet, tcool = curve()
    if lag and (wet or hyst):
        for _ in range(2):   # the lag again from the background: above 455 K under moisture, away from the steps
            low_end = 455.0 if wet else 330.0
            skip = (min(ti, tcool) - 5.0 * w, max(ti, tcool) + 5.0 * w) if hyst else (0.0, 0.0)
            tau = minimize_scalar(mismatch, bounds=(0.0, 8.0), method="bounded", options={"xatol": 1e-3}).x
            Ts = lagged(t, tc, tau, radiative)
            sel = (cyc == 1) & (heat | cool) & (Ts < 600)
            ti, w, c1, c2, s, wet, tcool = curve()
    lo, hi = 320.0, 440.0
    cte = c1 + c2 * ((hi - 300) ** 2 - (lo - 300) ** 2) / (hi - lo) + s * (ndtr((hi - ti) / w) - ndtr((lo - ti) / w)) / (hi - lo)
    return {"inversion_temperature_k": float(ti), "mean_expansion_320_440_ppm_per_k": float(cte)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_run(sys.argv[1])) + "\n")
