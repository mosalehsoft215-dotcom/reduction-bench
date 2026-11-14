"""Second reduction, written apart from the reference.

    python3 second_notch.py <station folder> <answer.json>

- Barometer lag: both pressure records on a 5-minute grid; for each trial lag (whole minutes over +-300) the well
  pressure minus the shifted air pressure is high-passed (35-minute running mean removed) and its variance taken in
  12-hour windows stepped by 6 hours; from cumulative sums of that table a step function of one, two or three
  segments (at least two windows each) is fitted, each segment at the lag of least summed variance, a further segment
  kept when it cuts the total by more than 3 % and every step falls by 20 minutes or more. Where two segments meet, the barometer clock of the change
  (the clock stood still through a power cut) is the one, in 5-minute then 1-minute steps, that minimises the
  high-passed variance over the day around it; the air pressure is interpolated across the cut.
- Gaugings: each conductivity excess is fitted (scipy least_squares) with the stirred-pool passage, a Gaussian of
  mean mu and spread sigma convolved with an exponential of time constant tau, scaled by the salt mass over the flow;
  the fitted model is integrated to infinity, so a record cut short is completed by the model itself.
- Rating: Q(t) = C s(t) + a l(t), s = (h1 + 0.85 mm)^2.5 x Villemonte factor, l = sqrt(h1 - pool level, or + sill
  depth while the pool is below the sill); C and a minimise the squared log misfit of the five gaugings, each the
  salt-share harmonic mean of Q(t) over its model passage (Nelder-Mead from three starts, a kept >= 0).
- Leak step: at each of the four highest maxima of the compensated head (12 h apart) the leak may step from a1 to a2
  (C shared, Nelder-Mead); the best step is kept when it brings the misfit under 0.3 of the one-leak misfit and moves
  the leak by more than 0.3 L/s per sqrt(m).
- Stilling well: the channel level behind the well from the intake law, with k and the pool line from every
  stretch the pool stands above the crest (function channel below).
- Week: discharge C s(t) + a l(t) on the logged minutes, trapezoid volume, 15-minute peak on a 5-minute clock.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize


def rows(path, cols):
    with open(path, newline="") as fh:
        return np.array([[float(r[c]) for c in cols] for r in csv.DictReader(fh)])


def passage(t_s, mu, sigma, tau):
    grid = np.arange(-600.0, t_s.max() + 12 * tau + 1.0, 1.0)
    pulse = np.exp(-0.5 * ((grid - mu) / sigma) ** 2)
    pulse /= pulse.sum()
    decay = np.exp(-np.arange(0.0, 12 * tau, 1.0) / tau)
    decay /= decay.sum()
    dens = np.convolve(pulse, decay)[: len(grid)]
    return grid, dens


def shape_of(h, h2):
    f = np.where((h2 > 0) & (h > 0), np.clip(1 - np.clip(h2 / np.where(h > 0, h, 1), 0, 1) ** 2.5, 0, None) ** 0.385, 1.0)
    return np.clip(h + 0.00085, 0, None) ** 2.5 * f


def level_of(s, h2):
    grid = np.linspace(0.0, 0.8, 1601)
    out = np.empty_like(s)
    for i, (si, xi) in enumerate(zip(s, h2)):
        sg = shape_of(grid, np.full(grid.size, xi))
        j = int(np.clip(np.searchsorted(sg, si), 1, grid.size - 1))
        out[i] = grid[j - 1] + (si - sg[j - 1]) / max(sg[j] - sg[j - 1], 1e-15) * (grid[j] - grid[j - 1])
    return out


def channel(m, hw, h2):
    """the channel level behind the stilling well (dh_w/dt = k sign(h_c - h_w) sqrt|h_c - h_w|). Rates are 3-point
    averages of centred differences. Every stretch with the pool at least 5 mm over the crest for 40 minutes or more is
    scanned over 40 values of k (0.002-0.07): the channel level from the law, its weir shape against the pool level,
    and the scatter about a straight line; the k of least scatter is kept when it cuts the scatter of a clean pipe
    (0.07) by 40 % or more, then k and the line are refined together (Nelder-Mead) so that the well, stepped by Heun's
    method in 0.5-minute steps from the line's channel level, follows the logged well level. That k holds from the lowest well level of the 18 h before the stretch on; inside the stretch
    the line gives the channel level, elsewhere the law with k (0.07 where no stretch tells)"""
    rate = np.gradient(hw, m)
    rate = np.convolve(np.concatenate([[rate[0]], rate, [rate[-1]]]), np.ones(3) / 3.0, mode="valid")
    up = h2 > 0.005
    edges = np.flatnonzero(np.diff(np.concatenate([[0], up.astype(int), [0]])))
    runs = [(a, b) for a, b in zip(edges[::2], edges[1::2]) if m[b - 1] - m[a] >= 40.0]
    k = np.full(m.size, 0.07)
    lines = []
    for a, b in runs:
        x, y, r = h2[a:b], hw[a:b], rate[a:b]

        def scatter(kk):
            s = shape_of(y + np.sign(r) * (r / kk) ** 2, x)
            z = np.polyfit(x, s, 1)
            return float(np.sum((s - np.polyval(z, x)) ** 2)), z
        ks = np.geomspace(0.002, 0.07, 40)
        vals = [scatter(kk)[0] for kk in ks]
        j = int(np.argmin(vals))
        if vals[j] > 0.6 * vals[-1]:
            continue
        kk, z = ks[j], scatter(ks[j])[1]
        tm = m[a:b]

        def sse(v):                                          # Heun steps of 0.5 min driven by the line's channel level
            kq, zz = math.exp(v[0]), v[1:]
            hcl = level_of(np.polyval(zz, x), x)
            w_ = y[0]
            err = 0.0
            for i in range(1, y.size):
                n_ = max(int(round((tm[i] - tm[i - 1]) / 0.5)), 1)
                dt = (tm[i] - tm[i - 1]) / n_
                for s_ in range(n_):
                    c1 = hcl[i - 1] + (hcl[i] - hcl[i - 1]) * s_ / n_
                    c2 = hcl[i - 1] + (hcl[i] - hcl[i - 1]) * (s_ + 1) / n_
                    f1 = kq * np.sign(c1 - w_) * math.sqrt(abs(c1 - w_))
                    wp = w_ + dt * f1
                    f2 = kq * np.sign(c2 - wp) * math.sqrt(abs(c2 - wp))
                    w_ = w_ + 0.5 * dt * (f1 + f2)
                err += (w_ - y[i]) ** 2
            return err
        best = minimize(sse, np.concatenate([[math.log(kk)], z]), method="Nelder-Mead",
                        options={"xatol": 1e-7, "fatol": 1e-14, "maxiter": 1500})
        kk, z = math.exp(best.x[0]), best.x[1:]
        if not 0.0 < kk < 0.07:
            continue
        start = (m >= m[a] - 18 * 60.0) & (m <= m[a])
        t0 = m[start][int(np.argmin(hw[start]))] if start.any() else m[a]
        k[m >= t0] = kk
        lines.append((a, b, z))
    hc = hw + np.sign(rate) * (rate / k) ** 2
    for a, b, z in lines:
        hc[a:b] = level_of(np.polyval(z, h2[a:b]), h2[a:b])
    return hc


def reduce(folder):
    folder = Path(folder)
    sh = json.loads((folder / "station_sheet.json").read_text())
    w = rows(folder / "well_log.csv", ("minute", "p_kPa", "tailwater_mm", "water_C"))
    b = rows(folder / "barometer_log.csv", ("minute", "p_air_hPa"))
    b[:, 1] /= 10.0
    grid = np.arange(w[0, 0], w[-1, 0] + 1.0, 5.0)
    pw = np.interp(grid, w[:, 0], w[:, 1])
    hp = lambda x: np.r_[0.0, np.diff(x)]  # noqa: E731
    lags = np.arange(-300, 301)
    win = [(i, i + 144) for i in range(4, grid.size - 148, 72)]
    var = np.array([[float(np.var(r[a:z])) for a, z in win]
                    for r in (hp(pw - np.interp(grid + lag, b[:, 0], b[:, 1])) for lag in lags)])
    cum = np.concatenate([np.zeros((lags.size, 1)), np.cumsum(var, axis=1)], axis=1)
    nw = len(win)
    cost = lambda a, z: float(np.min(cum[:, z] - cum[:, a]))           # noqa: E731
    lag_of = lambda a, z: float(lags[np.argmin(cum[:, z] - cum[:, a])])  # noqa: E731
    plans = {1: (cost(0, nw), [0, nw])}
    two = [(cost(0, i) + cost(i, nw), [0, i, nw]) for i in range(2, nw - 1)]
    three = [(cost(0, i) + cost(i, j) + cost(j, nw), [0, i, j, nw]) for i in range(2, nw - 3) for j in range(i + 2, nw - 1)]
    plans[2], plans[3] = min(two), min(three)
    pick = plans[1]
    for n in (2, 3):
        edges = plans[n][1]
        lg = [lag_of(edges[q], edges[q + 1]) for q in range(n)]
        if plans[n][0] < 0.97 * pick[0] and all(lg[q] - lg[q + 1] >= 20 for q in range(n - 1)):
            pick = plans[n]
    edges = pick[1]
    merged = [list(range(edges[q], edges[q + 1])) for q in range(len(edges) - 1)]
    seg_lag = [lag_of(edges[q], edges[q + 1]) for q in range(len(edges) - 1)]

    def air(cuts):
        clock = b[:, 0]
        true = clock - seg_lag[0]
        for (c, lg) in cuts:
            true = np.where(clock >= c, clock - lg, true)
        o = np.argsort(true)
        return np.interp(w[:, 0], true[o], b[o, 1])

    cuts = []
    for k in range(1, len(seg_lag)):
        t_mid = grid[win[merged[k][0]][0]]
        near = (w[:, 0] > t_mid - 900.0) & (w[:, 0] < t_mid + 900.0)

        def score(c):
            r = w[near, 1] - air(cuts + [(c, seg_lag[k])])[near]
            return float(np.var(r - np.convolve(r, np.ones(7) / 7.0, mode="same")))
        base = t_mid + seg_lag[k - 1]
        c5 = min(np.arange(base - 720.0, base + 720.0, 5.0), key=score)
        cuts.append((float(min(np.arange(c5 - 5.0, c5 + 5.1, 1.0), key=score)), seg_lag[k]))
    t = w[:, 3]
    rho = 999.97495 * (1 - (t - 3.983035) ** 2 * (t + 301.797) / (522528.9 * (t + 69.34881)))
    h1 = (w[:, 1] - air(cuts)) * 1000.0 / (rho * 9.81) - sh["well_logger"]["sensor_below_crest_m"]
    h2 = w[:, 2] / 1000.0
    raw = channel(w[:, 0], h1, h2)

    h = raw
    drowned = (h2 > 0) & (h > 0)
    factor = np.where(drowned, np.clip(1 - np.clip(h2 / np.where(h > 0, h, 1), 0, 1) ** 2.5, 0, None) ** 0.385, 1.0)
    shape = np.clip(h + 0.00085, 0, None) ** 2.5 * factor
    sill = sh["plate"]["sill_below_crest_m"]
    leak = np.sqrt(np.clip(h - np.maximum(h2, -sill), 0, None))

    passes = []
    for g in sh["gaugings"]:
        tr = rows(folder / g["file"], ("minute", "ec_uS_cm"))
        t0 = g["injection_minute"]
        bg = float(tr[tr[:, 0] < t0, 1].mean())
        sel = tr[:, 0] >= t0
        ts = (tr[sel, 0] - t0) * 60.0
        conc = (tr[sel, 1] - bg) / g["meter_uS_cm_per_mg_L"] / 1000.0

        def model(x):
            area, mu, sigma, tau = x
            gr, dens = passage(ts, mu, sigma, tau)
            per_s = np.interp(ts, gr, dens)
            return area * per_s

        x0 = [float(conc.sum() * 10.0), float(ts[np.argmax(conc)]), 50.0, 200.0]
        fit = least_squares(lambda x: model(x) - conc, x0, bounds=([1e-6, 0.0, 5.0, 20.0], [1e6, 3000.0, 400.0, 2000.0]))
        area, mu, sigma, tau = fit.x
        q = g["salt_kg"] / area
        gr, dens = passage(ts, mu, sigma, tau)
        passes.append((q, t0 + gr / 60.0, dens))

    at = [(q, np.interp(mins, w[:, 0], shape), np.interp(mins, w[:, 0], leak), dens / dens.sum()) for q, mins, dens in passes]

    def misfit(x):
        c, a = x[0], abs(x[1]) / 1000.0
        return sum((math.log(q) + math.log(float(np.sum(d / np.maximum(c * s + a * l, 1e-12))))) ** 2 for q, s, l, d in at)

    c0 = float(np.median([q * float(np.sum(d / np.maximum(s, 1e-12))) for q, s, l, d in at]))
    fits = [minimize(misfit, [c0 * f, a0], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-12, "maxiter": 4000})
            for f, a0 in ((1.0, 0.0), (0.95, 1.0), (0.9, 2.5))]
    best = min(fits, key=lambda r: r.fun)
    c, a = float(best.x[0]), abs(float(best.x[1])) / 1000.0
    qw = c * shape + a * leak
    order = np.argsort(h1)[::-1]
    peaks = []
    for i in order:
        if 0 < i < h1.size - 1 and all(abs(w[i, 0] - m) > 720.0 for m in peaks):
            peaks.append(float(w[i, 0]))
        if len(peaks) == 4:
            break
    steps = []
    for m in peaks:
        after = [float(np.sum(d * (mins >= m))) for q, mins, d in passes]

        def misfit2(x, after=after):
            c2, a1, a2 = x[0], abs(x[1]) / 1000.0, abs(x[2]) / 1000.0
            return sum((math.log(q) + math.log(float(np.sum(d / np.maximum(c2 * s + (a1 + (a2 - a1) * f) * l, 1e-12))))) ** 2
                       for (q, s, l, d), f in zip(at, after))
        fit2 = min((minimize(misfit2, [c, a1, a2], method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-12,
                                                                                   "maxiter": 4000})
                    for a1, a2 in ((a * 1000.0, a * 1000.0), (2.0, 0.2), (0.2, 2.0))), key=lambda r: r.fun)
        steps.append((fit2.fun, fit2.x, m))
    fun2, x2, m = min(steps, key=lambda r: r[0])
    if fun2 < 0.3 * best.fun and abs(abs(x2[2]) - abs(x2[1])) > 0.3:
        c, a1, a2 = float(x2[0]), abs(float(x2[1])) / 1000.0, abs(float(x2[2])) / 1000.0
        qw = c * shape + np.where(w[:, 0] >= m, a2, a1) * leak
    secs = w[:, 0] * 60.0
    volume = float(np.sum((qw[1:] + qw[:-1]) / 2 * np.diff(secs)) + 150.0 * (qw[0] + qw[-1]))
    clock = np.arange(w[0, 0], w[-1, 0] + 1, 5.0)
    peak = float(np.max(np.convolve(np.interp(clock, w[:, 0], qw), np.ones(3) / 3, mode="valid")))
    return {"volume_m3": volume, "peak_15min_L_s": peak * 1000.0}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)
