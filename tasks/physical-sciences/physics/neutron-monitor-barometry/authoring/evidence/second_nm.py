"""Second, separately written reduction of a neutron-monitor station (same call as the reference).

    python second_nm.py <station directory> <result.json>

- Monitor pressure from an airport office: its sea-level value back to the office's station pressure (QNH with the
  standard atmosphere, QFF through the column at the office's air temperature, as the notes define them), then the
  hydrostatic equation integrated exactly through a layer whose temperature falls linearly at 0.0065 K/m
  from the office's air temperature, p = p_a (T(h) / T_a)^(g / (R_d 0.0065)) (the mean-temperature hypsometric form
  of the notes agrees with it to a few parts in a million over these layers).  A hut barometer's reading is raised by
  c rho v^2 / 2 with rho from the reading and the station's air temperature.
- The station's checks first: the reference variation from the reference's own hourly rates, (R / 421500 - 1) 100,
  and a hut barometer's readings mapped onto the travelling standard by the line through the visit pairs (solved as
  a 2 x 2 normal-equation system).
- Counts: a Poisson generalised linear model with log link, ln mu = ln N0 - beta (p - p0) / 100 + a cos + b sin of the
  hour + offset ln(1 + gamma r / 100), fitted for fixed (gamma, c) by iteratively reweighted least squares to the
  maximum of the Poisson likelihood (the daily wave enters to first order; its square, a few 1e-6, is dropped).
- gamma and c: the profile likelihood over them maximised by Nelder-Mead from gamma = 0.8, c = 0 (c held at 0 where
  the pressure is an airport office's, whose barometer sees no wind).
- Clock and reduction of the readings: not taken from the sheet but from the data, allowing one change within the
  record.  A choice is a slide of the source's rows against the counts (-12 to +12 hours) and, for an office, QNH or
  QFF; a stretch of hours is scored by the median deviance per hour of a Poisson fit (gamma 0.8; for a hut a free
  rho v^2 / 2 term) of the stretch alone (the median, so that the few hours of an enhancement do not decide it).
  The choice at each end of the record is the most frequent winner of its five outermost 48-hour stretches (stepped
  by a day).  Where the two ends agree, the best single choice scored on the whole record holds for all of it.
  Where they differ, 72-hour stretches stepped by a day are scored, and the record splits halfway between the end of
  the last stretch the first choice wins and the start of the first stretch the second choice wins.  The split is
  then refined by one Poisson fit of the 240 hours around it (ln mu linear in the pressure, the reference variation,
  the daily wave and, for a hut, rho v^2 / 2): split hours within 36 of it (every second hour, then the neighbours
  of the best) and the clocks next to each side's choice, the plan with the smallest mean deviance per hour (each
  hour capped at 50) winning.  The refined split is kept only when its mean deviance per hour over the whole record
  (one fit) is at least 3 % below that of the best single choice; otherwise the single choice holds.
  Count hours whose source row the two sides would both claim (a clock reset) are dropped.
- Counters: the 18 counters' sum drops by 1/18 for every counter out.  Under the plan, a Poisson fit of ln N on the
  pressure, the reference, the daily wave (and a hut's wind term) with one intercept per 72 hours gives the level of
  ln N hour by hour; hours more than 5 sigma above the median of the 25 around them are set aside (enhancements);
  the level is cut where the mean of the 12 kept hours after an hour less that of the 12 before has its largest
  extremes beyond 3.5 % (at least a day apart); each stretch counts 18 exp(mean - mean of the highest stretch two days long or
  more) counters, rounded, at most 18.  Plan and counters are found three times (each plan with the counters before
  it, their share ln(k / 18) in its fits), then the plan once more; the fits below carry the share as an offset.
- Enhancements: after the fit, hours whose Pearson residual exceeds +5 are taken out with the 3 hours before and the
  24 after each, and the fit is made again (twice at most, never leaving out more than a fifth of the record).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

GRAV, R_AIR, GAMMA_T = 9.80665, 287.05, 0.0065


def station_pressure(sheet, reading, temp_c, kind):
    src = sheet["pressure_source"]
    if src["kind"].startswith("airport"):
        h, ha = sheet["altitude_m"], src["elevation_m"]
        ta = temp_c + 273.15
        if kind == "QFF":
            pa = reading / np.exp(GRAV * ha / (R_AIR * (ta + GAMMA_T * ha / 2.0)))
        else:
            pa = reading * (1.0 - 2.25577e-5 * ha) ** 5.25588
        return pa * ((ta - GAMMA_T * (h - ha)) / ta) ** (GRAV / (R_AIR * GAMMA_T)), False
    return reading, True


def glm(X, n, offset, iters=60):
    """Poisson log-link fit by IRLS; returns the coefficients and the log-likelihood (without the n! term)."""
    beta = np.linalg.lstsq(X, np.log(n) - offset, rcond=None)[0]
    for _ in range(iters):
        eta = X @ beta + offset
        mu = np.exp(eta)
        z = eta - offset + (n - mu) / mu
        wts = mu
        new = np.linalg.solve(X.T @ (X * wts[:, None]), X.T @ (wts * z))
        if np.max(np.abs(new - beta)) < 1e-12:
            beta = new
            break
        beta = new
    eta = X @ beta + offset
    return beta, float(np.sum(n * eta - np.exp(eta)))


def slide(rows, k):
    """Counts of row i paired with the source readings of row i + k (rows beyond the record dropped)."""
    hour, n, reading, temp_c, wind, ref = rows.T
    idx = np.arange(hour.size)
    ok = (idx + k >= 0) & (idx + k < hour.size)
    j = idx[ok] + k
    return hour[ok], n[ok], reading[j], temp_c[j], wind[j], ref[ok]


def reduce_station(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "station.json").read_text())
    rows = np.loadtxt(folder / "hourly.csv", delimiter=",", skiprows=1)
    refr = np.genfromtxt(folder / "reference.csv", delimiter=",", names=True)
    by_hour = {int(h): v for h, v in zip(refr["hour"], refr["reference_rate_per_h"])}
    rows[:, 5] = np.array([100.0 * (by_hour[int(h)] / 421500.0 - 1.0) for h in rows[:, 0]])
    if (folder / "barometer_checks.csv").exists():
        vis = np.genfromtxt(folder / "barometer_checks.csv", delimiter=",", names=True)
        xr, ys = vis["hut_barometer_hpa"], vis["travelling_standard_hpa"]
        nrm = np.array([[xr.size, xr.sum()], [xr.sum(), (xr * xr).sum()]])
        a0, a1 = np.linalg.solve(nrm, np.array([ys.sum(), (xr * ys).sum()]))
        rows[:, 2] = a0 + a1 * rows[:, 2]
    p0 = sheet["standard_pressure_hpa"]

    size = rows.shape[0]
    kinds = ("QNH", "QFF") if sheet["pressure_source"]["kind"].startswith("airport") else ("QNH",)
    choices = [(k, kind) for k in range(-12, 13) for kind in kinds]
    pressure = {}
    for k, kind in choices:
        idx = np.arange(size)
        idx = idx[(idx + k >= 0) & (idx + k < size)]
        pressure[(k, kind)] = (idx, station_pressure(sheet, rows[idx + k, 2], rows[idx + k, 3], kind)[0])

    def build(plan):
        hours = np.array([h for h in range(size) if 0 <= h + plan[h][0] < size])
        used = hours + np.array([plan[h][0] for h in hours])
        once = np.bincount(used, minlength=size)[used] == 1             # rows claimed twice across a reset go
        hours, used = hours[once], used[once]
        p = np.empty(hours.size)
        for kind in kinds:
            here = np.array([plan[h][1] == kind for h in hours], bool)
            if here.any():
                p[here] = station_pressure(sheet, rows[used[here], 2], rows[used[here], 3], kind)[0]
        return hours, used, p

    def plan_for(offset):
        """Clock and reduction of every count hour (one change allowed), with the counters' share ln(k / 18) of every
        count hour taken into the fits."""
        def stretch(choice, lo, hi):
            """Median deviance per hour of the Poisson fit of count hours lo..hi-1 under one choice, times their number."""
            idx, p = pressure[choice]
            sel = (idx >= lo) & (idx < hi)
            if sel.sum() < 24:
                return np.inf
            hour, n, ref = rows[idx[sel], 0], rows[idx[sel], 1], rows[idx[sel], 5]
            ang = 2 * np.pi * hour / 24.0
            cols = [np.ones_like(hour), -(p[sel] - p0) / 100.0, np.cos(ang), np.sin(ang)]
            if not sheet["pressure_source"]["kind"].startswith("airport"):       # a hut's wind term, free
                j = idx[sel] + choice[0]
                cols.append(0.5 * (rows[j, 2] * 100.0 / (R_AIR * (rows[j, 3] + 273.15))) * rows[j, 4] ** 2 / 100.0)
            X = np.column_stack(cols)
            off = np.log(np.clip(1.0 + 0.8 * ref / 100.0, 1e-9, None)) + offset[idx[sel]]
            coef, _ = glm(X, n, off)
            mu = np.exp(X @ coef + off)
            return float(np.median(2 * (n * np.log(n / mu) - (n - mu)))) * sel.sum()

        def best(lo, hi):
            return min(((stretch(c, lo, hi), c) for c in choices), key=lambda t: t[0])

        whole, single = best(0, size)
        # the choice at each end of the record: the most frequent winner of the first and of the last five 48-hour stretches
        edge = lambda starts: max(set(w := [best(a, a + 48)[1] for a in starts]), key=w.count)  # noqa: E731
        left, right = edge((0, 24, 48, 72, 96)), edge((size - 144, size - 120, size - 96, size - 72, size - 48))
        if left != right:
            # the split between the last stretch still won by the left choice and the first won by the right one
            wins = [(a, best(a, a + 72)[1]) for a in range(0, size - 71, 24)]
            lo = max([a + 72 for a, c in wins if c == left and a + 72 <= next((b for b, d in wins if d == right), size)],
                     default=72)
            hi = min([a for a, c in wins if c == right and a >= lo - 72], default=size - 72)
            coarse = (None, (lo + hi) // 2)

        def joint(plan, lo=0, hi=None):
            """Mean deviance per hour of one Poisson fit of hours lo..hi-1 under a plan (ln mu linear in the pressure, the
            reference variation, the daily wave and, for a hut, rho v^2 / 2), each hour's deviance capped at 50 (an
            enhancement weighs alike under every plan, a misplaced hour still counts)."""
            hours, used, p = build(plan)
            sel = (hours >= lo) & (hours < (size if hi is None else hi))
            hours, used, p = hours[sel], used[sel], p[sel]
            hour, n, ref = rows[hours, 0], rows[hours, 1], rows[hours, 5]
            ang = 2 * np.pi * hour / 24.0
            cols = [np.ones_like(hour), -(p - p0) / 100.0, ref / 100.0, np.cos(ang), np.sin(ang)]
            if not sheet["pressure_source"]["kind"].startswith("airport"):       # a hut's wind term, free
                cols.append(0.5 * (rows[used, 2] * 100.0 / (R_AIR * (rows[used, 3] + 273.15))) * rows[used, 4] ** 2 / 100.0)
            X = np.column_stack(cols)
            coef, _ = glm(X, n, offset[hours])
            mu = np.exp(X @ coef + offset[hours])
            return float(np.mean(np.minimum(2 * (n * np.log(n / mu) - (n - mu)), 50.0)))

        if left != right:
            at = coarse[1]
            # the split hour and the neighbouring clocks of each side refined by one fit of the whole record
            near = lambda c: [(k, c[1]) for k in (c[0] - 1, c[0], c[0] + 1) if -12 <= k <= 12]  # noqa: E731
            lo, hi = max(at - 120, 0), min(at + 120, size)        # one fit of the stretch around the split
            trials = [(joint([lc if h < hh else rc for h in range(size)], lo, hi), hh, lc, rc)
                      for lc in near(left) for rc in near(right) if lc != rc
                      for hh in range(max(at - 36, 24), min(at + 37, size - 23), 2)]
            _, at, left, right = min(trials, key=lambda t: t[0])
            _, at, left, right = min(((joint([left if h < hh else right for h in range(size)], lo, hi), hh, left, right)
                                      for hh in (at - 1, at, at + 1)), key=lambda t: t[0])
            plan = [left if h < at else right for h in range(size)]
            if joint(plan) > 0.97 * joint([single] * size):      # a split that hardly helps (QNH and QFF alike) is dropped
                plan = [single] * size
        else:
            plan = [single] * size
        return plan

    def counters(plan):
        """ln(k / 18) for every count hour: a Poisson fit of ln N on the pressure, the reference, the daily wave (and
        a hut's wind term) with one intercept per 72 hours gives the level of ln N hour by hour; hours more than 5
        sigma above the median of the 25 around them (an enhancement) are set aside; the
        level where the mean of the 12 kept hours after an hour less that of the 12 before has a local extreme beyond
        3.5 % (steps at least a day apart); each stretch's k is 18 exp(mean - mean of the highest stretch two days long or more), rounded,
        at most 18."""
        hours, used, p = build(plan)
        hour, n, ref = rows[hours, 0], rows[hours, 1], rows[hours, 5]
        ang = 2 * np.pi * hour / 24.0
        cols = [-(p - p0) / 100.0, ref / 100.0, np.cos(ang), np.sin(ang)]
        if not sheet["pressure_source"]["kind"].startswith("airport"):
            cols.append(0.5 * (rows[used, 2] * 100.0 / (R_AIR * (rows[used, 3] + 273.15))) * rows[used, 4] ** 2 / 100.0)
        k_cols = len(cols)
        cols += [(hours // 72 == b).astype(float) for b in range(size // 72 + 1)]
        X = np.column_stack(cols)
        keepc = X.any(axis=0)
        coef, _ = glm(X[:, keepc], n, np.zeros_like(hour))
        full = np.zeros(X.shape[1])
        full[keepc] = coef
        level = np.log(n) - X[:, :k_cols] @ full[:k_cols]
        sig = 1 / np.sqrt(n)
        med = np.array([np.median(level[max(i - 12, 0):i + 13]) for i in range(level.size)])
        ok = (level - med) / sig <= 5.0
        h, x = hours[ok], level[ok]
        # two-window detector: the mean of the 12 kept hours after each hour less that of the 12 before; its local
        # extremes beyond 3.5 % in size, at least a day apart, are the steps
        jump = np.zeros(h.size)
        for i in range(12, h.size - 11):
            jump[i] = np.mean(x[i:i + 12]) - np.mean(x[i - 12:i])
        cuts = []
        for i in np.argsort(-np.abs(jump)):
            if abs(jump[i]) < 0.035:
                break
            if all(abs(h[i] - h[c]) >= 24 for c in cuts):
                cuts.append(int(i))
        edges = [0] + sorted(cuts) + [h.size]
        segs = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
        means = [float(np.mean(x[a:b])) for a, b in segs]
        top = max((m for m, (a, b) in zip(means, segs) if h[b - 1] - h[a] >= 48), default=max(means))
        out = np.zeros(size)
        for (a, b), m in zip(segs, means):
            k = min(max(int(round(18 * np.exp(m - top))), 1), 18)
            lo = 0 if a == 0 else h[a]
            hi = size if b == h.size else h[b]
            out[lo:hi] = np.log(k / 18.0)
        return out

    offset = np.zeros(size)
    for _ in range(3):
        plan = plan_for(offset)
        offset = counters(plan)
    plan = plan_for(offset)
    hour_list, row_list, p_list = [], [], []
    for h in range(size):
        k, kind = plan[h]
        if 0 <= h + k < size:
            hour_list.append(h)
            row_list.append(h + k)
    rows_used = np.array(row_list)
    once = np.bincount(rows_used, minlength=size)[rows_used] == 1          # rows claimed twice across a reset go
    hours = np.array(hour_list)[once]
    rows_used = rows_used[once]
    p = np.empty(hours.size)
    for kind in kinds:
        here = np.array([plan[h][1] == kind for h in hours], bool)
        if here.any():
            p[here] = station_pressure(sheet, rows[rows_used[here], 2], rows[rows_used[here], 3], kind)[0]
    hut = not sheet["pressure_source"]["kind"].startswith("airport")
    hour, n, ref = rows[hours, 0], rows[hours, 1], rows[hours, 5]
    reading, temp_c, wind = rows[rows_used, 2], rows[rows_used, 3], rows[rows_used, 4]
    dyn = 0.5 * (reading * 100.0 / (R_AIR * (temp_c + 273.15))) * wind ** 2 / 100.0      # hPa per unit c
    ang = 2 * np.pi * hour / 24.0
    use = np.ones(hour.size, bool)

    def fit(g, c, sel):
        ps = p + c * dyn
        X = np.column_stack([np.ones_like(hour), -(ps - p0) / 100.0, np.cos(ang), np.sin(ang)])
        off = np.log(np.clip(1.0 + g * ref / 100.0, 1e-9, None)) + offset[hours]
        coef, ll = glm(X[sel], n[sel], off[sel])
        return coef, ll, np.exp(X @ coef + off)

    def solve(sel):
        if hut:
            res = minimize(lambda z: -fit(z[0], z[1], sel)[1], [0.8, 0.0], method="Nelder-Mead",
                           options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 4000})
            return res.x[0], res.x[1]
        res = minimize(lambda z: -fit(z[0], 0.0, sel)[1], [0.8], method="Nelder-Mead",
                       options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 2000})
        return res.x[0], 0.0

    g, c = solve(use)
    for _ in range(2):
        _, _, mu = fit(g, c, use)
        hot = np.flatnonzero((n - mu) / np.sqrt(mu) > 5.0)
        if hot.size == 0:
            break
        new = use.copy()
        for i in hot:
            new[max(i - 3, 0):i + 25] = False
        if np.array_equal(new, use) or new.sum() < 0.8 * hour.size:
            break
        use = new
        g, c = solve(use)
    coef, _, _ = fit(g, c, use)
    return {"barometric_coefficient_pct_per_hpa": float(coef[1]),
            "rate_at_standard_pressure_per_h": float(np.exp(coef[0]))}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_station(sys.argv[1])))
