"""Reduces 24 days of hourly neutron-monitor data to the barometric coefficient and the rate at the station's standard
pressure for the undisturbed galactic flux.

    python3 nm_reduce.py <station directory> <result.json>

Before everything, the station's own checks replace two delivered values.  The reference variation is computed from
the reference monitor's own hourly rates (reference.csv), r = (R / R0 - 1) x 100 with the network's undisturbed rate
R0 = 421500 per hour, not taken from the data service's column (which is computed against the base the service held).
A hut barometer's readings are carried to the travelling standard by the least-squares line standard = A + B reading
through the visits in barometer_checks.csv (the barometer is linear in pressure over its range).
0. Clock and reduction, read from the data hour by hour (the register's entries in the sheet are not used, and a
   source can change either within the record).  Candidates: the source's rows slid against the counts by -12 to +12
   hours and, for an office, its readings taken as QNH and as QFF.  (a) In 96-hour windows stepped by a day, every
   candidate is fitted by a linear count-weighted least-squares model (ln N on the pressure, the reference, the daily
   wave and, for a hut, rho v^2 / 2; hours beyond 5 sigma refitted out, three passes) and the one with the smallest
   mean squared standardised residual (each hour capped at 25, so that an enhancement weighs no more than a few
   ordinary hours) wins the window; a candidate whose slide leaves fewer than half the window's hours does not
   compete.  A window that holds a change, a counter's step or a large enhancement fits no candidate well,
   so only windows whose winner scores within three times the 20th percentile of all the windows' winning scores
   propose candidates, and every hour starts with the winner of the nearest such window.  (b) The model is fitted
   once for the whole record on the hours as assigned (the station's constants do not change), which gives every hour
   a capped squared residual under each winning candidate.  (c) Dynamic programming assigns one candidate to every
   hour with the smallest total of these residuals plus 60 per change of candidate (in this fit the level has one
   intercept per 48 hours, so that a counter's step, step 5, does not blind the comparison); (b) and (c) are repeated until
   the assignment stops changing (three times at most).  Every candidate other than the register's clock and
   reduction carries 0.15 per hour (in the windows' mean scores and in the hourly residuals): the register holds
   wherever the data cannot tell the candidates apart, as for an office a few metres from the monitor.  Where a
   change of clock makes two count hours claim the same source row (around a reset), both hours are left out.
1. Static pressure at the monitor.  An office's sea-level value goes back to its station pressure p_a, as QNH with the
   standard atmosphere, p_a = QNH (1 - 2.25577e-5 h_a)^5.25588, as QFF through the column at the office's air
   temperature, p_a = QFF exp(-g h_a / (R_d (T_a + 0.0065 h_a / 2))); the monitor's pressure follows from the
   hypsometric equation p = p_a exp(-g (h - h_a) / (R_d T_m)), T_m the mean of T_a and T_a - 0.0065 (h - h_a).
   A hut barometer reads the static pressure less c rho v^2 / 2; rho = p / (R_d T) from the station's air temperature,
   and the hut's c is fitted with the rest, between 0 and 1.5.
2. Model of the hourly counts: ln N = ln N0 - beta (p - p0) / 100 + ln(1 + gamma r / 100) + ln(1 + a cos + b sin of
   the hour), r the reference monitor's variation (0 % = undisturbed flux) and gamma the station's coupling to it.
3. Fit: Poisson-weighted least squares on ln N (weights sqrt(counts)), trust-region from a linear start
   (the linear model in ln N0, beta, gamma, a, b with c = 0); the coupling and the wind coefficient are fitted jointly
   with the barometric coefficient because storms bring both low pressure and wind, and Forbush decreases may fall
   on the lows.  N0 = exp(ln N0): the rate at p0, r = 0, no wind and the daily mean of the wave.
4. Enhancements: hours whose counts stand more than 6 sigma above the fit are left out with the 2 hours before and
   the 18 after each, and the fit repeated (up to four times; never more than a fifth of the record).
5. Counters: the counts are the sum of 18 counters, and one that fails takes its 1/18 with it until it is repaired.
   Under the clock and reduction of step 0, a linear count-weighted fit whose level has one intercept per 48 hours
   (so that the pressure, reference and wave terms are fitted through any step) gives the level of ln N hour by
   hour; hours more than 6 sigma above the running median of the 25 around them (an enhancement) are set aside, the
   rest is cut into stretches of constant level by optimal partitioning (weighted squared deviations plus 500 per
   cut), neighbours less than 3.5 % apart merged (one counter is 5.6 %) and stretches under 12 hours joined to a
   neighbour, and each stretch's number of counters is 18 exp(level - highest level held a day or more) rounded
   (between 1 and 18).  Step 0 and this step are made twice, the second time step 0 with the counters found (their
   share taken out of ln N); step 0 is then made a third time with the final counters, and steps 1-4 fit the model
   with ln(counters / 18) as an offset.  N0 is the rate of all 18 counters.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

G, RD, LAPSE = 9.80665, 287.05, 0.0065
R0 = 421500.0                           # the reference monitor's undisturbed rate per hour (the network's base)


def monitor_pressure(sheet, p_log, t_c, reduction="QNH", temperature=True):
    """Static pressure at the monitor from the logged pressure; True when the source is the hut barometer."""
    src = sheet["pressure_source"]
    if not src["kind"].startswith("airport"):
        return p_log, True
    h, h_a = sheet["altitude_m"], src["elevation_m"]
    t_a = t_c + 273.15
    if reduction == "QFF":
        p_a = p_log * np.exp(-G * h_a / (RD * (t_a + LAPSE * h_a / 2)))
    else:
        p_a = p_log * (1 - 2.25577e-5 * h_a) ** 5.25588
    if not temperature:                                   # standard atmosphere from the office up to the monitor
        return p_a * ((1 - 2.25577e-5 * h) / (1 - 2.25577e-5 * h_a)) ** 5.25588, False
    t_m = t_a - LAPSE * (h - h_a) / 2
    return p_a * np.exp(-G * (h - h_a) / (RD * t_m)), False


def design(sheet, rows, k, red, temperature, offset=None):
    """Count hours u whose reading row u + k exists, with the columns of the linear model of ln N under clock k and
    reduction red (1, the pressure term, the reference, the daily wave and, for a hut, rho v^2 / 2), ln N and the
    count weights sqrt(N)."""
    p0 = sheet["standard_pressure_hpa"]
    t, n, p_log, t_c, v, r = rows.T
    u = np.arange(t.size)
    u = u[(u + k >= 0) & (u + k < t.size)]
    j = u + k
    p, hut = monitor_pressure(sheet, p_log[j], t_c[j], red, temperature)
    cols = [np.ones(u.size), -(p - p0) / 100, r[u] / 100, np.cos(2 * np.pi * t[u] / 24), np.sin(2 * np.pi * t[u] / 24)]
    if hut:
        cols.append(0.5 * (p * 100 / (RD * (t_c[j] + 273.15))) * v[j] ** 2 / 100)
    y = np.log(n[u]) - (0.0 if offset is None else offset[u])      # counts brought to all 18 counters
    return u, np.column_stack(cols), y, np.sqrt(n[u])


def robust(X, y, w):
    """Count-weighted least squares with hours beyond 5 sigma refitted out (three passes): coefficients."""
    keep = np.ones(y.size, bool)
    for _ in range(3):
        coef, *_ = np.linalg.lstsq(X[keep] * w[keep, None], y[keep] * w[keep], rcond=None)
        z = (y - X @ coef) * w
        keep = np.abs(z) < 5 * max(np.std(z[keep]), 1.0)
    return coef


def schedule(sheet, rows, temperature, segments=True, window=96, penalty=60.0, prior=0.15, offset=None):
    """Clock offset and reduction of every count hour (see step 0); segments=False: one candidate for the whole record,
    the one with the smallest capped mean squared residual over all of it."""
    office = sheet["pressure_source"]["kind"].startswith("airport")
    cands = [(k, red) for red in (("QNH", "QFF") if office else ("QNH",)) for k in range(-12, 13)]
    size = rows.shape[0]
    built = {c: design(sheet, rows, c[0], c[1], temperature, offset) for c in cands}
    src = sheet["pressure_source"]
    register = (int(src.get("clock_utc_offset_h", 0)), src.get("reduction", "QNH"))
    doubt = {c: 0.0 if c == register else prior for c in cands}       # the register holds unless the data say otherwise

    def capped(c, lo, hi):
        u, X, y, w = built[c]
        sel = (u >= lo) & (u < hi)
        if sel.sum() < (hi - lo) // 2:
            return np.inf
        coef = robust(X[sel], y[sel], w[sel])
        return float(np.mean(np.minimum(((y[sel] - X[sel] @ coef) * w[sel]) ** 2, 25.0))) + doubt[c]

    if not segments:
        best = min(cands, key=lambda c: capped(c, 0, size))
        return np.full(size, best[0]), [best[1]] * size
    starts = list(range(0, max(size - window, 0) + 1, 24))
    scored = [min((capped(c, a, a + window), c) for c in cands) for a in starts]
    # a window that holds a change of clock or reduction, or a large enhancement, fits no candidate well: only windows
    # whose winner fits within three times the lowest fifth of the winning scores propose candidates and start the
    # assignment (with changes and counter steps, half the windows or more can fit badly)
    level = 3.0 * float(np.percentile([sc for sc, _ in scored], 20))
    good = [(a + window / 2, c) for a, (sc, c) in zip(starts, scored) if sc <= level]
    kept = sorted({c for _, c in good}, key=[c for _, c in good].index)
    centres = np.array([x for x, _ in good])
    path = np.array([kept.index(good[int(np.argmin(np.abs(centres - h)))][1]) for h in range(size)])
    blocks = lambda u: (u[:, None] // 48 == np.arange(size // 48 + 1)[None, :]).astype(float)  # noqa: E731
    for _ in range(3):
        # one set of station constants for the whole record, fitted on the hours as currently assigned; the level
        # gets one intercept per 48 hours, so that a counter's step does not blind the comparison of the candidates
        parts = []
        for m, c in enumerate(kept):
            u, X, y, w = built[c]
            sel = path[u] == m
            parts.append((np.c_[X[sel][:, 1:], blocks(u[sel])], y[sel], w[sel]))
        coef = robust(*(np.concatenate(q) for q in zip(*parts)))
        cost = np.full((size, len(kept)), 25.0)
        for m, c in enumerate(kept):
            u, X, y, w = built[c]
            cost[u, m] = np.minimum(((y - np.c_[X[:, 1:], blocks(u)] @ coef) * w) ** 2, 25.0) + doubt[c]
        total = cost[0].copy()
        back = np.zeros((size, len(kept)), int)
        for h in range(1, size):
            move = total.min() + penalty
            back[h] = np.where(total <= move, np.arange(len(kept)), int(np.argmin(total)))
            total = np.minimum(total, move) + cost[h]
        new = np.empty(size, int)
        new[-1] = int(np.argmin(total))
        for h in range(size - 1, 0, -1):
            new[h - 1] = back[h, new[h]]
        if np.array_equal(new, path):
            break
        path = new
    return np.array([kept[m][0] for m in path]), [kept[m][1] for m in path]


def counters(hours, level, w, size, penalty=500.0, step=0.035):
    """Number of the 18 counters at work in every count hour, from the level of ln N against the full model: the
    record is cut into stretches of constant level by optimal partitioning (weighted squared deviations from each
    stretch's mean plus `penalty` per cut), neighbouring stretches less than `step` apart in level are merged, a
    stretch shorter than 12 hours joins its neighbour, and a
    stretch's count is 18 exp(level - highest level) rounded (a share
    below half a counter is not a counter), the highest level being that of the highest stretch a day long or more
    (an enhancement's few hours are no level) and no stretch counting more than 18.  Returns ln(counters / 18) for every hour of the record."""
    # an enhancement is hours far above the level around them: hours more than 6 sigma above the running median of
    # the 25 hours around them are left out before the partition (a counter's step moves the median with it)
    med = np.array([np.median(level[max(i - 12, 0):i + 13]) for i in range(level.size)])
    calm = (level - med) * w <= 6.0
    hours, x, ww = hours[calm], level[calm], (w * w)[calm]
    m = x.size
    cw, cx, cxx = (np.r_[0.0, np.cumsum(a)] for a in (ww, ww * x, ww * x * x))
    best = np.zeros(m + 1)
    cut = np.zeros(m + 1, int)
    for j in range(1, m + 1):
        i = np.arange(j)
        sw, sx, sxx = cw[j] - cw[i], cx[j] - cx[i], cxx[j] - cxx[i]
        cost = best[i] + sxx - sx * sx / sw + penalty
        cut[j] = int(np.argmin(cost))
        best[j] = cost[cut[j]]
    bounds, j = [], m
    while j > 0:
        bounds.append((cut[j], j))
        j = cut[j]
    bounds.reverse()
    mean = lambda a, b: float(np.sum(ww[a:b] * x[a:b]) / np.sum(ww[a:b]))  # noqa: E731
    merged = [bounds[0]]                    # neighbours closer than `step` (half a counter is 2.8 %) are one stretch
    for a, b in bounds[1:]:
        if abs(mean(a, b) - mean(*merged[-1])) < step:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    bounds = []                             # a stretch shorter than 12 hours is no counter's outage: it joins its
    for a, b in merged:                     # neighbour before it (or after it, at the start)
        if bounds and hours[b - 1] - hours[a] < 12:
            bounds[-1] = (bounds[-1][0], b)
        elif bounds and hours[bounds[-1][1] - 1] - hours[bounds[-1][0]] < 12:
            bounds[-1] = (bounds[-1][0], b)
        else:
            bounds.append((a, b))
    means = [mean(a, b) for a, b in bounds]
    # all 18 at work: the highest level held for a day at least (an enhancement's few hours are not a level)
    top = max((mu for mu, (a, b) in zip(means, bounds) if hours[b - 1] - hours[a] >= 24), default=max(means))
    tubes = [min(max(int(round(18 * np.exp(mu - top))), 1), 18) for mu in means]
    out = np.zeros(size)
    marks = [hours[a] for a, _ in bounds] + [size]
    for (a, b), k, lo, hi in zip(bounds, tubes, marks[:-1], marks[1:]):
        out[(np.arange(size) >= (0 if a == 0 else lo)) & (np.arange(size) < hi)] = np.log(k / 18.0)
    return out


def checked_rows(folder, rows, reference="rates", barometer="checked"):
    """The hourly rows with the reference variation from the reference's own rates and a hut barometer's readings
    carried to the travelling standard (reference='delivered' / barometer='as read' keep the columns as they are)."""
    rows = rows.copy()
    if reference == "rates":
        ref = np.loadtxt(Path(folder) / "reference.csv", delimiter=",", skiprows=1, ndmin=2)
        rate = dict(zip(ref[:, 0].astype(int), ref[:, 1]))
        rows[:, 5] = [(rate[int(h)] / R0 - 1) * 100 for h in rows[:, 0]]
    chk = Path(folder) / "barometer_checks.csv"
    if barometer == "checked" and chk.exists():
        c = np.loadtxt(chk, delimiter=",", skiprows=1, ndmin=2)
        slope, inter = np.polyfit(c[:, 1], c[:, 2], 1)
        rows[:, 2] = inter + slope * rows[:, 2]
    return rows


def reduce_station(folder, temperature=True, wind=True, primary=True, coupling=True, qnh=True, sheet=False,
                   segments=True, enhancement=True, tubes=True, reference="rates", barometer="checked"):
    """Options other than the defaults are the shortcuts tabulated by authoring/evidence/bench.py (sheet=True: the
    register's clock and reduction for the whole record; segments=False: one clock and one reduction read from the
    data for the whole record; qnh=False: the delivered pressure as the monitor's, rows as listed; tubes=False: all 18
    counters taken to count throughout; reference='delivered': the service's variation column; barometer='as read':
    the hut barometer's readings as logged)."""
    folder = Path(folder)
    st = json.loads((folder / "station.json").read_text())
    rows = checked_rows(folder, np.loadtxt(folder / "hourly.csv", delimiter=",", skiprows=1), reference, barometer)
    size = rows.shape[0]
    offset = np.zeros(size)
    if tubes and qnh and not sheet:
        # the counters at work, from ln N against a linear fit whose level has one intercept per 48 hours (so that the
        # pressure, the reference and the wave are fitted through any counter's step), then the full fit
        for _ in range(2):                  # clock and reduction, then counters, then both again with them
            u, X, y, w = aligned(st, rows, temperature, offset)
            cols = np.c_[X[:, 1:], (u[:, None] // 48 == np.arange(size // 48 + 1)[None, :]).astype(float)]
            coef = robust(cols, y, w)
            level = y - X[:, 1:] @ coef[:X.shape[1] - 1]
            offset = counters(u, level, w, size)
    x, hours, level, w, keep = solve(st, rows, offset, temperature, wind, primary, coupling, qnh, sheet, segments,
                                     enhancement)
    ln0, beta = x[0], x[1]
    return {"barometric_coefficient_pct_per_hpa": float(beta), "rate_at_standard_pressure_per_h": float(np.exp(ln0))}


def aligned(st, rows, temperature, offset=None):
    """Count hours, linear design (1, pressure term, reference, daily wave, hut wind term), ln N and weights under the
    hour-by-hour clock and reduction of step 0."""
    size = rows.shape[0]
    ks, reds = schedule(st, rows, temperature, offset=offset)
    parts = [design(st, rows, k, red, temperature) for k, red in sorted(set(zip(ks.tolist(), reds)))]
    want = {(h, int(ks[h]), reds[h]) for h in range(size)}
    keys = sorted(set(zip(ks.tolist(), reds)))
    us, Xs, ys, ws = [], [], [], []
    for (k, red), (u, X, y, w) in zip(keys, parts):
        sel = np.array([(int(h), k, red) in want for h in u], bool)
        us.append(u[sel]); Xs.append(X[sel]); ys.append(y[sel]); ws.append(w[sel])
    u, X, y, w = (np.concatenate(q) for q in (us, Xs, ys, ws))
    j = u + ks[u]
    once = np.bincount(j, minlength=size)[j] == 1        # rows two count hours would claim (around a reset) go
    u, X, y, w = u[once], X[once], y[once], w[once]
    order = np.argsort(u)
    return u[order], X[order], y[order], w[order]


def solve(st, rows, offset, temperature, wind, primary, coupling, qnh, sheet, segments, enhancement):
    """One round: clock and reduction of every hour, then the fit of the model with the counters' share as an offset.
    Returns the parameters, the count hours used, ln N less the model without the counters' share, the weights and
    the hours kept (outside enhancements)."""
    src = st["pressure_source"]
    size = rows.shape[0]
    if sheet or not qnh:
        ks, reds = np.full(size, int(src.get("clock_utc_offset_h", 0)) if qnh else 0), [src.get("reduction", "QNH")] * size
    else:
        ks, reds = schedule(st, rows, temperature, segments, offset=offset)
    u = np.arange(size)
    j = u + ks
    ok = (j >= 0) & (j < size)
    # a clock reset leaves readings of some hours unrecorded and rows that two count hours would claim: those pairs go
    ok &= np.bincount(np.clip(j, 0, size - 1)[ok], minlength=size)[np.clip(j, 0, size - 1)] == 1
    u, j, reds = u[ok], j[ok], np.array(reds)[ok]
    t, n, r = rows[u, 0], rows[u, 1], rows[u, 5]
    p_log, t_c, v = rows[j, 2], rows[j, 3], rows[j, 4]
    p0 = st["standard_pressure_hpa"]
    hut = src["kind"].startswith("hut")
    if qnh:
        p = np.empty(u.size)
        for red in set(reds.tolist()):
            here = reds == red
            p[here] = monitor_pressure(st, p_log[here], t_c[here], red, temperature)[0]
    else:
        p = p_log
    q = 0.5 * (p * 100 / (RD * (t_c + 273.15))) * v ** 2 / 100        # rho v^2 / 2 in hPa
    use_wind = wind and hut
    r_in = r if primary else np.zeros_like(r)
    y, w = np.log(n) - offset[u], np.sqrt(n)
    cs, sn = np.cos(2 * np.pi * t / 24), np.sin(2 * np.pi * t / 24)

    def model(x, sel):
        ln0, beta, gam, a, b, c = x
        ps = p[sel] + (c * q[sel] if use_wind else 0.0)
        g = gam if coupling else 1.0
        return (ln0 - beta * (ps - p0) / 100 + np.log(np.clip(1 + g * r_in[sel] / 100, 1e-6, None))
                + np.log(np.clip(1 + a * cs[sel] + b * sn[sel], 1e-6, None)))

    def fit(sel, x0=None):
        if x0 is None:
            lin = np.column_stack([np.ones_like(t), -(p - p0) / 100, r_in / 100, cs, sn])[sel]
            start, *_ = np.linalg.lstsq(lin * w[sel, None], y[sel] * w[sel], rcond=None)
            x0 = np.r_[start[0], start[1], start[2] if coupling and primary else 1.0, start[3], start[4], 0.0]
        lo = np.r_[-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 0.0]          # the hut's c: 0 to about 1 (notes)
        hi = np.r_[np.inf, np.inf, np.inf, np.inf, np.inf, 1.5]
        return least_squares(lambda x: (model(x, sel) - y[sel]) * w[sel], np.clip(x0, lo, hi), bounds=(lo, hi),
                             x_scale="jac", max_nfev=2000).x

    keep = np.ones(t.size, bool)
    x = fit(keep)
    if enhancement:
        # a ground-level enhancement is an excess the reference does not explain: hours whose counts stand more than
        # 6 sigma above the fit, with the 2 hours before and 18 after each, are left out and the fit repeated (at
        # most a fifth of the record: a fit that misses everywhere is not an enhancement)
        for _ in range(4):
            z = (y - model(x, np.ones(t.size, bool))) * w
            hot = np.flatnonzero(z > 6.0)
            if hot.size == 0:
                break
            new = keep.copy()
            for i in hot:
                new[max(i - 2, 0):i + 19] = False
            if np.array_equal(new, keep) or new.sum() < 0.8 * t.size:     # an excess, not a misfit of the record
                break
            keep = new
            x = fit(keep, x)
    return x, u, y + offset[u] - model(x, np.ones(t.size, bool)), w, keep


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_station(sys.argv[1])))
