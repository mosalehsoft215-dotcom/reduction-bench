"""Reference reduction of one gauging station: volume passed and largest 15-minute mean discharge.

    python3 notch_flow.py <station folder> <answer.json>

1. Barometer onto the logger clock. The barometer keeps its own clock, which stands still while it has no power, so
   its lag behind the logger is a step function that falls by the length of every cut. The well pressure carries every
   air-pressure jump and wave on top of the slower head changes, so for a trial lag the compensated head is formed and
   its short-period roughness (absolute first difference on the 5-minute clock) is measured in 12-hour windows
   stepped by 4 hours, for every whole-minute lag over +-300. A step function with no, one or two falls is fitted
   to that table (each segment at the lag that is smoothest over its windows; a fall kept when it cuts the summed
   roughness by more than 3 %), and each fall is then placed on the barometer's own clock (5-minute steps, then
   minutes) where the head over the day around it is smoothest. Barometer samples go onto the logger clock with their
   segment's lag; the air pressure is interpolated across each cut.
2. Heads: h1 = (p_well - p_air) / (rho(T) g) - sensor depth, CIPM density at the logged temperature; h2 from the
   tailwater gauge.
3. Each dilution gauging: background = mean conductivity before the injection; excess concentration
   c = (EC - background) / k (mg/L, then kg/m3); the integral of c over the record. When the record stops before the
   trace is back at background, the whole observed passage is fitted (scipy curve_fit, three starts) with a Gaussian
   slug smeared by the stirred pool (exponentially modified Gaussian) and the model's area beyond the last sample is
   added; Q = M / integral. The concentration weights of the passage are kept.
4. Rating: Q = C (h1 + kh)^2.5 x Villemonte factor (1 - (h2/h1)^2.5)^0.385 while h2 > 0 (kh = 0.85 mm), plus the
   plate leak a sqrt(H), H = h1 minus the pool level, or plus the sill depth while the pool is below the sill. A
   gauging returns the concentration-weighted harmonic mean of Q(t) over its passage; C and a are fitted to the five
   gaugings by least squares on log discharge (scipy least_squares, four starts of a). A storm can change the gap
   under the plate: at each storm peak of the head record (local maxima 12 h apart, the four highest) the leak may
   step from a1 to a2 (C kept); the step is kept when its best fit cuts the misfit to under a quarter and moves the
   leak by more than 0.3 L/s per sqrt(m).
5. The stilling well (intake law dh_w/dt = k sign(h1 - h_w) sqrt|h1 - h_w|): storms are the head record's peaks
   that rise 3 cm or more over the week's lower quartile (12 h apart), each starting at the lowest level since the
   previous one. For each storm whose pool stands 4 mm or more above the crest within 6 h of its peak, k is the value
   (90-point grid 0.0015-0.3 then golden section) that makes the weir shape of the channel level, from the law with a
   25-minute Savitzky-Golay rate, a straight line against the pool level (the rising and the falling limb close into
   one line); k and the line are then refined together by least squares so that the well, stepped exactly through
   the law from the line's channel level, follows the logged well level. Each storm's k holds from its start on
   (a clean pipe, 0.07, where the pool tells nothing); the channel level follows from the law, and while the pool
   stands up, from the line.
6. The rating on every sample of the week, the 5-minute clock with gaps bridged linearly in discharge; volume = sum
   x 300 s; peak = largest mean of three consecutive values, in L/s.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit, least_squares
from scipy.special import erfc, erfcx

G = 9.81
KH = 0.00085                       # m, Kindsvater-Shen kh of a 90-degree notch


def density(t):
    return 999.97495 * (1.0 - (t - 3.983035) ** 2 * (t + 301.797) / (522528.9 * (t + 69.34881)))


def read_csv(path, cols):
    with open(path, newline="") as fh:
        return np.array([[float(r[c]) for c in cols] for r in csv.DictReader(fh)])


def load(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "station_sheet.json").read_text())
    well = read_csv(folder / "well_log.csv", ("minute", "p_kPa", "tailwater_mm", "water_C"))
    baro = read_csv(folder / "barometer_log.csv", ("minute", "p_air_hPa"))
    baro[:, 1] /= 10.0
    traces = [read_csv(folder / g["file"], ("minute", "ec_uS_cm")) for g in sheet["gaugings"]]
    return sheet, well, baro, traces


def head_series(well, baro, lag, depth):
    p_air = np.interp(well[:, 0] + lag, baro[:, 0], baro[:, 1])
    return (well[:, 1] - p_air) * 1000.0 / (density(well[:, 3]) * G) - depth


def roughness(h, minutes):
    clock = np.arange(minutes[0], minutes[-1] + 1, 5.0)
    return float(np.mean(np.abs(np.diff(np.interp(clock, minutes, h)))))


def barometer_lag(well, baro, depth):
    score = lambda lag: roughness(head_series(well, baro, lag, depth), well[:, 0])  # noqa: E731
    best = min(range(-180, 181), key=score)
    return min((best + k / 10.0 for k in range(-10, 11)), key=score)


def air_on_logger(well, baro, segments):
    """air pressure at the logger's minutes from barometer samples placed with their segment's lag:
    segments = [(clock where the segment starts, lag)], the clock standing still between them"""
    clock = baro[:, 0]
    true = np.empty_like(clock)
    for k, (c0, lag) in enumerate(segments):
        c1 = segments[k + 1][0] if k + 1 < len(segments) else np.inf
        part = (clock >= c0) & (clock < c1)
        true[part] = clock[part] - lag
    order = np.argsort(true)
    return np.interp(well[:, 0], true[order], baro[order, 1])


def stepped_lag(well, baro, depth):
    """barometer segments [(start clock, lag)] from a 12-h window table of head roughness against whole-minute lags"""
    minutes = well[:, 0]
    clock5 = np.arange(minutes[0], minutes[-1] + 1, 5.0)
    rho = density(np.interp(clock5, minutes, well[:, 3]))
    p_w = np.interp(clock5, minutes, well[:, 1])
    lags = np.arange(-300, 301)
    starts = np.arange(clock5[0], clock5[-1] - 720.0 + 1, 240.0)
    table = np.empty((lags.size, starts.size))
    for i, lag in enumerate(lags):
        h = (p_w - np.interp(clock5 + lag, baro[:, 0], baro[:, 1] / 1.0)) * 1000.0 / (rho * G)
        step = np.abs(np.diff(h))
        csum = np.r_[0.0, np.cumsum(step)]
        lo = np.searchsorted(clock5, starts)
        hi = np.searchsorted(clock5, starts + 720.0)
        table[i] = (csum[np.minimum(hi, step.size)] - csum[lo]) / np.maximum(hi - lo, 1)
    nw = starts.size

    def seg(a, b):                       # best single lag over windows a..b-1 and its summed roughness
        tot = table[:, a:b].sum(axis=1)
        k = int(np.argmin(tot))
        return float(tot[k]), int(lags[k])

    plans = [(seg(0, nw)[0], [(0, seg(0, nw)[1])])]
    for i in range(2, nw - 1):
        a, b = seg(0, i), seg(i, nw)
        plans.append((a[0] + b[0], [(0, a[1]), (i, b[1])]))
        for j in range(i + 2, nw - 1):
            c, d = seg(i, j), seg(j, nw)
            plans.append((a[0] + c[0] + d[0], [(0, a[1]), (i, c[1]), (j, d[1])]))
    best = {n: min((pl for pl in plans if len(pl[1]) == n), key=lambda q: q[0]) for n in (1, 2, 3)}
    choice = best[1]
    for n in (2, 3):
        steps = [best[n][1][k][1] - best[n][1][k + 1][1] for k in range(n - 1)]
        if best[n][0] < 0.97 * choice[0] and all(st >= 20 for st in steps):
            choice = best[n]
    plan = choice[1]
    segments = [(-np.inf, float(plan[0][1]))]
    for k in range(1, len(plan)):
        w0, lag_before, lag_after = plan[k][0], plan[k - 1][1], plan[k][1]
        t_mid = starts[w0] + 120.0                       # the fall lies about here on the logger clock
        near = (minutes > t_mid - 960.0) & (minutes < t_mid + 960.0)

        def rough(c):
            trial = segments + [(c, float(lag_after))]
            h = head_from(well, air_on_logger(well, baro, trial), depth)
            return roughness(h[near], minutes[near])

        cands = np.arange(t_mid - 720.0 + lag_before, t_mid + 720.0 + lag_before, 5.0)
        c5 = min(cands, key=rough)
        c1 = min(np.arange(c5 - 5.0, c5 + 5.1, 1.0), key=rough)
        segments.append((float(c1), float(lag_after)))
    return segments


def head_from(well, p_air, depth):
    return (well[:, 1] - p_air) * 1000.0 / (density(well[:, 3]) * G) - depth


def gauged_flow(entry, trace):
    """returns Q (m3/s), the sample minutes and the concentration weights of the observed passage"""
    t, ec = trace[:, 0], trace[:, 1]
    t0 = entry["injection_minute"]
    background = float(np.mean(ec[t < t0]))
    after = t >= t0
    ts, conc = t[after], (ec[after] - background) / entry["meter_uS_cm_per_mg_L"] / 1000.0   # kg/m3
    dt_s = float(np.median(np.diff(ts))) * 60.0
    area = float(np.sum(conc) * dt_s)
    peak = float(conc.max())
    k = int(np.argmax(conc))
    end = conc[-12:].mean()
    tail_area = 0.0
    if end > 0.005 * peak:
        # the record stopped early: fit the whole passage with a Gaussian slug smeared by the stirred pool
        # (exponentially modified Gaussian) and take the part of its area beyond the last sample
        x = (ts - t0) * 60.0
        def emg(x, area, mu, sigma, tau):
            z = (mu + sigma * sigma / tau - x) / (math.sqrt(2.0) * sigma)
            ahead = np.exp(-0.5 * ((x - mu) / sigma) ** 2) * erfcx(np.maximum(z, 0.0))
            behind = np.exp(np.minimum((mu - x) / tau + 0.5 * (sigma / tau) ** 2, 700.0)) * erfc(np.minimum(z, 0.0))
            return area / (2.0 * tau) * np.where(z >= 0.0, ahead, behind)
        mu0 = float(x[k]) - 60.0
        best = None
        for tau0 in (120.0, 220.0, 360.0):
            try:
                prm, _ = curve_fit(emg, x, conc, p0=(area * 1.1, mu0, 50.0, tau0),
                                   bounds=([0.0, 0.0, 5.0, 20.0], [np.inf, 3600.0, 600.0, 3000.0]), maxfev=4000)
            except RuntimeError:
                continue
            sse = float(np.sum((emg(x, *prm) - conc) ** 2))
            if best is None or sse < best[0]:
                best = (sse, prm)
        prm = best[1]
        grid = np.arange(x[-1] + 0.5 * dt_s, x[-1] + 12 * prm[3] + 10 * prm[2], 1.0)
        tail_area = float(np.sum(emg(grid, *prm)))
    weights = np.clip(conc, 0.0, None) * dt_s
    weights = np.append(weights, tail_area)
    minutes = np.append(ts, ts[-1] + 2.0)
    return entry["salt_kg"] / (area + tail_area), minutes, weights / weights.sum()


def gauging_cost(gaugings, well, weir, leak, c, a):
    tot = 0.0
    for q, minutes, w in gaugings:
        s = np.interp(minutes, well[:, 0], weir)
        l = np.interp(minutes, well[:, 0], leak)
        tot += (math.log(q) + math.log(float(np.sum(w / np.maximum(c * s + a * l, 1e-12))))) ** 2
    return 0.5 * tot


def drown(h1, h2):
    ratio = np.where((h2 > 0) & (h1 > 0), np.clip(h2 / np.where(h1 > 0, h1, 1.0), 0.0, 1.0), 0.0)
    return np.where((h2 > 0) & (h1 > 0), np.clip(1.0 - ratio ** 2.5, 0.0, None) ** 0.385, 1.0)


def flow_parts(h, h2, sill):
    """weir shape (h + kh)^2.5 x Villemonte factor, and the leak shape sqrt(drop): drop from the well level to the pool,
    or to the sill while the pool is below it"""
    weir = np.clip(h + KH, 0.0, None) ** 2.5 * drown(h, h2)
    leak = np.sqrt(np.clip(h - np.where(h2 > -sill, h2, -sill), 0.0, None))
    return weir, leak


def gauging_misfit(gaugings, well, h1, h2, sill, leak_on, split=None):
    """weir and leak shapes on the week, and C and a fitted to the five gaugings (each a concentration-weighted
    harmonic mean of the flow over the salt's passage) on log discharge; with split (a logger minute) the leak is a1
    before it and a2 after it"""
    weir, leak = flow_parts(h1, h2, sill)
    if split is not None:
        return split_fit(gaugings, well, weir, leak, split)
    rows = [(q, np.interp(minutes, well[:, 0], weir), np.interp(minutes, well[:, 0], leak), w) for q, minutes, w in gaugings]

    def res(z):
        return np.array([math.log(q) + math.log(float(np.sum(w / np.maximum(z[0] * s + z[1] * l, 1e-12))))
                         for q, s, l, w in rows])
    c0 = math.exp(np.mean([math.log(q * float(np.sum(w / np.maximum(s, 1e-12)))) for q, s, l, w in rows]))
    if not leak_on:
        return (c0, 0.0), weir, leak
    best = None
    for a0 in (0.0, 0.001, 0.002, 0.003):
        got = least_squares(res, [min(max(c0, 0.6), 2.9), a0], bounds=([0.5, 0.0], [3.0, 0.01]), x_scale=[0.1, 0.001])
        if best is None or got.cost < best.cost:
            best = got
    return tuple(best.x), weir, leak


def split_fit(gaugings, well, weir, leak, split):
    rows = [(q, np.interp(minutes, well[:, 0], weir), np.interp(minutes, well[:, 0], leak), w,
             float(np.sum(w * (minutes >= split)))) for q, minutes, w in gaugings]

    def res(z):
        out = []
        for q, s, l, w, after in rows:
            a = z[1] * (1.0 - after) + z[2] * after
            out.append(math.log(q) + math.log(float(np.sum(w / np.maximum(z[0] * s + a * l, 1e-12)))))
        return np.array(out)
    best = None
    for a0 in (0.0, 0.0015, 0.003):
        for a1 in (0.0, 0.0015, 0.003):
            got = least_squares(res, [1.35, a0, a1], bounds=([0.5, 0.0, 0.0], [3.0, 0.01, 0.01]),
                                x_scale=[0.1, 0.001, 0.001])
            if best is None or got.cost < best.cost:
                best = got
    return tuple(best.x), weir, leak, float(best.cost)


def storm_peaks(minutes, h1, count=4, apart=720.0):
    order = np.argsort(h1)[::-1]
    peaks = []
    for i in order:
        if 0 < i < h1.size - 1 and all(abs(minutes[i] - m) > apart for m in peaks):
            peaks.append(float(minutes[i]))
        if len(peaks) == count:
            break
    return peaks


CLEAN_K = 0.07                     # m^0.5/min, a clean intake: the well within a millimetre of the channel
K_GRID = np.geomspace(0.0015, 0.3, 90)


def channel_level(minutes, h_well, k):
    """invert the intake law dh_w/dt = k sign(h_c - h_w) sqrt|h_c - h_w|: h_c = h_w + sign(dh_w/dt) (dh_w/dt / k)^2,
    with the well's rate from a 25-minute quadratic Savitzky-Golay derivative (k may vary along the record)"""
    rate = savgol_filter(h_well, 5, 2, deriv=1, delta=5.0)
    return h_well + np.sign(rate) * (rate / k) ** 2


def level_from_shape(s, h2):
    """the channel level whose weir shape (h + kh)^2.5 x Villemonte factor is s (bisection, the shape rises with h)"""
    lo, hi = np.zeros_like(s), np.full_like(s, 0.8)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        up = np.clip(mid + KH, 0.0, None) ** 2.5 * drown(mid, h2) < s
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return 0.5 * (lo + hi)


def storms(minutes, h1, rise=0.03, apart=720.0):
    """storm peaks (rises of more than 3 cm over the week's lower quartile, 12 h apart) and the minute each storm
    starts (the lowest level since the previous peak)"""
    base = float(np.percentile(h1, 25))
    cands = sorted(m for m in storm_peaks(minutes, h1, count=6, apart=apart)
                   if float(np.interp(m, minutes, h1)) > base + rise)
    peaks, starts, prev = [], [], minutes[0]
    for m in cands:
        span = (minutes >= prev) & (minutes <= m)
        low = int(np.argmin(h1[span])) if span.any() else 0
        if span.any() and h1[span][low] < float(np.interp(m, minutes, h1)) - rise:   # a new rise, not a recession
            peaks.append(m)
            starts.append(float(minutes[span][low]))
            prev = m
    return peaks, starts


def intake_constant(minutes, h1, h2, peak):
    """k for one storm from the pool gauge: the pool level follows the discharge with no pipe in between, so with the
    right k the channel level of the rising and of the falling limb lie on one smooth curve against the pool level
    (the loop closes); None when the pool stays below the crest around this peak"""
    near = (np.abs(minutes - peak) < 6 * 60.0) & (h2 > 0.004)
    if near.sum() < 8:
        return None

    x = h2[near]
    design = np.column_stack([np.ones(x.size), x])

    def line(k):
        hc = channel_level(minutes, h1, k)[near]
        s = np.clip(hc + KH, 0.0, None) ** 2.5 * drown(hc, x)
        coef = np.linalg.lstsq(design, s, rcond=None)[0]
        return coef, float(np.sum((s - design @ coef) ** 2))

    def spread(k):
        return line(k)[1]
    vals = [spread(k) for k in K_GRID]
    if min(vals) > 0.6 * vals[-1]:                           # the loop does not close: the intake is clean
        return None
    j = int(np.argmin(vals))
    lo, hi = math.log(K_GRID[max(j - 1, 0)]), math.log(K_GRID[min(j + 1, len(K_GRID) - 1)])
    g = (math.sqrt(5) - 1) / 2
    for _ in range(30):
        a, b = hi - g * (hi - lo), lo + g * (hi - lo)
        if spread(math.exp(a)) < spread(math.exp(b)):
            hi = b
        else:
            lo = a
    k = math.exp(0.5 * (lo + hi))
    coef = line(k)[0]
    # forward refinement: inside the window drive the well with the channel level of the pool line and fit k and the
    # line to the logged well level (the noise sits in the observations, so the fit is not biased towards a large k)
    idx = np.where(near)[0]
    run = idx[np.argmax(np.diff(np.concatenate([[0], np.cumsum(np.diff(idx) != 1)])) == 0)] if idx.size else idx
    seg = np.arange(idx[0], idx[-1] + 1)
    seg = seg[h2[seg] > 0.004]
    tm, xw, hw = minutes[seg], h2[seg], h1[seg]

    def simulate(z):
        kk, a, b = math.exp(z[0]), z[1], z[2]
        hc = level_from_shape(a + b * xw, xw)
        out = np.empty_like(hw)
        out[0] = hw[0]
        for i in range(1, hw.size):
            dt = tm[i] - tm[i - 1]
            w = out[i - 1]
            for j in range(10):                              # exact orifice steps, channel level interpolated
                c = hc[i - 1] + (hc[i] - hc[i - 1]) * (j + 1) / 10.0
                d = c - w
                r = max(math.sqrt(abs(d)) - 0.5 * kk * dt / 10.0, 0.0)
                w = c - math.copysign(r * r, d)
            out[i] = w
        return out

    if seg.size >= 8:
        z0 = [math.log(k), coef[0], coef[1]]
        got = least_squares(lambda z: (simulate(z) - hw)[1:] / 0.0008, z0, x_scale=[0.3, abs(coef[0]) + 1e-4, abs(coef[1]) + 1e-3])
        k, coef = math.exp(got.x[0]), got.x[1:]
    return k, near, coef


def intake_series(minutes, h1, h2):
    """k along the record: each storm's own k (from the pool gauge) from the start of that storm on, a clean intake
    where the pool gives nothing"""
    peaks, starts = storms(minutes, h1)
    k = np.full(minutes.size, CLEAN_K)
    found, curves = [], []
    for peak, start in zip(peaks, starts):
        got = intake_constant(minutes, h1, h2, peak)
        kk = CLEAN_K if got is None else min(got[0], CLEAN_K)
        k[minutes >= start] = kk
        found.append((peak, kk))
        if got is not None and got[0] < CLEAN_K:
            curves.append(got[1:])
    return k, found, curves


def reduce(folder, leak_on=True, cuts=True, leak_steps=True, intake=True):
    """intake=True is the reference; intake=False reads the well as the channel, "first" takes the first storm's k
    for the week and "gaugings" the one k that fits the gaugings best: routes of authoring/evidence"""
    """leak_on=False leaves the plate leak out (C alone from the gaugings); cuts=False takes one barometer lag for the
    week; leak_steps=False one leak for the week: routes of authoring/evidence"""
    sheet, well, baro, traces = load(folder)
    depth = sheet["well_logger"]["sensor_below_crest_m"]
    sill = sheet["plate"]["sill_below_crest_m"]
    if cuts:
        segments = stepped_lag(well, baro, depth)
        h1 = head_from(well, air_on_logger(well, baro, segments), depth)
        lag = segments[0][1]
    else:
        lag = barometer_lag(well, baro, depth)
        h1 = head_series(well, baro, lag, depth)
    h2 = well[:, 2] / 1000.0
    found = []
    if intake is True:
        kser, found, curves = intake_series(well[:, 0], h1, h2)
        h1 = channel_level(well[:, 0], h1, kser)
        for near, coef in curves:                           # while the pool is up, its line gives the channel level
            h1[near] = level_from_shape(coef[0] + coef[1] * h2[near], h2[near])
    elif intake == "first":                                 # route: the first storm's k for the whole week
        found = intake_series(well[:, 0], h1, h2)[1]
        h1 = channel_level(well[:, 0], h1, next((k for _, k in found if k < CLEAN_K), CLEAN_K))
    elif intake == "gaugings":                              # route: one k for the week, the one the gaugings like best
        g0 = [gauged_flow(entry, trace) for entry, trace in zip(sheet["gaugings"], traces)]

        def fit_k(k):
            hc = channel_level(well[:, 0], h1, k)
            (cc, aa), wr, lk = gauging_misfit(g0, well, hc, h2, sill, leak_on)
            return gauging_cost(g0, well, wr, lk, cc, aa)
        kbest = min(np.geomspace(0.003, CLEAN_K, 16), key=fit_k)
        if fit_k(kbest) > 0.75 * fit_k(CLEAN_K):            # only a clear gain over a clean pipe moves k
            kbest = CLEAN_K
        h1 = channel_level(well[:, 0], h1, kbest)
        found = [(0.0, float(kbest))]
    gaugings = [gauged_flow(entry, trace) for entry, trace in zip(sheet["gaugings"], traces)]
    (c, a), weir, leak = gauging_misfit(gaugings, well, h1, h2, sill, leak_on)
    q = c * weir + a * leak
    if leak_on and leak_steps:
        base = gauging_cost(gaugings, well, weir, leak, c, a)
        tries = [gauging_misfit(gaugings, well, h1, h2, sill, True, split=m) + (m,) for m in storm_peaks(well[:, 0], h1)]
        (c2, a1, a2), _, _, cost, m = min(tries, key=lambda r: r[3])
        if cost < 0.25 * base and abs(a2 - a1) > 0.0003:
            c, a = c2, a1
            q = c2 * weir + np.where(well[:, 0] >= m, a2, a1) * leak
    clock = np.arange(well[0, 0], well[-1, 0] + 1, 5.0)
    flows = np.interp(clock, well[:, 0], q)
    return {"volume_m3": float(flows.sum() * 300.0),
            "peak_15min_L_s": float(np.max(np.convolve(flows, np.ones(3) / 3.0, mode="valid")) * 1000.0),
            "c": float(c), "leak": float(a), "lag_min": float(lag), "intake": [[round(p), round(k, 5)] for p, k in found]}


def main():
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)


if __name__ == "__main__":
    main()
