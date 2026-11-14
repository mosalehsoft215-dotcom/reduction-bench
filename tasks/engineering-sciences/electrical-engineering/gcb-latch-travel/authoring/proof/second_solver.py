"""An independent reduction of one breaker's synchronising closings (shares no code with the reference).

python3 second_solver.py <unit folder> [answer file]

Time base from zero-crossing counts; coil current from log(di/dt); each side's EMFs from their own zero crossings; the
first touch from the collapse of the voltage across a pole; the second and third touches by root-finding on the
integral of the driving voltage; the pick-up current by a one-dimensional search with the pole offsets profiled out.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

A_CU = 3.93e-3


def columns(path):
    with open(path) as fh:
        rd = csv.reader(fh)
        head = next(rd)
        arr = np.array([[float(x) for x in r] for r in rd])
    return {h: arr[:, k] for k, h in enumerate(head)}


def crossings(t, y):
    """Times of upward zero crossings, by linear interpolation."""
    k = np.flatnonzero((y[:-1] < 0) & (y[1:] >= 0))
    return t[k] - y[k] * (t[k + 1] - t[k]) / (y[k + 1] - y[k])


def freq(t, y):
    z = crossings(t, y - np.median(y))
    return (z.size - 1) / (z[-1] - z[0]), z


def sinusoid(t, y, f):
    m = np.column_stack([np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)])
    return np.linalg.lstsq(m, y, rcond=None)[0]


class Emf:
    def __init__(self, terms):
        self.terms = terms  # (sin coef, cos coef, f)

    def at(self, t):
        return sum(a * np.sin(2 * np.pi * f * t) + b * np.cos(2 * np.pi * f * t) for a, b, f in self.terms)

    def area(self, t):
        return sum((b * np.sin(2 * np.pi * f * t) - a * np.cos(2 * np.pi * f * t)) / (2 * np.pi * f)
                   for a, b, f in self.terms)

    def minus(self, other, k=1.0):
        return Emf(self.terms + [(-k * a, -k * b, f) for a, b, f in other.terms])

    def times(self, k):
        return Emf([(k * a, k * b, f) for a, b, f in self.terms])


def coil_pure(t, i):
    """Energisation instant, final current and time constant: the rise is A (1 - exp(-(t - t0)/tau)), so its slope is
    exponential; tau from a line through log of the smoothed slope (armature dip and the cut left out)."""
    k0 = int(np.argmax(i > 0.2))
    d = np.gradient(np.convolve(i, np.ones(9) / 9, mode="same"), t)
    stop = k0 + 40 + int(np.argmin(d[k0 + 40:]))  # the b-contact cut: the steepest fall of all
    stop = k0 + 40 + int(np.argmax(d[k0 + 40:stop + 1] < -1500))
    kmin = k0 + 30 + int(np.argmin(d[k0 + 30:stop - 40]))  # deepest point of the armature dip
    fall = k0 + 30 + int(np.flatnonzero(d[k0 + 30:kmin] > 0)[-1])
    k_back = fall + int(round(0.015 / (t[1] - t[0])))
    keep = np.arange(k0 + 8, fall - 10)
    keep = keep[d[keep] > 0]
    slope, _ = np.polyfit(t[keep], np.log(d[keep]), 1)
    tau = -1.0 / slope
    sel = np.r_[k0 - 2:fall - 6, k_back:stop - 6]
    # with tau fixed the model is linear in A and A exp(t0/tau)
    m = np.column_stack([np.ones(sel.size), -np.exp(-t[sel] / tau)])
    for _ in range(3):
        a, b = np.linalg.lstsq(m, i[sel], rcond=None)[0]
        t0 = tau * np.log(b / a)
        sel = sel[t[sel] > t0]
        m = np.column_stack([np.ones(sel.size), -np.exp(-t[sel] / tau)])
    # refine tau and t0 together by a small grid around the slope estimate
    best = None
    for tt in tau * np.linspace(0.92, 1.08, 321):
        mm = np.column_stack([np.ones(sel.size), -np.exp(-t[sel] / tt)])
        coef, res, *_ = np.linalg.lstsq(mm, i[sel], rcond=None)
        r = float(np.sum((mm @ coef - i[sel]) ** 2))
        if best is None or r < best[0]:
            best = (r, tt, coef)
    _, tau, (a, b) = best
    t0 = tau * np.log(b / a)
    near = (t > t0 - 0.001) & (t < t0 + 0.006)
    fit = np.where(t[near] > t0, a - b * np.exp(-t[near] / tau), 0.0)
    return t0, a, tau, float(np.sum((i[near] - fit) ** 2))


LAGS = np.arange(0.0, 0.00305, 0.0001)


def coil_by_lag(t, i):
    """coil_pure after undoing each trial first-order lag of the clamp (x = i + lag di/dt), lag 0-3 ms."""
    slope = np.gradient(np.convolve(i, np.ones(5) / 5, mode="same"), t)
    out = []
    for lag in LAGS:
        with np.errstate(all="ignore"):
            t0, a, tau, r = coil_pure(t, i + lag * slope)
        out.append((t0, a, tau, r if np.isfinite(t0) and np.isfinite(r) else np.inf))
    return out


def touches(rec, t, l_loop):
    ph = "ABC"
    if "VG_A" in rec:
        g = {p: rec["VG_" + p] for p in ph}
    else:
        ab, bc = rec["VG_AB"], rec["VG_BC"]
        ca = -(ab + bc)
        g = {"A": (ab - ca) / 3, "B": (bc - ab) / 3, "C": (ca - bc) / 3}
    b = {p: rec["VB_" + p] for p in ph}
    x = {p: g[p] - b[p] for p in ph}
    # open-circuit EMFs from the pre-trigger record and the first 40 ms after the command (no pole can have touched)
    early = (t > t[0] + 0.002) & (t < 0.040)
    emf_g, emf_b = {}, {}
    for side, v, store in (("g", g, emf_g), ("b", b, emf_b)):
        f = np.mean([freq(t[early], v[p][early])[0] for p in ph])
        for p in ph:
            store[p] = Emf([(*sinusoid(t[early], v[p][early], f), f)])
    d = {p: emf_g[p].minus(emf_b[p]) for p in ph}
    noise = np.std(x["A"][early] - d["A"].at(t[early]))
    n = t.size
    k_from = int(np.searchsorted(t, 0.040))

    def collapse(p):
        m = (np.abs(x[p]) < 4 * noise) & (np.arange(n) >= k_from)
        run = np.convolve(m.astype(int), np.ones(6, int), mode="full")[5:]
        return int(np.argmax(run == 6))
    c = {p: collapse(p) for p in ph}
    first, second, third = sorted(ph, key=lambda p: c[p])
    # first touch: the bus zero sequence steps from 0 to d_first; choose the step sample near the collapse
    zero = (b["A"] + b["B"] + b["C"]) / 3
    lo, hi = c[first] - 6, c[second]
    dz = d[first].at(t[lo:hi])

    def sse(k):
        j = k - lo
        return float(np.sum(zero[lo:k] ** 2) + np.sum((zero[k:hi] - dz[j:]) ** 2))
    k1 = min(range(c[first] - 2, min(c[first] + 2, hi - 1) + 1), key=sse)
    dt = t[1] - t[0]
    t1 = t[k1] - 0.5 * dt
    i_loop = 0.5 * (rec["IG_" + first] - rec["IG_" + second])
    drive2 = d[first].minus(d[second]).times(1 / (2 * l_loop))
    drive3 = d[third].minus(Emf(d["A"].terms + d["B"].terms + d["C"].terms), 1 / 3).times(1 / l_loop)

    def start(i, drive, k, k_next):
        """Touch from the current: least absolute deviations of the integral model over the samples around the
        collapse, on a grid of 1/50 sample and then a bounded scalar search."""
        idx = np.arange(k - 8, min(k + 12, k_next))
        tt, ii = t[idx], i[idx]

        def cost(t0):
            return float(np.sum(np.abs(ii - np.where(tt > t0, drive.area(tt) - drive.area(t0), 0.0))))
        grid = np.arange(t[k] - 6 * dt, min(t[k] + 6 * dt, t[min(k + 12, k_next) - 1]), dt / 50)
        g = grid[int(np.argmin([cost(q) for q in grid]))]
        return float(minimize_scalar(cost, bounds=(g - dt / 50, g + dt / 50), method="bounded").x)
    t2 = start(i_loop, drive2, c[second], c[third])
    t3 = start(rec["IG_" + third], drive3, c[third], c[third] + 6)
    return {first: t1, second: t2, third: t3}


def solve(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "breaker.json").read_text())
    with open(folder / "closings.csv") as fh:
        rows = list(csv.DictReader(fh))
    r20, lc = sheet["closing_coil"]["resistance_ohm_at_20c"], sheet["closing_coil"]["inductance_h"]
    net = sheet["network"]
    l_loop = net["generator_subtransient_inductance_h"] + net["transformer_leakage_inductance_h_lv"]
    data = []
    for row in rows:
        rec = columns(folder / "dfr" / (row["closing"] + ".csv"))
        t = rec["t_ms"] * 1e-3
        f_rec, _ = freq(t[t < 0], rec["VB_A"][t < 0])
        k = f_rec / float(row["sync_bus_hz"])
        if abs(k - 1) > 0.003:
            t = t * k
        tt = touches(rec, t, l_loop)
        data.append(dict(by_lag=coil_by_lag(t, rec["ICOIL"]), touch=tt, v=float(np.median(rec["VDC"])),
                         rc=r20 * (1 + A_CU * (float(row["cubicle_temp_c"]) - 20))))
    # one clamp per unit: the lag that makes the onsets of all closings together most like a clean start, kept only
    # when it lowers their misfit by more than a fifth
    total = np.array([sum(q["by_lag"][k][3] for q in data) for k in range(LAGS.size)])
    k = int(np.argmin(total))
    if not total[k] < 0.8 * total[0]:
        k = 0
    for q in data:
        t0, _, q["tau"], _ = q["by_lag"][k]
        q["T"] = np.array([q["touch"][p] - t0 for p in "ABC"])
    rest = float(np.median([lc / q["tau"] - q["rc"] for q in data]))
    T = np.array([q["T"] for q in data])

    def latch(ip):
        out = []
        for q in data:
            r = q["rc"] + rest
            out.append(-(lc / r) * np.log(1 - ip * r / q["v"]))
        return np.array(out)

    keep = np.ones_like(T, bool)

    def cost(ip):
        res = np.where(keep, T - latch(ip)[:, None], np.nan)
        return float(np.nansum((res - np.nanmean(res, axis=0)) ** 2))
    hi = min(q["v"] / (q["rc"] + rest) for q in data) * 0.999
    # least squares, then once more without the touch times that sit more than 4 MAD (at least 0.15 ms) off
    for _ in range(2):
        ip = minimize_scalar(cost, bounds=(0.3 * hi, hi), method="bounded", options={"xatol": 1e-6}).x
        res = T - latch(ip)[:, None]
        dev = res - np.median(res, axis=0)
        mad = 1.4826 * np.median(np.abs(dev))
        keep = np.abs(dev) <= max(4 * mad, 1.5e-4)
    c = np.nanmean(np.where(keep, T - latch(ip)[:, None], np.nan), axis=0)
    return {"pickup_current_a": float(ip), "release_to_touch_ms": 1000 * float(np.mean(c))}


if __name__ == "__main__":
    ans = solve(sys.argv[1])
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(json.dumps(ans) + "\n")
    print(json.dumps(ans))
