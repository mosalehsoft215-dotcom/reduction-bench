"""Reduces one compact-tension fatigue crack to its growth rate at dK_eff = 15 MPa m^0.5 and its Paris exponent.

    python3 crack_reduce.py <crack directory> <result.json>

1. Loads: each reading's peak and valley loads from a straight line through its 60 loads (the load setting can change
   during a test).
2. Opening load and compliance: each unloading is a continuous two-line fit, straight above the opening load and
   stiffer below it; the break, scanned over the lower 80 % of the load range and refined by a parabola, is P_op and
   the upper slope is the open-crack compliance.
3. Effective modulus: the pre-crack front (nine-point average) is the crack at count 0 and the final front the crack
   at stop_cycles. Where recording began after count 0 or cycling ran on after the last reading, each front is carried
   to that reading by the Paris law of the lengths it implies (fitted as in steps 5-6, at the load range and closure
   trend of the nearest readings at that load, fourth-order Runge-Kutta), solved for the anchor length by Brent's
   method. The modulus that puts the first reading's compliance length on its anchor is averaged with the one for the
   last reading; with an actuator-measured opening the two anchors fix both the modulus and the load-train
   compliance. B_e = B - (B - B_N)^2 / B, B_N being the span of the nine measurement positions across the fracture
   surface (the net section as machined).
4. Crack length of every reading by inverting E B_e v/P = F(a/W), or a/W(U) of the front-face relation where the
   clip gauge sits at the front face.
5. Rates by the seven-point incremental polynomial within each stretch of constant load, leaving out the readings
   inside the plastic zone r_p = (K_max / sigma_y)^2 / (2 pi) left by a load reduction, pre-cracking above the
   test's first P_max counting as one at count 0;
   dK_eff = 31.623 (P_max - P_op) f(a/W) / (sqrt(B B_N) sqrt(W)) at the fitted length.
6. Paris law fitted in log-log over 9 <= dK_eff <= 30 MPa m^0.5; the rate at dK_eff = 15 in nm per cycle and the
   exponent.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

KSCALE = 31.6227766


def f_ct(al):
    return (2 + al) / (1 - al) ** 1.5 * (0.886 + 4.64 * al - 13.32 * al ** 2 + 14.72 * al ** 3 - 5.6 * al ** 4)


def F_ct(al):
    return ((1 + al) / (1 - al)) ** 2 * (2.1630 + 12.219 * al - 20.065 * al ** 2 - 0.9925 * al ** 3
                                         + 20.609 * al ** 4 - 9.9314 * al ** 5)


def readings(path):
    data = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            k = int(r["reading"])
            d = data.setdefault(k, {"n": int(r["cycles"]), "p": [], "v": []})
            d["p"].append(float(r["load_kn"]))
            d["v"].append(float(r["cod_mm"]))
    return [data[k] for k in sorted(data)]


def load_range(rd):
    """Peak and valley load of the cycles at this reading, from its own unloading."""
    p = np.array(rd["p"])
    line = np.polyval(np.polyfit(np.arange(p.size), p, 1), [0, p.size - 1])
    return line[0], line[1]


def two_branch(rd):
    """Open-crack compliance and opening load of one unloading: a continuous two-line fit, straight above the
    opening load and stiffer below it; the break is scanned over the load range and refined by a parabola."""
    p, v = np.array(rd["p"]), np.array(rd["v"])
    top, bottom = load_range(rd)

    def fit(pc):
        a = np.column_stack([np.ones_like(p), np.maximum(p, pc), np.minimum(p, pc) - pc])
        coef, res, *_ = np.linalg.lstsq(a, v, rcond=None)
        return float(np.sum((a @ coef - v) ** 2)), coef[1]
    grid = bottom + (top - bottom) * np.linspace(0.0, 0.8, 161)
    sse = np.array([fit(x)[0] for x in grid])
    k = int(np.argmin(sse))
    pc = grid[k]
    if 0 < k < grid.size - 1:
        y0, y1, y2 = sse[k - 1], sse[k], sse[k + 1]
        den = y0 - 2 * y1 + y2
        if den > 0:
            pc = grid[k] + 0.5 * (y0 - y2) / den * (grid[1] - grid[0])
    return fit(pc)[1], max(pc, bottom)


def compliance(rd):
    return two_branch(rd)[0]


FRONT_FACE = (1.0010, -4.6695, 18.460, -236.82, 1214.9, -2143.6)   # a/W from U at X = -0.25 W (ASTM E647)


def front_face_alpha(g):
    u = 1 / (np.sqrt(g) + 1)
    return sum(c * u ** i for i, c in enumerate(FRONT_FACE))


def g_of(al, front):
    """E B_e v / P for a / W = al at the load line or at the front face."""
    if not front:
        return F_ct(al)
    return brentq(lambda g: front_face_alpha(g) - al, 1.0, 1e5)


def length(c, e, be, w, front=False):
    target = c * e * be
    if front:
        return w * front_face_alpha(target)
    return w * brentq(lambda al: F_ct(al) - target, 0.05, 0.9)


def nine_point(values):
    v = np.asarray(values, float)
    return ((v[0] + v[-1]) / 2 + v[1:-1].sum()) / 8


def fronts(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return (nine_point([r["precrack_front_mm"] for r in rows]), nine_point([r["final_front_mm"] for r in rows]))


def incremental_rates(n, a):
    out_n, out_a, out_r = [], [], []
    for i in range(3, len(n) - 3):
        nn, aa = n[i - 3:i + 4], a[i - 3:i + 4]
        c1, c2 = nn.mean(), (nn[-1] - nn[0]) / 2
        x = (nn - c1) / c2
        b2, b1, b0 = np.polyfit(x, aa, 2)
        xi = (n[i] - c1) / c2
        out_n.append(n[i])
        out_a.append(b0 + b1 * xi + b2 * xi ** 2)
        out_r.append(b1 / c2 + 2 * b2 * xi / c2)
    return np.array(out_a), np.array(out_r)


DK_LO, DK_HI = 9.0, 30.0           # Paris window in dK_eff, MPa m^0.5


def reduce_crack(folder, anchors="fronts", front_relation=True, effective=True, transients=True, index_lengths=False,
                 net="surface", counts=True, precrack=True, anchors_only=False):
    folder = Path(folder)
    sh = json.loads((folder / "sheet.json").read_text())
    w, b, bn = sh["width_mm"], sh["thickness_mm"], sh["net_thickness_mm"]
    if net == "surface":                                   # the net section as machined: the span of the nine-point
        with open(folder / "fracture_fronts.csv", newline="") as fh:   # positions across the fracture surface
            across = [float(r["position_mm"]) for r in csv.DictReader(fh)]
        bn = max(across) - min(across)
    be = b - (b - bn) ** 2 / b
    rds = readings(folder / "unloadings.csv")
    fits = [two_branch(r) for r in rds]
    comp = np.array([f[0] for f in fits])
    p_op = np.array([f[1] for f in fits])
    lr = np.array([load_range(r) for r in rds])
    dp_read = lr[:, 0] - lr[:, 1]
    eff_share = (lr[:, 0] - p_op) / dp_read if effective else np.ones(len(rds))
    n = np.array([r["n"] for r in rds], float)
    front = front_relation and "front face" in sh.get("opening_measured_by", "")
    a_lo, a_hi = (fronts(folder / "fracture_fronts.csv") if anchors == "fronts"
                  else (sh["optical_crack_start_mm"], sh["optical_crack_end_mm"]))
    f_start, f_end = g_of(a_lo / w, front), g_of(a_hi / w, front)
    # the pre-crack front is the crack at count 0 and the final front the crack at stop_cycles; recording can start
    # after count 0 and cycling can run on after the last stored reading, so the compliance is taken at those counts
    stop = float(sh.get("stop_cycles", n[-1]))
    late_start = counts and anchors == "fronts" and n[0] > 0
    late_stop = counts and anchors == "fronts" and stop > n[-1]
    actuator = sh.get("opening_measured_by", "").startswith("actuator")
    kfac = KSCALE / (np.sqrt(b * bn) * np.sqrt(w))
    zone = lambda k: (k / sh["yield_strength_mpa"]) ** 2 / (2 * np.pi) * 1e3      # plane-stress plastic zone, mm
    # the seven-point windows stay inside stretches of constant load: a load change bends a(N)
    cut = np.flatnonzero(np.abs(np.diff(dp_read)) > 0.02 * dp_read[:-1]) + 1
    bounds = np.r_[0, cut, len(n)]
    p_pre = sh.get("precrack_max_load_kn", 0.0)

    def lengths(l_first, l_last):
        """Crack length of every reading, the compliance calibrated so that the first reading sits at l_first and the
        last at l_last; with an actuator-measured opening the two also fix the load-train compliance C_m."""
        g_lo, g_hi = g_of(l_first / w, front), g_of(l_last / w, front)
        if actuator:
            e_eff = (g_hi - g_lo) / (be * (comp[-1] - comp[0]))
            c_m = comp[0] - g_lo / (e_eff * be)
        else:
            e_eff, c_m = 0.5 * (g_lo / (be * comp[0]) + g_hi / (be * comp[-1])), 0.0
        return np.array([length(c - c_m, e_eff, be, w, front) for c in comp])

    def paris(a):
        """Paris law (m, log C with the rate in m per cycle) from the seven-point rates of steady growth."""
        keep = np.ones(len(n), bool)
        if precrack and transients and p_pre > lr[0, 0]:
            # pre-cracking at a higher load is a load reduction at count 0: its plastic zone lies ahead of the
            # pre-crack front
            keep &= a >= a_lo + zone(kfac * p_pre * f_ct(a_lo / w))
        if transients:
            for i in cut:
                if dp_read[i] < dp_read[i - 1]:               # a load reduction: the old plastic zone retards growth
                    r_old = zone(kfac * lr[i - 1, 0] * f_ct(a[i - 1] / w))
                    keep &= ~((np.arange(len(n)) >= i) & (a < a[i - 1] + r_old))
        af, rate, dps, shares = [], [], [], []
        for lo_i, hi_i in zip(bounds[:-1], bounds[1:]):
            seg = np.arange(lo_i, hi_i)
            seg = seg[keep[seg]]
            if len(seg) < 7:
                continue
            sa, sr = incremental_rates(n[seg], a[seg])
            af.append(sa); rate.append(sr); dps.append(dp_read[seg][3:-3]); shares.append(eff_share[seg][3:-3])
        af, rate, dps, shares = map(np.concatenate, (af, rate, dps, shares))
        dk = kfac * dps * shares * f_ct(af / w)
        sel = (dk >= DK_LO) & (dk <= DK_HI) & (rate > 0)
        return np.polyfit(np.log(dk[sel]), np.log(rate[sel] * 1e-3), 1)

    def carry(a_from, cycles, a, end):
        """Length after `cycles` more cycles from a_from by the Paris law of these lengths, at the load range and the
        closure trend of the readings nearest to it at the load of reading `end` (fourth-order Runge-Kutta)."""
        m, logc = paris(a)
        same = np.flatnonzero(np.searchsorted(cut, np.arange(len(n)), side="right")
                              == np.searchsorted(cut, end, side="right"))   # readings at that reading's load
        near = same[np.argsort(np.abs(a[same] - a_from))[:20]]
        trend = np.polyfit(a[near], eff_share[near], 1)
        dp = float(np.median(dp_read[near]))
        rate = lambda x: np.exp(logc) * (kfac * dp * np.polyval(trend, x) * f_ct(x / w)) ** m * 1e3
        x, h = a_from, cycles / 200.0
        rate_safe = rate
        rate = lambda y: rate_safe(min(y, 0.9 * w))
        for _ in range(200):
            k1 = rate(x)
            k2 = rate(x + 0.5 * h * k1)
            k3 = rate(x + 0.5 * h * k2)
            k4 = rate(x + h * k3)
            x += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return x

    # the pre-crack front is the crack at count 0 and the final front the crack at stop_cycles; recording can begin
    # after count 0 and cycling can run on after the last stored reading, so each front is carried to its reading by
    # the growth law, solved together with the calibration it implies
    l_first, l_last = a_lo, a_hi
    for _ in range(2 if late_start and late_stop else 1):
        if late_stop:
            l_last = brentq(lambda x: carry(x, stop - n[-1], lengths(l_first, x), len(n) - 1) - a_hi, a_hi - 6.0, a_hi)
        if late_start:
            l_first = brentq(lambda x: carry(a_lo, n[0], lengths(x, l_last), 0) - x, a_lo, a_lo + 6.0)
    if anchors_only:                                   # for the evidence routes: the two anchor lengths
        return l_first, l_last
    a = lengths(l_first, l_last)
    if index_lengths:                                  # a shortcut kept for the evidence: lengths spread evenly by index
        a = np.linspace(a_lo, a_hi, len(rds))
    m, logc = paris(a)
    return {"rate_at_dkeff15_nm_per_cycle": float(np.exp(logc) * 15.0 ** m * 1e9), "paris_exponent": float(m)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_crack(sys.argv[1])))
