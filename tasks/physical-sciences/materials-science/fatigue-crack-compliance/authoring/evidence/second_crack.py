"""Second, separately written reduction of a fatigue crack (same call as the reference).

Differences from solution/crack_reduce.py: each unloading is fitted by nonlinear least squares as a continuous
two-line curve with the opening load as a free parameter (no scan); B_N as the span of the nine positions across the
fracture surface (not the sheet's drawing value); the effective modulus from a joint fit at the first and the last
reading, whose lengths are the travelling-microscope surface lengths of those readings each raised by its broken
front's lead of the through-thickness average (trapezoidal) over the two surface points, rather than the fronts
carried to those readings by the growth law; the front-face relation used directly as a/W(U) where the gauge sits
there; and rates by the secant method from each reading to the first later reading of the same load level at least
0.3 mm further on, with dK_eff at the mid length from the mean opening load of the readings spanned, instead of the
seven-point polynomial. Readings inside the plastic zone left by a load reduction, or by pre-cracking above the test
load, are left out, and the Paris law is fitted over 9 <= dK_eff <= 30.
"""
import json
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq, least_squares, minimize_scalar


def k_shape(x):
    return (2 + x) / (1 - x) ** 1.5 * (0.886 + 4.64 * x - 13.32 * x ** 2 + 14.72 * x ** 3 - 5.6 * x ** 4)


def c_shape(x):
    return ((1 + x) / (1 - x)) ** 2 * (2.1630 + 12.219 * x - 20.065 * x ** 2 - 0.9925 * x ** 3
                                       + 20.609 * x ** 4 - 9.9314 * x ** 5)


# ASTM E647 compliance at the front face (X = -0.25 W): a/W as a polynomial in U = 1 / (sqrt(E B_e v / P) + 1)
FACE_COEF = np.array([1.0010, -4.6695, 18.460, -236.82, 1214.9, -2143.6])


def face_length(g, w):
    u = 1 / (np.sqrt(g) + 1)
    return w * np.polyval(FACE_COEF[::-1], u)


def average_front(column):
    # nine equally spaced points across the thickness, end points at half weight
    return float(np.trapezoid(column, dx=1.0) / (column.size - 1))


def reduce_crack(folder):
    folder = Path(folder)
    meta = json.loads((folder / "sheet.json").read_text())
    tab = np.genfromtxt(folder / "unloadings.csv", delimiter=",", names=True)
    w, thick = meta["width_mm"], meta["thickness_mm"]
    pos = np.genfromtxt(folder / "fracture_fronts.csv", delimiter=",", names=True)["position_mm"]
    net = float(np.ptp(pos))                    # net section measured on the fracture surface, not the drawing
    b_eff = thick - (thick - net) ** 2 / thick
    ids = np.unique(tab["reading"])
    slopes, counts, ranges, tops, opens = [], [], [], [], []
    for i in ids:
        sel = tab["reading"] == i
        p, v = tab["load_kn"][sel], tab["cod_mm"][sel]
        ends = np.polyval(np.polyfit(np.arange(p.size), p, 1), [0, p.size - 1])
        ranges.append(ends[0] - ends[1])
        tops.append(ends[0])
        up = p >= ends[1] + 0.8 * (ends[0] - ends[1])
        s0 = np.polyfit(p[up], v[up], 1)

        def kink(q):
            v0, s1, s2, pc = q
            return v0 + s1 * np.maximum(p, pc) + s2 * (np.minimum(p, pc) - pc) - v
        sol = least_squares(kink, [s0[1], s0[0], 0.7 * s0[0], ends[1] + 0.4 * (ends[0] - ends[1])],
                            bounds=([-np.inf, 0, 0, ends[1]], [np.inf, np.inf, np.inf, ends[0] - 0.15 * (ends[0] - ends[1])]),
                            x_scale="jac", xtol=1e-12, ftol=1e-12)
        slopes.append(sol.x[1])
        opens.append(sol.x[3])
        counts.append(tab["cycles"][sel][0])
    slopes, counts, ranges = np.array(slopes), np.array(counts, float), np.array(ranges)
    tops, opens = np.array(tops), np.array(opens)

    fr = np.genfromtxt(folder / "fracture_fronts.csv", delimiter=",", names=True)
    pre_avg = average_front(fr["precrack_front_mm"])
    # anchors at the first and the last reading: the travelling-microscope surface lengths taken at those readings,
    # raised by how far each broken front's through-thickness average leads its two surface points (the fronts
    # themselves belong to count 0 and to the stop of cycling, not to the first and the last reading)
    lead = [average_front(fr[k]) - 0.5 * (fr[k][0] + fr[k][-1]) for k in ("precrack_front_mm", "final_front_mm")]
    opt = np.array([meta["optical_crack_start_mm"] + lead[0], meta["optical_crack_end_mm"] + lead[1]])
    at_face = "front face" in meta.get("opening_measured_by", "")
    extra = 0.0
    if meta.get("opening_measured_by", "").startswith("actuator"):
        # slopes = F(a/W) / (E B_e) + C_m at both optical lengths: a 2 x 2 linear system in (1/E, C_m)
        lhs = np.array([[c_shape(opt[0] / w) / b_eff, 1.0], [c_shape(opt[1] / w) / b_eff, 1.0]])
        inv_e, extra = np.linalg.solve(lhs, slopes[[0, -1]])
        mod = 1 / inv_e

    def lengths(mod):
        if at_face:
            return np.array([face_length((s - extra) * mod * b_eff, w) for s in slopes])
        return np.array([w * brentq(lambda x: c_shape(x) - (s - extra) * mod * b_eff, 0.05, 0.9) for s in slopes])
    if extra == 0.0:
        mod = minimize_scalar(lambda e: np.sum((lengths(e)[[0, -1]] - opt) ** 2),
                              bounds=(20, 400), method="bounded", options={"xatol": 1e-6}).x
    a = lengths(mod)
    geo = lambda x, dp: 31.6227766 * dp * k_shape(x / w) / np.sqrt(thick * net * w)
    eff = (tops - opens) / ranges
    change = np.flatnonzero(np.abs(np.diff(ranges)) > 0.02 * ranges[:-1]) + 1
    usable = np.ones(a.size, bool)
    p_pre = meta.get("precrack_max_load_kn", 0.0)
    if p_pre > tops[0]:                         # pre-cracked above the test load: its plastic zone retards the start
        usable &= a >= pre_avg + (geo(pre_avg, p_pre) / meta["yield_strength_mpa"]) ** 2 / (2 * np.pi) * 1e3
    for c in change:
        if ranges[c] < ranges[c - 1]:
            r_old = (geo(a[c - 1], tops[c - 1]) / meta["yield_strength_mpa"]) ** 2 / (2 * np.pi) * 1e3
            usable[c:] &= a[c:] >= a[c - 1] + r_old
    seg_id = np.searchsorted(change, np.arange(a.size), side="right")
    rates, dks = [], []
    for i0 in np.flatnonzero(usable):
        # secant to the first later reading of the same stretch at least 0.3 mm further on (readings are logged at
        # a fixed cycle interval, so neighbours can be only a few hundredths of a millimetre apart)
        later = np.flatnonzero((np.arange(a.size) > i0) & usable & (seg_id == seg_id[i0]) & (a >= a[i0] + 0.3))
        if later.size == 0:
            continue
        i1 = later[0]
        if not usable[i0:i1 + 1].all():
            continue
        rates.append((a[i1] - a[i0]) / (counts[i1] - counts[i0]))           # secant rate, mm per cycle
        dks.append(geo(0.5 * (a[i0] + a[i1]), ranges[i1]) * np.mean(eff[i0:i1 + 1]))
    rates, dks = np.array(rates), np.array(dks)
    win = (dks >= 9.0) & (dks <= 30.0) & (rates > 0)
    m, logc = np.polyfit(np.log(dks[win]), np.log(rates[win] * 1e-3), 1)
    return {"rate_at_dkeff15_nm_per_cycle": float(np.exp(logc) * 15.0 ** m * 1e9), "paris_exponent": float(m)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_crack(sys.argv[1])))
