"""Closing-latch pick-up current and mechanism release-to-touch time of a generator circuit breaker, from its
synchronising records.

python3 close_timing.py <unit folder> <answer file>

Per closing (fault-recorder file dfr/Cnn.csv, synchroniser and station values in closings.csv):
1. Time base: the bus voltage's frequency before the command, as the recorder's time stamps show it, against the
   synchroniser's bus frequency; when they differ by more than 0.3 % the time stamps are rescaled.
2. Coil energisation t_e and time constant tau: A (1 - exp(-(t - t_e)/tau)) fitted to the closing-coil current before
   the armature dip and after it, up to the b-contact cut, with a first-order clamp lag where the onset asks for one.
   Only the shape of the clamp signal is used.
3. Contact touch of each pole (pole_touches): change points of the bus zero-sequence and across-pole voltages for the
   order and sample of each touch, and the currents the second and third poles let in for the instant within it.
4. Across closings: the battery bus voltage V (DFR channel VDC) drives the whole closing circuit, whose resistance is
   L/tau; the coil's part is copper at the cubicle temperature and the rest is the same on every closing. The latch
   releases at a fixed coil current I_pu, at t_e - (L/R) ln(1 - I_pu R / V); each pole then touches after its own
   fixed release-to-touch time c_p. I_pu and c_A, c_B, c_C are fitted to all closings; the answer is I_pu and the
   mean of c_p.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, least_squares

ALPHA_CU = 0.00393
PHASES = "ABC"


def load(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "breaker.json").read_text())
    with (folder / "closings.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    return folder, sheet, rows


def read_dfr(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: np.asarray(data[name], float) for name in data.dtype.names}


def sine_frequency(t_s, y):
    def model(t, a, f, ph, c):
        return a * np.sin(2 * np.pi * f * t + ph) + c
    a0 = 0.5 * (y.max() - y.min())
    best = None
    for ph0 in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        try:
            p, _ = curve_fit(model, t_s, y, p0=[a0, 50.0, ph0, 0.0], maxfev=4000)
        except RuntimeError:
            continue
        r = np.sum((model(t_s, *p) - y) ** 2)
        if best is None or r < best[0]:
            best = (r, p)
    return abs(best[1][1])


def coil_fit(t, i, allow_lag=True):
    """Energisation time, I_inf and tau from the rising coil current (allow_lag=False: clamp taken as lag-free)."""
    on = np.flatnonzero(i > 0.15)[0]
    # the armature dip: first fall of the smoothed current after it has risen
    smooth = np.convolve(i, np.ones(5) / 5, mode="same")
    n = np.arange(len(smooth) - 1)
    dip = int(np.flatnonzero((np.diff(smooth) < -0.02) & (n > on + 20))[0])
    # after the dip the current resumes the same rise until the b-contact cuts it (a fast fall)
    resume = int(np.searchsorted(t, t[dip] + 0.014))
    cut = int(np.flatnonzero((np.diff(smooth) < -0.05) & (n > resume))[0])
    sel = np.r_[max(on - 4, 0):dip - 6, resume:cut - 4]

    def model(tt, i_inf, t_e, tau):
        return np.where(tt > t_e, i_inf * (1 - np.exp(-(tt - t_e) / tau)), 0.0)
    p, _ = curve_fit(model, t[sel], i[sel], p0=[i[cut - 4] * 1.1, t[on] - 0.0003, 0.02], maxfev=20000)

    # the clamp may follow the current through a first-order lag tau_f: the rise it shows is then
    # A (1 - (tau exp(-x/tau) - tau_f exp(-x/tau_f)) / (tau - tau_f)), x = t - t_e, which starts with zero slope
    def lagged(tt, i_inf, t_e, tau, tau_f):
        x = np.clip(tt - t_e, 0.0, None)
        return np.where(tt > t_e, i_inf * (1 - (tau * np.exp(-x / tau) - tau_f * np.exp(-x / tau_f)) / (tau - tau_f)),
                        0.0)
    q, _ = curve_fit(lagged, t[sel], i[sel], p0=[p[0], p[1] - 0.0005, p[2], 0.0005],
                     bounds=([0.5 * p[0], p[1] - 0.004, 0.5 * p[2], 1e-6], [2 * p[0], p[1] + 0.001, 2 * p[2], 0.005]),
                     maxfev=20000)
    # the lag is kept only when the data clearly ask for it; otherwise it trades against t_e on a clean clamp
    r_pure = float(np.sum((i[sel] - model(t[sel], *p)) ** 2))
    r_lag = float(np.sum((i[sel] - lagged(t[sel], *q)) ** 2))
    if allow_lag and q[3] > 2e-4 and (r_pure - r_lag) / (r_lag / (sel.size - 4)) > 30:
        return q[1], q[0], q[2]
    return p[1], p[0], p[2]


class Wave:
    """A sum of sinusoids, each (a, b, w): a sin(wt) + b cos(wt)."""

    def __init__(self, parts):
        self.parts = list(parts)

    def __sub__(self, other):
        return Wave(self.parts + [(-a, -b, w) for a, b, w in other.parts])

    def scaled(self, k):
        return Wave([(k * a, k * b, w) for a, b, w in self.parts])

    def __call__(self, t):
        return sum(a * np.sin(w * t) + b * np.cos(w * t) for a, b, w in self.parts)

    def integral(self, t):
        return sum((-a * np.cos(w * t) + b * np.sin(w * t)) / w for a, b, w in self.parts)


def side_emfs(t, sel, v):
    """The three EMFs of one side, from its voltages while that side is still open: a common angular frequency from the
    space vector, then each phase's own sine and cosine."""
    s = v["A"][sel] + np.exp(2j * np.pi / 3) * v["B"][sel] + np.exp(-2j * np.pi / 3) * v["C"][sel]
    w = abs(np.polyfit(t[sel], np.unwrap(np.angle(s)), 1)[0])
    basis = np.vstack([np.sin(w * t[sel]), np.cos(w * t[sel])]).T
    return {p: Wave([(*np.linalg.lstsq(basis, v[p][sel], rcond=None)[0], w)]) for p in PHASES}


def pole_touches(rec, t, l_loop):
    """Contact touch of each pole. Before the first touch both sides run on their own EMFs. With pole a closed the
    floating bus side shifts by d_a = e_gen,a - e_bus,a; with a and b closed it sits at (d_a + d_b)/2; with all three at
    the mean of d. The voltage across an open pole p is then d_p minus that shift, and zero across a closed pole. The
    order and the three change samples are fitted to the zero-sequence and across voltages together; the second and
    third touches are then placed inside their sample interval from the currents they let in."""
    pre = t < 0.0
    dt = t[1] - t[0]
    if "VG_A" in rec:
        vg = {p: rec["VG_" + p] for p in PHASES}
    else:  # phase-to-phase VT on the generator side; the earthed generator carries no zero sequence
        ab, bc = rec["VG_AB"], rec["VG_BC"]
        vg = {"A": (2 * ab + bc) / 3, "B": (bc - ab) / 3, "C": -(ab + 2 * bc) / 3}
    vb = {p: rec["VB_" + p] for p in PHASES}
    zero = (vb["A"] + vb["B"] + vb["C"]) / 3
    # both sides are undisturbed before the command plus the shortest closing time; fit up to the first sign of contact
    sz = np.std(zero[pre])
    level = np.convolve(np.abs(zero), np.ones(4) / 4, mode="same")
    k_any = int(np.flatnonzero((level > 4 * sz) & (t > 0.02))[0])
    sel = (t > t[0] + 0.004) & (t < t[k_any] - 0.004)
    eg, eb = side_emfs(t, sel, vg), side_emfs(t, sel, vb)
    d = {p: eg[p] - eb[p] for p in PHASES}
    var_g = np.mean([np.var(vg[p][sel] - eg[p](t[sel])) for p in PHASES])
    var_b = np.mean([np.var(vb[p][sel] - eb[p](t[sel])) for p in PHASES])
    win = np.flatnonzero((t > 0.02) & (t < t[-1] - 0.01))
    tw = t[win]
    D = {p: d[p](tw) for p in PHASES}
    base = (eb["A"](tw) + eb["B"](tw) + eb["C"](tw)) / 3
    Z = zero[win]
    X = {p: (vg[p] - vb[p])[win] for p in PHASES}

    def seg_cost(shift, closed):
        c = (Z - base - shift) ** 2 / (var_b / 3)
        for p in PHASES:
            model = 0.0 if p in closed else D[p] - shift
            c = c + (X[p] - model) ** 2 / (var_g + var_b)
        return np.concatenate([[0.0], np.cumsum(c)])

    best = None
    for order in [(a, b, c) for a in PHASES for b in PHASES for c in PHASES if len({a, b, c}) == 3]:
        a, b, _ = order
        mean_d = (D["A"] + D["B"] + D["C"]) / 3
        P = [seg_cost(0.0 * Z, ()), seg_cost(D[a], (a,)), seg_cost((D[a] + D[b]) / 2, (a, b)),
             seg_cost(mean_d, tuple(PHASES))]
        n = Z.size
        # minimise P0[k1] + P1[k2]-P1[k1] + P2[k3]-P2[k2] + P3[n]-P3[k3] over k1 <= k2 <= k3
        A = P[0] - P[1]
        B = P[1] - P[2]
        C = P[2] - P[3]
        ia = np.zeros(n + 1, int)
        for k in range(1, n + 1):
            ia[k] = k if A[k] < A[ia[k - 1]] else ia[k - 1]
        AB = A[ia] + B
        ib = np.zeros(n + 1, int)
        for k in range(1, n + 1):
            ib[k] = k if AB[k] < AB[ib[k - 1]] else ib[k - 1]
        tot = AB[ib] + C
        k3 = int(np.argmin(tot))
        k2 = int(ib[k3])
        k1 = int(ia[k2])
        cost = tot[k3] + P[3][n]
        if best is None or cost < best[0]:
            best = (cost, order, (win[0] + k1, win[0] + k2, win[0] + k3), P)
    _, (first, second, third), (k1, k2, k3), P = best
    t1 = t[k1] - 0.5 * dt

    def voltage_cost(before, after):
        # chi-square of the voltages, apart from a constant, when this touch moves to sample k
        return lambda k: before[k - win[0]] - after[k - win[0]]
    var_i = np.var(rec["IG_A"][pre])
    loop = (rec["IG_" + first] - rec["IG_" + second]) / 2
    t2 = start_from_drive(t, loop, (d[first] - d[second]).scaled(1 / (2 * l_loop)), k2, k3, var_i / 2,
                          voltage_cost(P[1], P[2]), t[k1 - 1])
    mean_w = Wave(d["A"].parts + d["B"].parts + d["C"].parts).scaled(1 / 3)
    t3 = start_from_drive(t, rec["IG_" + third], (d[third] - mean_w).scaled(1 / l_loop), k3, k3 + 12, var_i,
                          voltage_cost(P[2], P[3]), t2)
    return {first: t1, second: t2, third: t3}


def start_from_drive(t, i, drive, k, k_end, var_i, voltage_cost, t_prev):
    """A pole's touch instant from the current it lets in and the voltages together. From the touch on, the current is
    the integral of the driving voltage over the loop inductance; the voltages change at the first sample after the
    touch. Each candidate instant near the change point k, after the previous touch t_prev and before k_end, is scored
    by both chi-squares; the current fit runs to k_end. The voltages decide when the voltage across the closing pole is
    large, the current when the drive is."""
    dt = t[1] - t[0]
    tt, ii = t[k - 8:k_end], i[k - 8:k_end]
    best = None
    for t0 in np.arange(max(t[k - 6], t_prev), t[min(k + 4, k_end - 1)], 0.005 * dt):
        model = np.where(tt > t0, drive.integral(tt) - drive.integral(t0), 0.0)
        chi = float(np.sum((ii - model) ** 2)) / var_i + voltage_cost(int(np.searchsorted(t, t0, side="right")))
        if best is None or chi < best[0]:
            best = (chi, t0)
    return best[1]


def reduce_unit(folder):
    folder, sheet, rows = load(folder)
    r20 = sheet["closing_coil"]["resistance_ohm_at_20c"]
    l_coil = sheet["closing_coil"]["inductance_h"]
    net = sheet["network"]
    l_loop = net["generator_subtransient_inductance_h"] + net["transformer_leakage_inductance_h_lv"]
    obs = []
    for row in rows:
        rec = read_dfr(folder / "dfr" / (row["closing"] + ".csv"))
        t = rec["t_ms"] / 1000.0
        pre = (t < -0.005)
        f_app = sine_frequency(t[pre], rec["VB_A"][pre])
        scale = f_app / float(row["sync_bus_hz"])
        if abs(scale - 1) > 0.003:
            t = t * scale
        t_e, i_inf, tau = coil_fit(t, rec["ICOIL"])
        touches = pole_touches(rec, t, l_loop)
        r_coil = r20 * (1 + ALPHA_CU * (float(row["cubicle_temp_c"]) - 20.0))
        obs.append(dict(v=float(np.mean(rec["VDC"])), tau=tau, r_coil=r_coil, i_inf=i_inf,
                        times={p: touches[p] - t_e for p in PHASES}))
    # the coil's time constant gives the whole loop resistance; what is not the coil (copper at the cubicle temperature)
    # is the rest of the circuit, the same on every closing
    r_rest = float(np.mean([l_coil / o["tau"] - o["r_coil"] for o in obs]))
    for o in obs:
        o["r_loop"] = o["r_coil"] + r_rest

    def resid(x):
        i_pu, c = x[0], dict(zip(PHASES, x[1:]))
        out = []
        for o in obs:
            latch = -(l_coil / o["r_loop"]) * math.log(max(1 - i_pu * o["r_loop"] / o["v"], 1e-9))
            out += [latch + c[p] - o["times"][p] for p in PHASES]
        return np.array(out)
    i_max = min(o["v"] / o["r_loop"] for o in obs) * 0.999
    # a touch made while the voltage across the pole is near zero is poorly fixed; a robust loss keeps one such closing
    # from pulling the fit
    fit = least_squares(resid, [0.6 * i_max, 0.035, 0.035, 0.035], bounds=([0.05, 0, 0, 0], [i_max, 0.2, 0.2, 0.2]),
                        loss="soft_l1", f_scale=1e-4)
    i_pu, c = fit.x[0], dict(zip(PHASES, fit.x[1:]))
    clamp = float(np.mean([o["i_inf"] * o["r_loop"] / o["v"] for o in obs]))
    answer = {"pickup_current_a": round(float(i_pu), 4),
              "release_to_touch_ms": round(1000 * float(np.mean([c[p] for p in PHASES])), 4)}
    notes = dict(r_rest=r_rest, clamp_gain=clamp, rms_ms=1000 * float(np.sqrt(np.mean(fit.fun ** 2))),
                 release_ms={p: 1000 * c[p] for p in PHASES})
    return answer, notes


def main():
    answer, notes = reduce_unit(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(answer, indent=1) + "\n", newline="\n")
    print(json.dumps(notes))


if __name__ == "__main__":
    main()
