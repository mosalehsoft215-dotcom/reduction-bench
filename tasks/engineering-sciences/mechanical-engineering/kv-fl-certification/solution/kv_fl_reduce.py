"""Reference reduction of one valve specimen: Kv and FL of the valve itself at rated travel.

    python3 kv_fl_reduce.py <specimen folder> <answer.json>

1. Pump-off spells. With the valve shut (travel below 1 %, at the start, at the end and after any trip) the pump is
   off: both taps carry the tank's head, so the differential transmitter reads its zero, and the gauge transmitter reads
   the tank's head (held by the overflow) plus its own column. Air enters or leaves an impulse line when the pump
   stops, so the median dp of each spell is the dp zero in force until the next spell, and the change of the median p1
   from the first spell is a change of the upstream column, subtracted from p1 until the next spell. The meter cannot
   mark these spells: its damped figure is still running down there.
2. Smoothing. The shown flow, its time derivative, dp and p1 are taken through a 20 s Savitzky-Golay filter (second
   order), which also removes the pump pulsation.
3. The surge vessel keeps the loop moving through every hold, so no hold is averaged. Every rated-travel sample with
   the pump on (away from the edges of each stretch, every second sample) is one equation: the converter shows q with
   q + tau dq/dt = Q + q0, Q the installed relations Qu = Fp Kv sqrt(dp/100/G), Qc = FLP Kv sqrt((p1 + e - FF pv)/100/G)
   (IEC 60534-2-1, N2 = 0.0016, zeta1 = 0.5 (1 - (d/D)^2)^2, zeta2 = (1 - (d/D)^2)^2, zetaB1 = 1 - (d/D)^4) joined by
   the trim's transition Q = Qu (1 + (Qu/Qc)^n)^(-1/n); p1 is the gauge reading plus the ISO 2533 atmosphere at the lab
   elevation plus the unknown impulse-line column e; G, pv and FF come from the water temperature.
4. Kv, FL, n, e, q0 and tau are fitted together by least squares from eighteen starts; e is fixed by the choked flows
   of the three upstream levels, q0 by the unchoked flows, n by the bend between them and tau by the transients that
   every throttle move leaves in the record.
The functions holds, zero_track, corrected_holds and fit keep the hold-mean routes of authoring/evidence (the lab's
rule of settled holds and their last 45 s, and the same with each window mean lifted by tau).
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import savgol_filter


def density(t):
    return 999.97495 * (1.0 - (t - 3.983035) ** 2 * (t + 301.797) / (522528.9 * (t + 69.34881)))


def vapour(t):
    return 0.133322 * 10.0 ** (8.07131 - 1730.63 / (233.426 + t))


def load(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "specimen_sheet.json").read_text())
    with open(folder / "daq_record.csv", newline="") as fh:
        rows = np.array([[float(v) for v in r.values()] for r in csv.DictReader(fh)])
    return sheet, rows        # columns: t, travel, flow, p1 gauge, dp, temp


def smooth(x, n):
    pad = np.concatenate([np.full(n, x[0]), x, np.full(n, x[-1])])
    return np.convolve(pad, np.ones(2 * n + 1) / (2 * n + 1), mode="same")[n:-n]


def holds(rows, rate, whole=False, settled=False):
    """settled=True: the lab's rule, runs of at least 50 s in which dp and travel are steady (smoothed rates below
    0.15 kPa/s and 0.05 %/s), each represented by its last 45 s. Otherwise the record is cut at the throttle and travel
    moves (smoothed dp rate above 0.3 kPa/s, or travel moving), and every stretch between two moves that lasts 30 s or
    more with the pump on is a hold whether or not it has settled, with its last 45 s (whole=True: the whole stretch
    after its first 5 s)"""
    t, travel, flow, dp = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 4]
    w = int(5 * rate)
    ddp = np.abs(np.gradient(smooth(dp, w), t))
    dtr = np.abs(np.gradient(smooth(travel, w), t))
    if settled:
        steady, least = (ddp < 0.15) & (dtr < 0.05) & (flow > 0.0), 50.0
    else:
        steady, least = (ddp < 0.3) & (dtr < 0.05) & (flow > 0.0), 30.0
    out, start = [], None
    for i, ok in enumerate(np.append(steady, False)):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if t[i - 1] - t[start] >= least:
                if whole:
                    out.append(np.arange(start + int(5 * rate), i))
                else:
                    last = t[i - 1] - 45.0
                    out.append(np.where((t >= last) & (t <= t[i - 1]) & (np.arange(len(t)) >= start))[0])
            start = None
    return out


def zero_track(rows, rate):
    """dp zero in force at every sample for the hold routes: the median dp of the latest pump-off spell before it (the
    first spell's for samples ahead of it), taken as it was read with the pump off"""
    index, dzero, _, _ = spells(rows, rate)
    return dzero[index]


def meter_damping(rows, rate, spans, zero, patm):
    """the converter's first-order damping tau (s): inside every hold of a part travel and in the lower half of the
    rated-travel holds the valve is far from choking, so the true flow is K_h sqrt(dp / G) with one K_h per hold while
    the meter shows q with q + tau dq/dt = true flow + its zero; one linear least-squares fit over all those samples
    (a K_h per hold, a shared zero and tau) with dq/dt from a 10 s Savitzky-Golay derivative"""
    w = int(10 * rate) | 1
    dq = savgol_filter(rows[:, 2], w, 2, deriv=1, delta=1.0 / rate)
    cols, ys, blocks = [], [], []
    for idx in spans:
        m = rows[idx]
        dp = m[:, 4] - zero[idx]
        if m[:, 1].mean() > 99.0 and dp.mean() > 0.45 * (m[:, 3].mean() + patm):
            continue
        g = density(m[:, 5]) / 999.10
        blocks.append((np.sqrt(np.clip(dp, 0.0, None) / 100.0 / g), -dq[idx], m[:, 2]))
    n = len(blocks)
    for j, (s, d, y) in enumerate(blocks):
        part = np.zeros((len(s), n + 2))
        part[:, j], part[:, n], part[:, n + 1] = s, d, 1.0
        cols.append(part)
        ys.append(y)
    coef = np.linalg.lstsq(np.vstack(cols), np.concatenate(ys), rcond=None)[0]
    return max(float(coef[n]), 0.0)


def corrected_holds(sheet, rows, damping=True):
    rate = sheet["daq"]["rate_Hz"]
    zero = zero_track(rows, rate)
    patm = 101.325 * (1.0 - 2.25577e-5 * sheet["lab"]["elevation_m"]) ** 5.25588
    tau = meter_damping(rows, rate, holds(rows, rate, whole=True), zero, patm) if damping else 0.0
    smooth_q = savgol_filter(rows[:, 2], int(10 * rate) | 1, 1)
    pts = []
    for idx in holds(rows, rate, settled=not damping):
        m = rows[idx].mean(axis=0)
        # mean of the true flow over the window = mean of the shown flow + tau (shown at the end - at the start) / T
        lift = tau * (smooth_q[idx[-1]] - smooth_q[idx[0]]) / max(rows[idx[-1], 0] - rows[idx[0], 0], 1.0)
        dp_zero = float(zero[idx].mean())
        pv = vapour(m[5])
        pts.append({"travel": m[1], "q": m[2] + lift, "p1": m[3] + patm, "dp": m[4] - dp_zero,
                    "g": density(m[5]) / 999.10, "fpv": (0.96 - 0.28 * math.sqrt(pv / 22064.0)) * pv})
    return pts, tau


def model(x, d, big, a):
    kv, fl, n, e, q0 = x[0], x[1], math.exp(x[2]), x[3], x[4]
    r2 = (d / big) ** 2
    load_ = (kv / d ** 2) ** 2 / 0.0016
    fp = 1.0 / math.sqrt(1.0 + 1.5 * (1 - r2) ** 2 * load_)
    flp = fl / math.sqrt(1.0 + fl * fl * (0.5 * (1 - r2) ** 2 + 1.0 - r2 ** 2) * load_)
    qu = fp * kv * np.sqrt(a["dp"] / 100.0 / a["g"])
    qc = flp * kv * np.sqrt(np.maximum(a["p1"] + e - a["fpv"], 1e-6) / 100.0 / a["g"])
    return qu * (1.0 + (qu / qc) ** n) ** (-1.0 / n) + q0


def fit(pts, d, big, free=("n", "e", "q0")):
    pts = [p for p in pts if p["dp"] > 0.5]
    a = {k: np.array([p[k] for p in pts]) for k in ("q", "p1", "dp", "g", "fpv")}
    qmax = float(a["q"].max())
    k0 = float(np.nanmax(a["q"] * np.sqrt(a["g"] / (np.maximum(a["dp"], 0.05) / 100.0))))
    lo = [1e-3, 0.4, math.log(1.2) if "n" in free else math.log(399.0), -80.0 if "e" in free else -1e-9,
          -0.15 * qmax if "q0" in free else -1e-12]
    hi = [1e4, 0.995, math.log(400.0), 80.0 if "e" in free else 1e-9, 0.15 * qmax if "q0" in free else 1e-12]
    best = None
    for f0 in (0.7, 0.85, 0.95):
        for n0 in ((3.0, 10.0, 60.0) if "n" in free else (399.5,)):
            x0 = np.clip([k0 * 1.1, f0, math.log(n0), 0.0, 0.0], np.array(lo) + 1e-12, np.array(hi) - 1e-12)
            r = least_squares(lambda x: (model(x, d, big, a) - a["q"]) / qmax, x0, bounds=(lo, hi),
                              x_scale=[k0, 0.1, 1.0, 10.0, 0.01 * qmax])
            if best is None or r.cost < best.cost:
                best = r
    return float(best.x[0]), float(best.x[1]), best.x


def spells(rows, rate):
    """closed-valve spells (travel below 1 %, at least 20 s): the pump is off there, both taps carry the tank's head, dp
    reads its zero and p1 the head plus the transmitter's column, while a damped meter is still showing the flow running
    down; a closed-valve spell with the pump running (dp far above the first spell's) stops the flow but leaves the
    pressures up. Returns the pump-off spell in force at every sample (the latest one that began before it, the first
    for samples ahead of the first), the median dp and p1 of each pump-off spell, and the pump run of every sample
    (runs end at every closed-valve spell, where the flow through the meter stops)"""
    shut = rows[:, 1] < 1.0
    found, start = [], None
    for i, ok in enumerate(np.append(shut, False)):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            if i - start >= 20 * rate:
                found.append((start, float(np.median(rows[start + int(5 * rate):i, 4])),
                              float(np.median(rows[start + int(5 * rate):i, 3]))))
            start = None
    # pump off: both taps carry the tank's head, so dp is near its zero; a spell whose dp stays far above the first
    # spell's is a valve shut with the pump running (the flow stops, the pressures stay up)
    off = [f for f in found if abs(f[1] - found[0][1]) < 20.0]
    index, run = np.zeros(len(rows), dtype=int), np.zeros(len(rows), dtype=int)
    for j, (a, _, _) in enumerate(off):
        index[a:] = j
    for j, (a, _, _) in enumerate(found):
        run[a:] = j
    return index, np.array([f[1] for f in off]), np.array([f[2] for f in off]), run


def tap_state(x_e, reading, dpread, spell, dzero, pread, patm, boyle=True):
    """absolute upstream tap pressure and true dp from the readings: the p1 transmitter reads the tap pressure less the
    atmosphere less e (its column, fitted) plus the column change of the air taken in since the first spell; the dp
    transmitter reads dp plus its own zero (first spell) plus the change of the air in its high-pressure leg. Both air
    changes are read at the pump-off tap pressure p_off and shrink with the tap pressure p as p_off / p (Boyle), so p
    solves p^2 - (R + patm + e) p + dc p_off = 0"""
    p_off = patm + pread[0] + x_e
    dc = (pread - pread[0])[spell]
    dz = (dzero - dzero[0])[spell]
    s_ = reading + patm + x_e
    if boyle:
        p = 0.5 * (s_ + np.sqrt(np.maximum(s_ * s_ - 4.0 * dc * p_off, 1e-6)))
        return p, dpread - dzero[0] - dz * p_off / p
    return s_ - dc, dpread - dzero[0] - dz


def sample_fit(sheet, rows, tau_free=True, free=("n", "e", "q0"), boyle=True, stretch_zero=True):
    """every rated-travel sample with the pump on: the meter shows q with q + tau dq/dt = Q(dp, p1) + q0 (first-order
    converter damping, q0 the meter's zero, a new one after every stop of the flow), Q the installed-valve law at that
    instant with p1 and dp from tap_state; Kv, FL, n, e, tau and one q0 per pump run by least squares"""
    rate = sheet["daq"]["rate_Hz"]
    patm = 101.325 * (1.0 - 2.25577e-5 * sheet["lab"]["elevation_m"]) ** 5.25588
    spell, dzero, pread, run = spells(rows, rate)
    w = int(20 * rate) | 1
    q = savgol_filter(rows[:, 2], w, 2)
    dq = savgol_filter(rows[:, 2], w, 2, deriv=1, delta=1.0 / rate)
    dpr = savgol_filter(rows[:, 4], w, 2)
    p1r = savgol_filter(rows[:, 3], w, 2)
    temp = rows[:, 5]
    on = (rows[:, 1] > 99.0) & (rows[:, 2] > 0.0) & (dpr - dzero[spell] > 3.0)
    on &= np.convolve(on, np.ones(2 * w + 1), mode="same") >= 2 * w + 1      # away from the edges of each stretch
    keep = np.where(on)[0][::2]
    runs = np.unique(run[keep]) if stretch_zero else np.array([0])
    run_of = np.searchsorted(runs, run[keep]) if stretch_zero else np.zeros(len(keep), dtype=int)
    pv = vapour(temp[keep])
    a = {"q": q[keep], "dq": dq[keep], "g": density(temp[keep]) / 999.10, "fpv": (0.96 - 0.28 * np.sqrt(pv / 22064.0)) * pv}
    d, big = sheet["valve_end_bore_mm"], sheet["loop_pipe"]["inside_diameter_mm"]
    qmax = float(a["q"].max())
    nr = len(runs)

    def state(e):
        p, dpt = tap_state(e, p1r[keep], dpr[keep], spell[keep], dzero, pread, patm, boyle)
        return dict(a, p1=p, dp=np.maximum(dpt, 0.05))

    a0 = state(0.0)
    k0 = float(np.percentile(a["q"] * np.sqrt(a["g"] / (a0["dp"] / 100.0)), 90))
    lo = [1e-3, 0.4, math.log(1.2) if "n" in free else math.log(399.0), -80.0 if "e" in free else -1e-6, 0.0]
    hi = [1e4, 0.995, math.log(400.0), 80.0 if "e" in free else 1e-6, 120.0 if tau_free else 1e-6]
    lo += [-0.15 * qmax if "q0" in free else -1e-6] * nr
    hi += [0.15 * qmax if "q0" in free else 1e-6] * nr

    def resid(x):
        st = state(x[3])
        return (model(np.r_[x[:4], 0.0], d, big, st) + x[5:][run_of] - a["q"] - x[4] * a["dq"]) / qmax

    best = None
    for f0 in (0.7, 0.85, 0.95):
        for n0 in ((3.0, 10.0, 60.0) if "n" in free else (399.5,)):
            for tau0 in ((2.0, 30.0) if tau_free else (0.0,)):
                x0 = np.clip([k0 * 1.1, f0, math.log(n0), 0.0, tau0] + [0.0] * nr, np.array(lo) + 1e-9, np.array(hi) - 1e-9)
                r = least_squares(resid, x0, bounds=(lo, hi), x_scale=[k0, 0.1, 1.0, 10.0, 5.0] + [0.01 * qmax] * nr)
                if best is None or r.cost < best.cost:
                    best = r
    return float(best.x[0]), float(best.x[1]), best.x


def reduce(folder, free=("n", "e", "q0"), damping=True, samples=True, boyle=True, stretch_zero=True):
    """the reference fits every rated-travel sample with the converter's damping; routes of the evidence:
    samples=True, damping=False the same fit with the shown flow taken as it is; samples=False the hold means
    (damping=False: the lab's convention of settled holds and last 45 s, the shown flow as it is; damping=True: every
    stretch between moves, each window mean lifted by tau times the shown flow's change over the window)"""
    sheet, rows = load(folder)
    if samples:
        kv, fl, _ = sample_fit(sheet, rows, damping, free, boyle, stretch_zero)
        return {"kv_rated_m3h": kv, "fl_rated": fl}
    pts, _ = corrected_holds(sheet, rows, damping)
    pts = [p for p in pts if p["travel"] > 99.0]
    kv, fl, _ = fit(pts, sheet["valve_end_bore_mm"], sheet["loop_pipe"]["inside_diameter_mm"], free)
    return {"kv_rated_m3h": kv, "fl_rated": fl}


def main():
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)


if __name__ == "__main__":
    main()
