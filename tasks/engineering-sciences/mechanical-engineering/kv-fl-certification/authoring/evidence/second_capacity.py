"""Second reduction, written apart from the reference: converter damping from the meter's run-downs, 1 s blocks, and a
Nelder-Mead fit.

    python3 second_capacity.py <specimen folder> <answer.json>

- dp zero and upstream column: closed-valve spells (travel below 1 %, 20 s or longer) whose dp lies within 20 kPa of
  the first spell's are pump-off spells (the others are valves shut with the pump running); the mean dp and p1 of
  each. The first spell's dp is the transmitter's zero; the change of a later spell's dp and p1 from the first is air
  taken into the dp high leg and the p1 line, read at the pump-off tap pressure; with the pump running it shrinks in
  the ratio of the absolute tap pressures (Boyle), the tap pressure solved inside the fit with the p1 correction.
- Damping: when the pump stops, the flow through the meter ends within a second, and the converter's figure runs
  down as A exp(-t / tau) + c until the cut-off; tau is the median over the run-downs long enough to fit (four samples
  or more) of a three-parameter fit (scipy curve_fit) to each. Below 5 s (or with no run-down long enough) the
  run-downs are too short to read tau from, and tau is refined in the fit below from that start.
- 1 s blocks: block means of every column; the flow through the meter is the shown flow plus tau times its central
  difference over two blocks.
- Fit: rated-travel blocks with the pump on (15 s away from any start or stop), averaged in chunks of ten consecutive
  blocks; Nelder-Mead on log Kv, FL, log n, the
  p1 correction (and tau when refined), the meter zero of each run (between closed-valve spells) solved by weighted least squares at every step,
  minimising the squared relative flow misfit of the installed relations with the
  trim transition, from nine starts.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit, minimize


def rho(t):
    return 999.97495 * (1 - (t - 3.983035) ** 2 * (t + 301.797) / (522528.9 * (t + 69.34881)))


def reduce(folder):
    folder = Path(folder)
    sh = json.loads((folder / "specimen_sheet.json").read_text())
    with open(folder / "daq_record.csv", newline="") as fh:
        a = np.array([[float(r[k]) for k in ("t_s", "travel_pct", "flow_m3h", "p1_gauge_kPa", "dp_kPa", "water_C")]
                      for r in csv.DictReader(fh)])
    hz = sh["daq"]["rate_Hz"]
    shut = a[:, 1] < 1.0
    spell = np.full(len(a), -1)
    runid = np.zeros(len(a), dtype=int)
    zs, cs = [], []
    i = 0
    while i < len(a):
        if shut[i]:
            j = i
            while j + 1 < len(a) and shut[j + 1]:
                j += 1
            if a[j, 0] - a[i, 0] >= 20.0:
                runid[i:] += 1                       # the flow through the meter stops here
                zm = float(a[i + int(5 * hz):j + 1, 4].mean())
                if not zs or abs(zm - zs[0]) < 20.0:     # pump off (else the valve is shut with the pump running)
                    zs.append(zm)
                    cs.append(float(a[i + int(5 * hz):j + 1, 3].mean()))
                    spell[i:] = len(zs) - 1
            i = j + 1
        else:
            i += 1
    spell[spell < 0] = 0
    zs, cs = np.array(zs), np.array(cs)
    # run-downs of the damped figure after each stop
    taus = []
    stops = np.where(shut[1:] & ~shut[:-1])[0] + 1
    for s in stops:
        j = s
        while j < len(a) and a[j, 2] > 0.0 and shut[j]:
            j += 1
        if j - s < 4:                      # too short to fit: the figure was below the cut-off almost at once
            continue
        tt, qq = a[s:j, 0] - a[s, 0], a[s:j, 2]
        try:
            prm, _ = curve_fit(lambda x, amp, tau, c: amp * np.exp(-x / tau) + c, tt, qq,
                               p0=(qq[0], max(tt[-1] / 2.0, 0.5), 0.0), maxfev=5000)
            taus.append(float(abs(prm[1])))
        except RuntimeError:
            continue
    tau = float(np.median(taus)) if taus else 0.0
    refine = tau < 5.0          # run-downs of a lightly damped converter are too short to read tau from
    nb = len(a) // int(hz)
    blk = a[: nb * int(hz)].reshape(nb, int(hz), 6).mean(axis=1)
    sb = spell[: nb * int(hz)].reshape(nb, int(hz)).max(axis=1)
    rb = runid[: nb * int(hz)].reshape(nb, int(hz)).max(axis=1)
    on = (blk[:, 1] > 99.0) & (a[: nb * int(hz), 2].reshape(nb, int(hz)).min(axis=1) > 0.0)
    near = np.convolve(~on, np.ones(31), mode="same") > 0
    keep = np.where(on & ~near)[0]
    keep = keep[(keep > 0) & (keep < nb - 1)]
    qshow = blk[keep, 2]
    qrate = (blk[keep + 1, 2] - blk[keep - 1, 2]) / 2.0
    qtrue = qshow + tau * qrate
    dp = blk[keep, 4] - zs[sb[keep]]
    ok = dp > 3.0
    keep, qtrue, dp, qshow, qrate = keep[ok], qtrue[ok], dp[ok], qshow[ok], qrate[ok]
    # 10-block chunks of consecutive kept blocks: means of the deconvolved flow and of every column
    breaks = np.where(np.diff(keep) > 1)[0] + 1
    groups = [c for run in np.split(np.arange(len(keep)), breaks) for c in np.array_split(run, max(len(run) // 10, 1))]
    qtrue = np.array([qtrue[c].mean() for c in groups])
    qshow = np.array([qshow[c].mean() for c in groups])
    qrate = np.array([qrate[c].mean() for c in groups])
    chunk = np.array([blk[keep[c]].mean(axis=0) for c in groups])
    sp = np.array([sb[keep[c][0]] for c in groups])
    rp = np.array([rb[keep[c][0]] for c in groups])
    run = np.searchsorted(np.unique(rp), rp)
    big = sh["loop_pipe"]["inside_diameter_mm"]
    d = sh["valve_end_bore_mm"]
    patm = 101.325 * (1 - 2.25577e-5 * sh["lab"]["elevation_m"]) ** 5.25588
    reading = chunk[:, 3] + patm
    dpread = chunk[:, 4]
    dc, dz = cs[sp] - cs[0], zs[sp] - zs[0]

    def taps(e):
        """the air taken into the p1 line and the dp high leg since the first spell, read at the pump-off tap pressure,
        shrinks as p_off / p with the tap pressure (Boyle): p from p^2 - (R + e) p + dc p_off = 0"""
        p_off = cs[0] + patm + e
        s_ = reading + e
        p = 0.5 * (s_ + np.sqrt(np.clip(s_ * s_ - 4 * dc * p_off, 1e-6, None)))
        return p, np.clip(dpread - zs[0] - dz * p_off / p, 0.05, None)
    g = rho(chunk[:, 5]) / 999.10
    pv = 0.133322 * 10 ** (8.07131 - 1730.63 / (233.426 + chunk[:, 5]))
    fpv = (0.96 - 0.28 * np.sqrt(pv / 22064)) * pv
    q = qtrue
    r2 = (d / big) ** 2
    z12, z1b = 1.5 * (1 - r2) ** 2, 0.5 * (1 - r2) ** 2 + 1 - r2 ** 2

    p1, dp = taps(0.0)
    x = np.sqrt(dp / 100 / g)
    low = dp < np.percentile(dp, 30)
    slope, q00 = np.polyfit(x[low], q[low], 1)
    e0 = 0.0

    def flow(v, q):
        kv, fl, nn, e = math.exp(v[0]), v[1], math.exp(v[2]), v[3]
        p1, dp = taps(e)
        load = (kv / d ** 2) ** 2 / 0.0016
        fp = 1 / math.sqrt(1 + z12 * load)
        flp = fl / math.sqrt(1 + fl * fl * z1b * load)
        qu = fp * kv * np.sqrt(dp / 100 / g)
        qc = flp * kv * np.sqrt(np.clip(p1 - fpv, 1e-6, None) / 100 / g)
        f = qu * (1 + (qu / qc) ** nn) ** (-1 / nn)
        # the meter's zero of each pump run (a new one after every stop) by weighted least squares at these values
        w = 1.0 / q ** 2
        q0 = np.bincount(run, w * (q - f)) / np.bincount(run, w)
        return f + q0[run]

    def loss(v):
        if not (0.3 < v[1] < 1.0) or not (0.1 < v[2] < 6.0) or not (0.0 <= v[4] < 120.0):
            return 1e9
        qt = qshow + (v[4] if refine else tau) * qrate     # the converter's lag
        return float(np.sum((flow(v, qt) / qt - 1) ** 2))

    best = None
    for f0 in (0.7, 0.85, 0.95):
        for n0 in (3.0, 8.0, 40.0):
            v0 = [math.log(max(slope, 1e-3) * 1.02), f0, math.log(n0), e0, tau]
            r = minimize(loss, v0, method="Nelder-Mead",
                         options={"xatol": 1e-9, "fatol": 1e-14, "maxiter": 6000, "maxfev": 6000, "adaptive": True})
            r = minimize(loss, r.x, method="Nelder-Mead",
                         options={"xatol": 1e-9, "fatol": 1e-14, "maxiter": 6000, "maxfev": 6000, "adaptive": True})
            if best is None or r.fun < best.fun:
                best = r
    return {"kv_rated_m3h": float(math.exp(best.x[0])), "fl_rated": float(best.x[1]), "tau_s": float(best.x[4]) if refine else tau}


if __name__ == "__main__":
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)
