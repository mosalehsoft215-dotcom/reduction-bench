"""Second reduction of one secondary-injection session, written apart from the reference.

    python3 second_found.py <relay folder> <answer.json>

- Cycle values: each 20-sample cycle (counted from the first sample) gets a least-squares 50 Hz sine and cosine plus
  a constant; the sine amplitude over sqrt(2) is its rms.
- Shots: from the table's start command, inception is the first sample above 10 % of the shot's peak and the trip the
  first touch of the binary input; the shot current is the 50 Hz amplitude of a least-squares sine, cosine and
  decaying exponential over inception to trip, the exponential's time constant taken from a grid of seven values
  (15-50 ms) by the smallest residual.
- Stop-input delay = mean high-set operate time minus 25 ms.
- Element: for each family and for instantaneous and dependent reset, Is, TMS and the definite minimum time are found
  by Nelder-Mead on the squared log-time misfit of the shots of 40 ms or longer (shorter ones still charge the
  element) from four starts.  Between shots the element is stepped once per cycle
  with that cycle's M (above pickup up, below it down by (1 - M^2) / (TMS tr) per second or to zero at once); the
  floor at zero is taken as G_n = S_n - min(0, min over k <= n of S_k) on the running sum S.  The smallest misfit wins.
- CT ratio from the primary-injection check by ordinary least squares with intercept (the intercept is dropped).
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

FAMILIES = {"SI": (0.14, 0.02, 13.5), "VI": (13.5, 1.0, 47.3), "EI": (80.0, 2.0, 80.0)}


def read(folder):
    folder = Path(folder)
    z = np.load(folder / "session_record.npz")
    fs = float(z["rate_hz"])
    x = z["current_counts"] * float(z["amps_per_count"])
    b = z["binary"]
    n = len(x) // 20
    tt = np.arange(20) / fs
    basis = np.column_stack([np.sin(2 * math.pi * 50 * tt), np.cos(2 * math.pi * 50 * tt), np.ones(20)])
    coef = np.linalg.lstsq(basis, x[:n * 20].reshape(n, 20).T, rcond=None)[0]
    cyc = np.hypot(coef[0], coef[1]) / math.sqrt(2)
    with open(folder / "test_set_table.csv", newline="") as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r["shot"]))
    shots = []
    for r in rows:
        s0 = int(round(float(r["start_clock_s"]) * fs))
        trip = s0 + int(np.argmax(b[s0:] == 1))
        seg = x[s0:trip]
        inc = s0 + int(np.argmax(np.abs(seg) > 0.1 * np.abs(seg).max()))
        tt = np.arange(trip - inc) / fs
        fits = []
        for tau in (0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05):
            lin = np.column_stack([np.sin(2 * math.pi * 50 * tt), np.cos(2 * math.pi * 50 * tt), np.exp(-tt / tau)])
            c, res = np.linalg.lstsq(lin, x[inc:trip], rcond=None)[:2]
            fits.append((float(res[0]) if len(res) else 0.0, float(np.hypot(c[0], c[1]))))
        amp = min(fits)[1]
        shots.append((r["element"], s0 // 20, trip // 20 + 1, amp / math.sqrt(2), (trip - inc) / fs))
    return cyc, shots


def predicted(cyc, shots, delay, fam, timed, isec, tms, tmin):
    k, a, tr = FAMILIES[fam]
    g, last, out = 0.0, 0, []
    for el, c0, c1, i, t in shots:
        m = cyc[last:c0] / isec
        if len(m):
            step = np.where(m > 1, (np.clip(m, 1, None) ** a - 1) / (k * tms), -(1 - m * m) / (tms * tr)) * 0.02
            if timed:
                s = g + np.cumsum(step)
                g = float(s[-1] - min(0.0, np.minimum.accumulate(s)[-1]))
            else:
                low = np.nonzero(m <= 1)[0]
                g = float(np.sum(step[low[-1] + 1:])) if len(low) else g + float(np.sum(step))
            g = min(g, 1.0)
        rate = max((i / isec) ** a - 1, 0.0) / (k * tms)
        if el == "inverse-time":
            out.append(max((1 - g) / max(rate, 1e-12), tmin))
            g = min(1.0, g + rate * (t - delay))
        else:
            g = min(1.0, g + rate * t)
        last = c1
    return np.array(out)


def reduce(folder):
    sheet = json.loads((Path(folder) / "relay_sheet.json").read_text())
    study = float(sheet["coordination_study"]["pickup_secondary_A"])
    cyc, shots = read(folder)
    hs = [s[4] for s in shots if s[0] == "high-set"]
    delay = max(0.0, float(np.mean(hs)) - 0.025) if hs else 0.0
    t_obs = np.array([s[4] for s in shots if s[0] == "inverse-time"]) - delay
    use = t_obs >= 0.04                                     # shots under two cycles: current not separable from the offset
    i_lo = min(s[3] for s in shots if s[0] == "inverse-time")
    best = None
    for fam in FAMILIES:
        for timed in (False, True):
            def J(z):
                isec, tms, tmin = math.exp(z[0]), math.exp(z[1]), z[2]
                if isec >= i_lo or tmin < 0 or tms > 5:
                    return 1e9
                p = predicted(cyc, shots, delay, fam, timed, isec, tms, tmin)
                return float(np.sum(((np.log(np.maximum(p, 1e-6)) - np.log(np.maximum(t_obs, 1e-4))) ** 2)[use]))
            for f in (0.45, 0.85):
                for tm in (0.2,):
                    for t0 in (0.0, 0.9 * float(t_obs[use].min())):
                        if study * f >= i_lo:
                            continue
                        r = minimize(J, [math.log(study * f), math.log(tm), t0], method="Nelder-Mead",
                                     options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 400})
                        if best is None or r.fun < best[0]:
                            best = (r.fun, r.x)
    with open(Path(folder) / "primary_check.csv", newline="") as fh:
        pc = np.array([[float(r["primary_A"]), float(r["relay_display_secondary_A"])] for r in csv.DictReader(fh)])
    slope = np.polyfit(pc[:, 1], pc[:, 0], 1)[0]
    return {"pickup_primary_A": float(math.exp(best[1][0]) * slope), "tms_found": float(math.exp(best[1][1]))}


def main(folder, out):
    with open(out, "w") as fh:
        json.dump(reduce(folder), fh)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
