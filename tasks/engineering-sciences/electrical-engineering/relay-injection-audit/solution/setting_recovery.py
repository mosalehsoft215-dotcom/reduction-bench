"""Reference reduction of one secondary-injection session: pickup of the inverse-time element in primary amperes and TMS.

    python3 setting_recovery.py <relay folder> <answer.json>

1. Session record (1 kHz from the first sample, cycles of 20 samples counted from there).  The 50 Hz rms of every
   cycle comes from a one-cycle discrete Fourier sum.  Each fault shot starts at its start command (test-set table,
   on the session clock); its current inception is the first sample above 3 % of the shot's peak and its trip the
   first touch of the binary input.  The fault current is the fundamental of a least-squares fit over inception to
   trip of a 50 Hz sine and cosine with a decaying DC offset (time constant searched), the third harmonic and a
   constant.  Every cycle outside the shots, load or breaker-open, keeps its own rms.
2. Stop-input delay: the median of the high-set shots' operate times minus 25 ms; the relay's operate time of an
   inverse-time shot is its inception-to-touch time minus that delay.
3. Replay of the inverse-time element over the whole session: from G = 0 at the first sample, each cycle outside the
   shots moves G by its own M = I / Is, (M^alpha - 1) / (k TMS) x 20 ms above pickup, and below pickup either
   resets it (instantaneous) or lowers it by (1 - M^2) / (TMS tr) x 20 ms (dependent reset), never below 0; the
   clamp at 0 is applied to whole stretches at once through the running minimum of the cumulative sum.  In a shot G
   grows linearly from its value at the start, the predicted operate time being max((1 - G0) / rate, t_min); after
   the shot G is what the observed shot time left.
4. For each family (SI, VI, EI) and each reset mode, Is (secondary), TMS and t_min are fitted on log operate time by
   least squares from eight starts (Is from 0.35 to 1.1 of the coordination study's pickup); the combination with
   the smallest residual is the one the relay runs.
5. CT ratio in service from the primary-injection check: slope of primary amperes on the relay's displayed secondary
   amperes through the origin.  Pickup in primary amperes = Is x that ratio.
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
from scipy.optimize import least_squares  # noqa: E402

CURVES = {"SI": (0.14, 0.02), "VI": (13.5, 1.0), "EI": (80.0, 2.0)}
RESET_TR = {"SI": 13.5, "VI": 47.3, "EI": 80.0}
HIGHSET_S = 0.025
CYC = 20
DT = 0.02


def fundamental(x, fs):
    """rms of the 50 Hz component between inception and trip: least squares of a sine and cosine at 50 Hz, a decaying
    DC offset exp(-t / tau) (tau by golden-section search, 5-120 ms), the third harmonic and a constant"""
    t = np.arange(len(x)) / fs
    w = 2 * math.pi * 50.0
    base = np.column_stack([np.sin(w * t), np.cos(w * t), np.sin(3 * w * t), np.cos(3 * w * t), np.ones_like(t)])

    def solve(tau):
        m = np.column_stack([base, np.exp(-t / tau)])
        c, *_ = np.linalg.lstsq(m, x, rcond=None)
        return float(np.sum((m @ c - x) ** 2)), c

    lo, hi = math.log(0.005), math.log(0.12)
    g = (math.sqrt(5) - 1) / 2
    for _ in range(28):
        a, b = hi - g * (hi - lo), lo + g * (hi - lo)
        if solve(math.exp(a))[0] < solve(math.exp(b))[0]:
            hi = b
        else:
            lo = a
    c = solve(math.exp((lo + hi) / 2))[1]
    return math.hypot(c[0], c[1]) / math.sqrt(2.0)


def session(folder):
    """cycle rms of the whole session, the shots (inception, trip, current, element) and the cycles each shot covers"""
    folder = Path(folder)
    z = np.load(folder / "session_record.npz")
    fs = float(z["rate_hz"])
    cur = z["current_counts"].astype(float) * float(z["amps_per_count"])
    b = z["binary"]
    ncyc = len(cur) // CYC
    blocks = cur[:ncyc * CYC].reshape(ncyc, CYC)
    ph = np.exp(-2j * math.pi * np.arange(CYC) / CYC)
    rms = np.abs(blocks @ ph) * math.sqrt(2.0) / CYC
    with open(folder / "test_set_table.csv", newline="") as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r["shot"]))
    shots = []
    for r in rows:
        s0 = int(round(float(r["start_clock_s"]) * fs))
        touch = s0 + int(np.argmax(b[s0:] > 0))
        seg = cur[s0:touch]
        inc = s0 + int(np.argmax(np.abs(seg) > 0.03 * np.abs(seg).max()))
        shots.append({"element": r["element"], "c0": s0 // CYC, "c1": touch // CYC + 1,
                      "i": fundamental(cur[inc:touch], fs), "t": (touch - inc) / fs})
    return rms, shots


def load_stretches(rms, shots):
    """the cycle rms between consecutive shots (and before the first)"""
    out, c = [], 0
    for s in shots:
        out.append(rms[c:s["c0"]])
        c = s["c1"]
    return out


def drift(m, a, k, tms, tr, timed, g0, mode):
    """G after a stretch of cycles with rms ratios m, from g0"""
    if len(m) == 0:
        return g0
    if mode == "static":
        return 0.0
    if mode == "zero":
        return max(0.0, g0 - len(m) * DT / (tms * tr)) if timed else 0.0
    if mode == "mean":
        mm = float(np.mean(m))
        span = len(m) * DT
        if mm > 1.0:
            return g0 + (mm ** a - 1.0) / (k * tms) * span
        return max(0.0, g0 - (1.0 - mm * mm) / (tms * tr) * span) if timed else 0.0
    up = np.where(m > 1.0, (np.maximum(m, 1.0) ** a - 1.0) / (k * tms), 0.0) * DT
    if not timed:
        below = np.flatnonzero(m <= 1.0)
        return float(up[below[-1] + 1:].sum()) if len(below) else g0 + float(up.sum())
    d = up - np.where(m <= 1.0, (1.0 - m * m) / (tms * tr), 0.0) * DT
    s = g0 + np.cumsum(d)
    return float(s[-1] - min(0.0, s.min()))


def replay(rms, shots, stretches, isec, tms, tmin, curve, timed, mode, delay):
    k, a = CURVES[curve]
    tr = RESET_TR[curve]
    g, pred = 0.0, []
    for s, st in zip(shots, stretches):
        g = min(1.0, drift(st / isec, a, k, tms, tr, timed, g, mode))
        rate = max((s["i"] / isec) ** a - 1.0, 0.0) / (k * tms)
        if s["element"] == "inverse-time":
            pred.append(max((1.0 - g) / max(rate, 1e-12), tmin))
            g = min(1.0, g + rate * (s["t"] - delay))
        else:
            g = min(1.0, g + rate * s["t"])
    return np.array(pred)


def fit(rms, shots, delay, curve, timed, mode, study):
    stretches = load_stretches(rms, shots)
    inv = [s for s in shots if s["element"] == "inverse-time"]
    t_obs = np.maximum(np.array([s["t"] for s in inv]) - delay, 1e-4)
    i_lo = min(s["i"] for s in inv)

    def resid(x):
        p = replay(rms, shots, stretches, math.exp(x[0]), math.exp(x[1]), x[2], curve, timed, mode, delay)
        return np.log(np.maximum(p, 1e-6)) - np.log(t_obs)

    best = None
    for frac in (0.35, 0.6, 0.9, 1.1):
        for tms0 in (0.2,):
            for tmin0 in (0.0, 0.97 * float(np.min(t_obs))):
                x0 = [math.log(min(study * frac, i_lo * 0.95)), math.log(tms0), tmin0]
                try:
                    r = least_squares(resid, x0, bounds=([math.log(1e-4), math.log(1e-3), 0.0],
                                                         [math.log(i_lo * 0.999), math.log(5.0), 10.0]))
                except ValueError:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    return math.exp(best.x[0]), math.exp(best.x[1]), best.x[2], best.cost


def ct_ratio(folder):
    with open(Path(folder) / "primary_check.csv", newline="") as fh:
        pairs = np.array([[float(r["primary_A"]), float(r["relay_display_secondary_A"])] for r in csv.DictReader(fh)])
    return float(np.dot(pairs[:, 0], pairs[:, 1]) / np.dot(pairs[:, 1], pairs[:, 1]))


def stop_delay(shots):
    hs = [s["t"] for s in shots if s["element"] == "high-set"]
    return max(0.0, float(np.median(hs)) - HIGHSET_S) if hs else 0.0


def reduce(folder, mode="replay", modes=(False, True)):
    """mode: "replay" (the reference); routes for authoring/evidence: "mean" (each stretch between shots at its mean
    current), "zero" (no current between shots), "static" (every shot from G = 0)"""
    sheet = json.loads((Path(folder) / "relay_sheet.json").read_text())
    study = float(sheet["coordination_study"]["pickup_secondary_A"])
    rms, shots = session(folder)
    delay = stop_delay(shots)
    fits = {(c, m): fit(rms, shots, delay, c, m, mode, study) for c in CURVES for m in modes}
    isec, tms, _, _ = min(fits.values(), key=lambda f: f[3])
    return {"pickup_primary_A": float(isec * ct_ratio(folder)), "tms_found": float(tms)}


def main():
    with open(sys.argv[2], "w") as fh:
        json.dump(reduce(sys.argv[1]), fh)


if __name__ == "__main__":
    main()
