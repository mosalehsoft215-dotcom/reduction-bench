"""Reduce one coupon session to the film's Gilbert damping, spectroscopic g and mu0 Meff.

    python3 fmr_reduce.py <coupon folder> <answer json>

1. Every sweep is fitted in gaussmeter-reading units with the full detection chain: a complex Lorentzian line
   (free phase), the first harmonic of the logged field modulation, a linear background, and the output filter's
   logged sections and time constant applied along the sweep in time (rate and direction from the register).
2. Gaussmeter: reading = s25 (1 + beta (T - 25)) B + o with beta from the method note; the six marker resonances
   (g 2.0036), each at the probe temperature of its resonance point, give s25 and o by least squares; each film
   resonance and width is converted at its own probe temperature (linear in time within the sweep).
3. Geometry: g, Meff and the tilt of the normal mount from all resonance fields at once (in-plane at 90 deg),
   with the Smit-Beljers resonance of an easy-plane film.
4. Damping: normal-mount half widths = alpha * Gilbert factor + dM * |dB_res / dMeff|, both factors from the fitted
   geometry; the in-plane widths are not used (two-magnon scattering adds to them in plane).
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lfilter

MU_B = 9.2740100783e-24
HBAR = 1.054571817e-34
H_PLANCK = 6.62607015e-34
G_MARKER = 2.0036
BETA = -4.5e-4
OVER = 4


# ---------------------------------------------------------------- physics
def w_field(freq, g):
    return 2 * math.pi * freq * HBAR / (g * MU_B)


def _branch(theta, th_b, meff):
    """Field, B1 and B2 on the saturated branch, parametrised by the magnetisation angle theta (> th_b)."""
    B = meff * np.sin(2 * theta) / (2 * np.sin(theta - th_b))
    c = np.cos(theta - th_b)
    return B, B * c - meff * np.cos(2 * theta), B * c - meff * np.cos(theta) ** 2


DELTA = np.geomspace(1e-12, 0.9, 3000)


def _theta_res(freq, g, meff, th_b):
    th_b = max(th_b, 1e-7)
    w2 = w_field(freq, g) ** 2
    B, b1, b2 = _branch(th_b + DELTA, th_b, meff)
    F = b1 * b2
    below = np.nonzero(F < w2)[0]
    if len(below) == 0 or below[0] == 0:      # no saturated-branch resonance: take the stiffness minimum
        return th_b + DELTA[int(np.argmin(F))], th_b
    k = int(below[0])
    lo, hi = DELTA[k - 1], DELTA[k]          # F(lo) >= w2 > F(hi)
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        _, c1, c2 = _branch(th_b + mid, th_b, meff)
        lo, hi = (mid, hi) if c1 * c2 >= w2 else (lo, mid)
    return th_b + math.sqrt(lo * hi), th_b


def b_res(freq, g, meff, th_b):
    if th_b >= math.pi / 2 - 1e-12:          # in plane: f = gamma/2pi sqrt(B (B + Meff))
        w = w_field(freq, g)
        return 0.5 * (-meff + math.sqrt(meff * meff + 4 * w * w))
    th, tb = _theta_res(freq, g, meff, th_b)
    return float(_branch(th, tb, meff)[0])


def width_factors(freq, g, meff, th_b):
    """Gilbert factor (half width per unit alpha) and |dB_res/dMeff| on the normal mount."""
    th, tb = _theta_res(freq, g, meff, th_b)
    eps = 1e-6 * (th - tb)
    Bp, p1, p2 = _branch(th + eps, tb, meff)
    Bm, m1, m2 = _branch(th - eps, tb, meff)
    B0, c1, c2 = _branch(th, tb, meff)
    dF_dB = (p1 * p2 - m1 * m2) / (Bp - Bm)
    gil = w_field(freq, g) * (c1 + c2) / abs(dF_dB)
    dm = 2e-5
    inh = abs(b_res(freq, g, meff + dm, th_b) - b_res(freq, g, meff - dm, th_b)) / (2 * dm)
    return float(gil), float(inh)


def ip_factors(freq, g, meff):
    """In-plane Gilbert factor and |dB_res/dMeff| (used only by the wrong-route switches)."""
    w = w_field(freq, g)
    br = b_res(freq, g, meff, math.pi / 2)
    return w * (2 * br + meff) / (2 * br + meff), br / (2 * br + meff)


# ---------------------------------------------------------------- one sweep
def harmonic(x, br, hw, a):
    u = br - x - 1j * hw
    if a <= 0:
        return 1.0 / u ** 2
    return (2.0 / a) * (u / (np.sqrt(u - a) * np.sqrt(u + a)) - 1.0)


def sections(y, step_s, tau, poles):
    if tau <= 0:
        return y
    p = math.exp(-step_s / tau)
    for _ in range(poles):
        y = lfilter([1 - p], [1, -p], y, zi=[p * y[0]])[0]
    return y


def fit_sweep(x, y, row, filt=True, mod=True):
    a = row["mod_mT"] * 1e-3 if mod else 0.0
    tau = row["tau_ms"] * 1e-3 if filt else 0.0
    poles = int(round(row["slope_dB_oct"] / 6))
    rate = row["rate_mT_s"] * 1e-3
    up = row["direction"] == "up"
    fine = np.linspace(x[0], x[-1], (len(x) - 1) * OVER + 1)
    step_s = abs(fine[1] - fine[0]) / rate
    span = x[-1] - x[0]
    xm = 0.5 * (x[0] + x[-1])

    def shape(br, hw, ph):
        core = np.real(np.exp(1j * ph) * harmonic(fine, br, hw, a))
        if up:
            return sections(core, step_s, tau, poles)[::OVER]
        return sections(core[::-1], step_s, tau, poles)[::-1][::OVER]

    def resid(p):
        br, hw, ph, amp, c0, c1 = p
        return amp * shape(br, hw, ph) + c0 + c1 * (x - xm) / span - y

    i_max, i_min = int(np.argmax(y)), int(np.argmin(y))
    br0 = 0.5 * (x[i_max] + x[i_min])
    hw0 = max(abs(x[i_max] - x[i_min]) * math.sqrt(3) / 2, 2 * abs(x[1] - x[0]))
    best = None
    for ph0 in (-1.2, -0.4, 0.4, 1.2, 2.0, -2.0):
        s0 = shape(br0, hw0, ph0)
        amp0 = float(np.dot(s0, y - y.mean()) / max(np.dot(s0, s0), 1e-30))
        p0 = [br0, hw0, ph0, amp0, float(np.median(y)), 0.0]
        try:
            r = least_squares(resid, p0, x_scale=[hw0, hw0, 1, abs(amp0) + 1e-9, 1, 1], max_nfev=400)
        except Exception:
            continue
        if r.x[1] > 0 and (best is None or r.cost < best.cost):
            best = r
    if best is None:
        raise RuntimeError("sweep fit failed")
    br, hw = best.x[0], abs(best.x[1])
    # probe temperature at the moment the resonance was passed (linear in time within the sweep)
    frac = (br - x[0]) / span
    if not up:
        frac = 1 - frac
    frac = min(max(frac, 0.0), 1.0)
    temp = row["probe_C_start"] + (row["probe_C_end"] - row["probe_C_start"]) * frac
    return br, hw, temp


# ---------------------------------------------------------------- session
def load(folder):
    folder = Path(folder)
    rows = []
    with open(folder / "register.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            for k in ("freq_GHz", "mod_mT", "tau_ms", "slope_dB_oct", "rate_mT_s", "start_min",
                      "probe_C_start", "probe_C_end"):
                r[k] = float(r[k])
            arr = np.loadtxt(folder / "sweeps" / f"{r['sweep']}.csv", delimiter=",", skiprows=1)
            r["x"], r["y"] = arr[:, 0] * 1e-3, arr[:, 1]
            rows.append(r)
    return rows


def reduce(folder, filt=True, mod=True, temp=True, tilt_free=True, widths="normal", gil_aligned=False):
    """Keyword switches exist for the wrong-route record (authoring/evidence); the defaults are the reduction."""
    rows = load(folder)
    beta = BETA if temp else 0.0
    for r in rows:
        r["br_read"], r["hw_read"], r["temp"] = fit_sweep(r["x"], r["y"], r, filt, mod)

    marks = [r for r in rows if r["geometry"] == "marker"]
    A = np.array([[(1 + beta * (r["temp"] - 25)) * H_PLANCK * r["freq_GHz"] * 1e9 / (G_MARKER * MU_B), 1.0]
                  for r in marks])
    s25, off = np.linalg.lstsq(A, np.array([r["br_read"] for r in marks]), rcond=None)[0]

    film = [r for r in rows if r["kind"] == "film"]
    for r in film:
        s = s25 * (1 + beta * (r["temp"] - 25))
        r["br"], r["hw"] = (r["br_read"] - off) / s, r["hw_read"] / s
        r["th"] = math.pi / 2 if r["geometry"] == "in-plane" else None

    ip = [r for r in film if r["geometry"] == "in-plane"]
    nm = [r for r in film if r["geometry"] == "normal"]

    def geo_res(p):
        g, meff, tilt = p if tilt_free else (p[0], p[1], 0.0)
        out = [b_res(r["freq_GHz"] * 1e9, g, meff, math.pi / 2) - r["br"] for r in ip]
        out += [b_res(r["freq_GHz"] * 1e9, g, meff, abs(tilt)) - r["br"] for r in nm]
        return np.array(out) * 1e4

    # start from the in-plane relation f = gamma/2pi sqrt(B (B + Meff)) and the normal one
    fs = np.array([r["freq_GHz"] * 1e9 for r in nm])
    sl, ic = np.polyfit([r["br"] for r in nm], fs, 1)
    g0 = sl * 2 * math.pi * HBAR / MU_B
    p0 = [min(max(g0, 1.5), 2.6), min(max(-ic / sl, 0.05), 2.9), 0.004]
    best = None
    for t0 in (0.002, 0.008, 0.015):
        p0[2] = t0
        r = least_squares(geo_res, p0, x_scale=[0.01, 0.01, 0.005], diff_step=1e-6,
                          bounds=([1.0, 0.01, 0.0], [3.0, 3.0, 0.035]))
        if best is None or r.cost < best.cost:
            best = r
    g, meff, tilt = best.x
    tilt = abs(tilt)

    if not tilt_free:
        tilt = 0.0
    use = nm if widths == "normal" else (ip if widths == "in-plane" else ip + nm)
    G = np.array([width_factors(r["freq_GHz"] * 1e9, g, meff, tilt) if r["geometry"] == "normal"
                  else ip_factors(r["freq_GHz"] * 1e9, g, meff) for r in use])
    if gil_aligned:              # wrong route: Gilbert width taken as alpha * omega / gamma whatever the angle
        G[:, 0] = [w_field(r["freq_GHz"] * 1e9, g) for r in use]
    alpha, dM = np.linalg.lstsq(G, np.array([r["hw"] for r in use]), rcond=None)[0]
    info = {"s25": s25, "offset_T": off, "tilt_deg": math.degrees(tilt), "dM_T": dM}
    return {"damping_alpha": float(alpha), "spectroscopic_g": float(g), "mu0_meff_T": float(meff)}, info


def main():
    res, _ = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
