"""Independent second measurement of shear-wave splitting: the eigenvalue method on a coarse-to-fine grid.

    python3 eigen_split.py <station folder> <answer file>

- Response: the calibration step's zero crossing and first trough give each channel's T0 and h (Nelder-Mead on the
  shape, gain from a linear fit); each spectrum is divided by the response from scipy.signal.freqs, with a 1 % water
  level, then band-limited to 0.03-0.8 Hz.
- Channel faults, in the time domain: each horizontal channel is scored by the sum over events of its squared
  zero-lag correlation with the vertical in the P window. A whole-sample shift, or a replacement transducer from a
  coarse grid, is taken when it lifts that score by more than a quarter.
- Orientation: for each event, the principal axis of the horizontal P motion (P - 3 s to P + 5 s, 0.1-1 Hz) gives
  a line; its sense comes from the sign of the projection's correlation with the vertical. The sensor rotation is
  a running circular mean over seven events of (catalogue azimuth + 180 - measured azimuth), and a running vote over
  the same events decides between a straight and a reversed-north wiring.
- Splitting: for trial (phi, dt) the slow component is advanced and the 2 x 2 covariance of the corrected horizontal
  motion (SKS - 12 s to + 22 s, 0.04-0.25 Hz) is formed. The smaller eigenvalue, divided by each event's own minimum,
  is summed over events: first on a 3 degree / 0.1 s grid, then on a 0.5 degree / 0.02 s grid around the best point.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.signal import freqs


def invert(x, dt, T0, h, G):
    f = np.fft.rfftfreq(len(x), dt)
    w0 = 2 * math.pi / T0
    _, H = freqs([G, 0, 0], [1, 2 * h * w0, w0 * w0], worN=2 * math.pi * f)
    level = 0.01 * np.max(np.abs(H))
    Hw = np.where(np.abs(H) < level, level * np.exp(1j * np.angle(H)), H)
    X = np.fft.rfft(x) / Hw
    X[(f < 0.03) | (f > 0.8)] = 0
    return np.fft.irfft(X, len(x))


def band(x, dt, lo, hi):
    f = np.fft.rfftfreq(len(x), dt)
    X = np.fft.rfft(x)
    X[(f < lo) | (f > hi)] = 0
    return np.fft.irfft(X, len(x))


def step_fit(y, t, at, v0, start):
    """T0, h, G of one channel from its calibration step."""
    tau = np.clip(t - at, 0, None)
    on = t >= at

    def shape(T0, h):
        w0 = 2 * math.pi / T0
        q = math.sqrt(1 - h * h)
        return np.where(on, np.exp(-h * w0 * tau) * (np.cos(w0 * q * tau) - h / q * np.sin(w0 * q * tau)), 0.0)

    def cost(p):
        T0, h = math.exp(p[0]), 1 / (1 + math.exp(-p[1]))
        s = shape(T0, h)
        return float(np.sum((y - s * (s @ y) / (s @ s)) ** 2))

    best = min((minimize(cost, [math.log(T), 0.5], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-3})
                for T in (start, 2.0, 6.0, 25.0, 120.0)), key=lambda r: r.fun)
    T0, h = math.exp(best.x[0]), 1 / (1 + math.exp(-best.x[1]))
    s = shape(T0, h)
    return T0, h, float(s @ y / (s @ s)) / v0


def score(x, z):
    """sum over events of the squared zero-lag correlation of x with z."""
    return float(np.sum(np.sum(x * z, 1) ** 2 / np.maximum(np.sum(x * x, 1) * np.sum(z * z, 1), 1e-30)))


def main(folder):
    meta = json.loads((folder / "station.json").read_text())
    with open(folder / "events.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    data = np.load(folder / "windows.npz")
    dt = 1.0 / meta["sample_rate_hz"]
    win = meta["windows"]
    cal = np.load(folder / "calibration.npz")["step"].astype(float)
    tcal = np.arange(cal.shape[1]) * dt
    info = meta["calibration"]
    resp = [step_fit(cal[c], tcal, info["step_at_s"], info["ground_velocity_step_m_s"], ch["natural_period_s"])
            for c, ch in enumerate(meta["channels"])]
    a0, a1 = int((win["p_arrival_s"] - 6) / dt), int((win["p_arrival_s"] + 10) / dt)
    raw_p, raw_s = data["p"].astype(float), data["sks"].astype(float)
    undo = lambda arr, c, r: np.array([invert(tr, dt, *r) for tr in arr[:, c]])
    narrow = lambda arr: np.array([band(tr, dt, 0.15, 0.8) for tr in arr])[:, a0:a1]
    z = narrow(undo(raw_p, 2, resp[2]))
    shift = [0, 0]
    for c in (0, 1):
        now = undo(raw_p, c, resp[c])
        plain = score(narrow(now), z)
        by_lag = max((score(narrow(np.roll(now, -k, axis=1)), z), k) for k in range(-30, 31))
        by_sensor = max((score(narrow(undo(raw_p, c, (T0, h, resp[c][2]))), z), T0, h)
                        for T0 in (2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 9.0, 12.0, 20.0, 40.0)
                        for h in (0.4, 0.55, 0.7, 0.85))
        if max(by_lag[0], by_sensor[0]) > 1.25 * plain:
            if by_lag[0] >= by_sensor[0]:
                shift[c] = by_lag[1]
            else:
                resp[c] = (by_sensor[1], by_sensor[2], resp[c][2])
    P = np.stack([undo(raw_p, c, resp[c]) for c in range(3)], axis=1)
    S = np.stack([undo(raw_s, c, resp[c]) for c in range(3)], axis=1)
    for c in (0, 1):
        if shift[c]:
            P[:, c] = np.roll(P[:, c], -shift[c], axis=1)
            S[:, c] = np.roll(S[:, c], -shift[c], axis=1)
    i0, i1 = int((win["p_arrival_s"] - 3) / dt), int((win["p_arrival_s"] + 5) / dt)
    meas, cat = [], []
    for k, r in enumerate(rows):
        n, e, z = (band(P[k, c], dt, 0.1, 1.0)[i0:i1] for c in range(3))
        cov = np.cov(np.vstack([e, n]))
        vals, vecs = np.linalg.eigh(cov)
        ve, vn = vecs[:, 1]
        if np.dot(ve * e + vn * n, z) < 0:
            ve, vn = -ve, -vn
        meas.append(math.atan2(ve, vn))
        cat.append(math.radians(float(r["back_azimuth_deg"]) + 180))
    meas, cat = np.array(meas), np.array(cat)
    # wiring by a running vote of seven events: which reading of each event lies nearer the catalogue
    plain = np.abs(np.angle(np.exp(1j * (cat - meas))))
    mirror = np.abs(np.angle(np.exp(1j * (cat - (math.pi - meas)))))
    votes = (mirror < plain).astype(float)
    reversed_north = np.array([np.mean(votes[max(0, k - 3):k + 4]) > 0.5 for k in range(len(votes))])
    d = np.angle(np.exp(1j * (cat - np.where(reversed_north, math.pi - meas, meas))))
    rots = [np.angle(np.exp(1j * d[max(0, k - 3):k + 4]).mean()) for k in range(len(d))]

    j0, j1 = int((win["sks_arrival_s"] - 12) / dt), int((win["sks_arrival_s"] + 22) / dt)
    freq = np.fft.rfftfreq(S.shape[2], dt)
    horiz = []
    for k in range(len(rows)):
        n = band(S[k, 0], dt, 0.04, 0.25)
        e = band(S[k, 1], dt, 0.04, 0.25)
        if reversed_north[k]:
            n = -n
        north = n * math.cos(rots[k]) - e * math.sin(rots[k])
        east = n * math.sin(rots[k]) + e * math.cos(rots[k])
        horiz.append((east, north))

    def surface(phis, lags):
        total = np.zeros((len(phis), len(lags)))
        shift = np.exp(2j * math.pi * freq[None, :] * np.asarray(lags)[:, None])
        for east, north in horiz:
            grid = np.zeros_like(total)
            for i, ph in enumerate(np.radians(phis)):
                fast = east * math.sin(ph) + north * math.cos(ph)
                slow = east * math.cos(ph) - north * math.sin(ph)
                adv = np.fft.irfft(np.fft.rfft(slow)[None, :] * shift, len(slow))[:, j0:j1]
                fw = fast[j0:j1][None, :]
                a, b, c = np.sum(fw * fw, 1), np.sum(adv * adv, 1), np.sum(fw * adv, 1)
                grid[i] = 0.5 * (a + b) - np.sqrt(0.25 * (a - b) ** 2 + c * c)
            total += grid / grid.min()
        return total

    phis, lags = np.arange(0, 180, 3.0), np.arange(0, 4.01, 0.1)
    s = surface(phis, lags)
    i, j = np.unravel_index(np.argmin(s), s.shape)
    phis2 = (phis[i] + np.arange(-4, 4.01, 0.5)) % 180
    lags2 = np.clip(lags[j] + np.arange(-0.2, 0.201, 0.02), 0, 4)
    s2 = surface(phis2, lags2)
    i, j = np.unravel_index(np.argmin(s2), s2.shape)
    return {"fast_direction_deg": float(phis2[i] % 180), "split_delay_s": float(lags2[j])}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(main(Path(sys.argv[1])), indent=1) + "\n")
