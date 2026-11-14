"""Independent coupon reduction using source-space spectral projectors."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv


def field(z, kind, f):
    raw = np.fft.rfft(z[kind].astype(float), axis=2)
    lag = z[kind + "_zero_us"]
    sign = z[kind + "_polarity"]
    aligned = np.exp(2j * np.pi * lag[:, :, None, :] * f[None, None, :, None]) * sign[:, :, None, :]
    return np.median(np.fft.irfft(raw * aligned[..., None], n=z[kind].shape[2], axis=2), axis=1)


def projectors(ratio):
    traces = np.trace(ratio, axis1=-2, axis2=-1)
    gaps = np.sqrt(traces**2 - 4 * np.linalg.det(ratio))
    i = int(np.argmax(abs(gaps)))
    a, b = (traces[i] + gaps[i]) / 2, (traces[i] - gaps[i]) / 2
    first = (ratio[i] - b * np.eye(2)) / (a - b)
    return np.array([first, np.eye(2) - first])


def lommel(f, s, z, a):
    """Coaxial equal-piston diffraction factor (Rogers and Van Buren), written for exp(i w t - i k z)."""
    v = 2 * np.pi * f * s * a ** 2 / z
    return 1 - np.exp(1j * v) * (jv(0, v) - 1j * jv(1, v))


def taper(trace, gates, t):
    out = trace.copy()
    for i, (a, b) in enumerate(gates):
        fraction = np.clip(np.minimum(t - a, b - t) / 1.65, 0, 1)
        out[i] *= np.sin(fraction * np.pi / 2)[:, None, None]**2
    return np.fft.rfft(out, axis=1)


def branches(turns, f, dl, low, high):
    """Unwrap each mode's turns over the baselines dl by enumerating the shortest one's branches."""
    d0 = int(np.argmin(dl))
    anchor = f * dl[d0]
    chosen = []
    for row in turns:
        cycles = np.arange(np.ceil(low * anchor - row[d0]), np.floor(high * anchor - row[d0]) + 1)
        choices = (cycles + row[d0]) / anchor
        delta = choices[:, None] * f * dl - row
        residual = delta - np.round(delta)
        chosen.append(row + np.round(delta[np.argmin((residual**2).sum(axis=1))]))
    return chosen


def invert(coupon, coef=None, films=True):
    g = json.loads((coupon / "geometry.json").read_text())
    radius = g.get("transducer_radius_mm", 0.0)
    with np.load(coupon / "waveforms.npz") as z:
        n = z["shot"].shape[2]
        dt = g["sample_interval_us"]
        f = np.fft.rfftfreq(n, dt)
        trace = field(z, "shot", f) - field(z, "lift_off", f)
    t = np.arange(n) * dt
    gates = np.asarray(g["direct_gate_us"])
    spectra = taper(trace, gates, t)
    length = np.asarray(g["thickness_mm"])
    order = np.argsort(length)
    length = length[order]
    spectra = spectra[order]
    dl = np.diff(length)
    bins = np.flatnonzero((f >= g["useful_band_mhz"][0]) & (f <= g["useful_band_mhz"][1]))
    fr = g["reference_frequency_mhz"]
    low, high = g["phase_slowness_bounds_us_per_mm"]
    echo = None
    if films:
        # Each placement carries its own couplant film. A first, film-blind pass only places the
        # triple-transit gate; inside one placement H_direct^-1 H_echo is free of both films.
        rough = invert(coupon, coef=None, films=False)["_slowness"]
        lo, hi = min(rough) * .93, max(rough) * 1.07
        acquired = np.asarray(g["thickness_mm"])
        egates = gates + np.column_stack((2 * lo * acquired, 2 * hi * acquired))
        echo = taper(trace, egates, t)[order]
        path, base = 2 * length, 2 * dl
    else:
        path, base = length, dl
    previous, angles, logs, raw, weights = None, [], [], [], []
    for j in bins:
        if echo is None:
            ratio = np.linalg.solve(spectra[:-1, j], spectra[1:, j])
        else:
            ratio = np.linalg.solve(spectra[:, j], echo[:, j])
        modes = projectors(ratio)
        if previous is not None:
            dist = abs(modes[:, None] - previous[None, :])**2
            if dist[0, 1].sum() + dist[1, 0].sum() < dist[0, 0].sum() + dist[1, 1].sum():
                modes = modes[::-1]
        previous = modes
        eigen = np.einsum("aij,bji->ab", modes, ratio)
        if echo is not None:
            if coef is not None and radius:
                for m in range(2):
                    sm = coef[0, m] + coef[1, m] * (f[j] / fr - 1)
                    eigen[m] = eigen[m] / (lommel(f[j], sm, 3 * length, radius) / lommel(f[j], sm, length, radius))
            pair = abs(eigen[:, 1:] * eigen[:, :-1])
            weights.append(pair**2 / (abs(eigen[:, 1:])**2 + abs(eigen[:, :-1])**2))
            eigen = eigen[:, 1:] / eigen[:, :-1]
        else:
            weights.append(np.ones_like(abs(eigen)))
        if echo is None and coef is not None and radius:      # second pass: remove each mode's beam spreading between neighbours
            for m in range(2):
                sm = coef[0, m] + coef[1, m] * (f[j] / fr - 1)
                eigen[m] = eigen[m] / (lommel(f[j], sm, length[1:], radius) / lommel(f[j], sm, length[:-1], radius))
        raw.append(-np.angle(eigen) / (2 * np.pi))
        angles.append(branches(raw[-1], f[j], base, low, high))
        logs.append(-np.log(abs(eigen)))
    fdl = f[bins, None] * base
    terms = np.column_stack((fdl.flatten(), (fdl * (f[bins, None] / fr - 1)).flatten()))
    w = np.asarray(weights).transpose(1, 0, 2).reshape(2, -1)
    first = np.asarray(angles).transpose(1, 0, 2).reshape(2, -1)
    turns = np.asarray(raw).transpose(1, 0, 2).reshape(2, -1)
    fitted = np.zeros((2, 2))
    for m in range(2):
        root = np.sqrt(w[m])[:, None]
        fitted[:, m] = np.linalg.lstsq(terms * root, first[m] * root[:, 0], rcond=None)[0]
        for _ in range(3):      # re-unwrap every bin against the band-wide line, which a few noisy bins cannot move
            model = terms @ fitted[:, m]
            fitted[:, m] = np.linalg.lstsq(terms * root, (turns[m] + np.round(model - turns[m])) * root[:, 0],
                                           rcond=None)[0]
    if coef is None and radius and films:
        return invert(coupon, fitted, films)
    s0 = fitted[0]
    lw = np.asarray(weights)
    a0 = (lw * np.asarray(logs) * fdl[:, None, :]).sum(axis=(0, 2)) * fr / (lw * fdl[:, None, :]**2).sum(axis=(0, 2))
    result = {"_slowness": [float(v) for v in s0]}
    for label, i in zip(("fast", "slow"), np.argsort(s0)):
        k = complex(2 * np.pi * fr * s0[i], -a0[i])
        modulus = g["density_kg_m3"] * (2 * np.pi * fr / k)**2
        result[f"{label}_storage_mpa"] = modulus.real
        result[f"{label}_loss_mpa"] = modulus.imag
    return result


if __name__ == "__main__":
    answer = {k: v for k, v in invert(Path(sys.argv[1])).items() if not k.startswith("_")}
    Path(sys.argv[2]).write_text(json.dumps(answer, allow_nan=False), encoding="utf8")
