"""Reference reduction of one flat-pair run: conversion gain and the area-kernel coefficient a00.

    python3 bf_reduce.py <sensor folder> <result file>

Steps
1. Offsets: each row's offset is the median of the last eight overscan reads; the zero frames give the residual
   fixed pattern of every column, which is removed as well.
2. Linearity: the mean signal s of each flat against its monitor charge m. The model is
   s (1 + c2 x^2 + c3 x^3) = k_r m with x = s / 60000 and one scale k_r per electrometer range, so a range that is
   calibrated differently cannot pose as a sensor response; c2 and c3 then undo the response pixel by pixel.
3. Deferred charge: in flats at 20-50 % of the top level the first overscan reads show the trail of each row. Charge that
   trails behind must be missing at the start of the row, so the correction is applied only when the deficit of
   the first imaging columns matches the trail; the held fraction and the release ratio are then undone row by row.
   The same is asked of the columns, where there is no overscan to show a trail: if the first rows of the
   mid-level flats fall short in a geometric profile, that is either charge held back in the parallel transfer
   (which also correlates neighbouring rows) or an edge of the die that collects less light (which does not).
   The reduction is done both ways and the one whose covariances the model fits better is kept.
4. Pair statistics: the second flat is scaled to the first, the difference has a quadratic surface removed
   (illumination that moved between the two flats), and its covariances at eight lags are taken.
5. Fit range: the run is used up to the first level whose flat has lost the upper tail of its pixel
   distribution (pixels at full well spill into their neighbours).
6. Model: for the area kernel a(k) in Fourier space, the covariance spectrum of the difference grows as
   P(k) = (exp(2 a(k) mu g) - 1) / (2 a(k)); its inverse transform at each lag, divided by g^2, plus white read
   noise at lag 0 and one constant shared by all lags within a row (the residual of the row offsets), is fitted
   to all lags and levels at once.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

LAGS = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 1), (1, 2)]
NEIGH = [(1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
ROW_LAGS = [i for i, (dy, dx) in enumerate(LAGS) if dy == 0]
GRID = 48


def load(folder):
    card = json.loads((folder / "sensor_card.json").read_text())
    frames = np.load(folder / "flat_pairs.npz")["frames"].astype(float)
    zeros = np.load(folder / "zero_frames.npz")["frames"].astype(float)
    with open(folder / "exposure_log.csv") as f:
        log = list(csv.DictReader(f))
    return card, frames, zeros, log


def remove_offsets(stack, ncol):
    return stack - np.median(stack[..., ncol + 8:], axis=-1, keepdims=True)


def linearise(img, log, ncol):
    sig = img[..., :ncol].mean(axis=(-1, -2))
    names = sorted({r["monitor_range"] for r in log})
    rows, s_all = [], []
    for r in log:
        s = sig[int(r["pair"]), int(r["frame"])]
        m = float(r["monitor_nC"])
        x = s / 6e4
        rows.append([m * (r["monitor_range"] == n) for n in names] + [-s * x * x, -s * x ** 3])
        s_all.append(s)
    a = np.array(rows) / np.array(s_all)[:, None]
    coef, *_ = np.linalg.lstsq(a, np.ones(len(s_all)), rcond=None)
    c2, c3 = coef[-2:]
    x = img / 6e4
    return img * (1 + c2 * x * x + c3 * x ** 3)


def mid_flats(img, ncol):
    """Flats between 20 % and 50 % of the brightest level: well lit, and below any full-well spill."""
    sig = img[..., :ncol].mean(axis=(-1, -2)).ravel()
    pick = (sig > 0.2 * sig.max()) & (sig < 0.5 * sig.max())
    return img.reshape(-1, *img.shape[-2:])[pick]


def undo_trailing(img, ncol):
    bright = mid_flats(img, ncol)
    inner = bright[:, :, 40:ncol - 4].mean()
    trail = bright[:, :, ncol:ncol + 3].mean(axis=(0, 1)) / bright[:, :, ncol - 1].mean()
    deficit = 1 - bright[:, :, :3].mean(axis=(0, 1)) / inner
    eps = trail[0]
    if not (eps > 3e-4 and abs(deficit[0] / eps - 1) < 0.3):
        return img[..., :ncol]
    ratio = min(max(trail[1] / trail[0], 0.0), 0.9)
    out = np.zeros(img.shape[:-1] + (ncol,))
    held = np.zeros(img.shape[:-1])
    for j in range(ncol):
        q = (img[..., j] - (1 - ratio) * held) / (1 - eps)
        out[..., j] = q
        held = held * ratio + eps * q
    return out


def row_start_profile(img):
    n = img.shape[-2]
    rows = mid_flats(img, img.shape[-1]).mean(axis=(0, 2))
    short = 1 - rows[:4] / rows[40:n - 8].mean()
    if short[0] > 2.5e-3 and 0 < short[1] / short[0] < 0.9 and short[2] < short[1]:
        return short[0], short[1] / short[0]
    return None


def undo_row_trailing(img, eps, ratio):
    n = img.shape[-2]
    out = np.zeros_like(img)
    held = np.zeros(img.shape[:-2] + img.shape[-1:])
    for i in range(n):
        q = (img[..., i, :] - (1 - ratio) * held) / (1 - eps)
        out[..., i, :] = q
        held = held * ratio + eps * q
    return out


def without_surface(d):
    n, m = d.shape
    yy, xx = np.mgrid[0:n, 0:m]
    yy = yy / n - 0.5
    xx = xx / m - 0.5
    basis = np.stack([np.ones_like(xx), xx, yy, xx * xx, yy * yy, xx * yy], -1).reshape(-1, 6)
    coef, *_ = np.linalg.lstsq(basis, d.ravel(), rcond=None)
    return d - (basis @ coef).reshape(n, m)


def pair_statistics(img):
    mu, cov, tail = [], [], []
    n, m = img.shape[-2:]
    for f1, f2 in img:
        m1, m2 = f1.mean(), f2.mean()
        d = without_surface(f1 - f2 * (m1 / m2))
        hi, mid, lo = np.quantile(f1, [0.995, 0.5, 0.005])
        tail.append((hi - mid) / (mid - lo))
        mu.append(0.5 * (m1 + m2))
        cov.append([(d[: n - dy, : m - dx] * d[dy:, dx:]).mean() / 2 for dy, dx in LAGS])
    return np.array(mu), np.array(cov), np.array(tail)


def model(p, mu):
    g, a00, *rest = p[:7]
    a = np.zeros((GRID, GRID))
    a[0, 0] = a00 * 1e-6
    for (dy, dx), v in zip(NEIGH, rest):
        for sy in ((1,) if dy == 0 else (1, -1)):
            for sx in ((1,) if dx == 0 else (1, -1)):
                a[(sy * dy) % GRID, (sx * dx) % GRID] += v * 1e-6
    ah = np.real(np.fft.fft2(a))
    ky = 2 * np.pi * np.fft.fftfreq(GRID)[:, None]
    kx = 2 * np.pi * np.fft.fftfreq(GRID)[None, :]
    waves = np.stack([np.cos(ky * dy + kx * dx) for dy, dx in LAGS])
    flat = np.abs(ah) < 1e-15
    safe = np.where(flat, 1.0, ah)
    e = (np.asarray(mu) * g)[:, None, None]
    spec = np.where(flat[None], e, np.expm1(2 * ah[None] * e) / (2 * safe[None]))
    out = np.einsum("mij,lij->ml", spec, waves) / GRID ** 2 / g ** 2
    out[:, 0] += p[7]
    out[:, ROW_LAGS] += p[8]          # what is left of each row's offset is shared along the row
    return out


def fit(mu, cov, npix):
    sd = np.maximum(cov[:, :1], 1.0) / np.sqrt(npix) * np.array([1.414] + [1.0] * (len(LAGS) - 1))
    res = lambda p: ((model(p, mu) - cov) / sd).ravel()
    start = [1.5, -2.5, 0.6, 0.4, 0.1, 0.03, 0.03, 10.0, 0.0]
    best = least_squares(res, start, x_scale=[0.1, 0.5, 0.1, 0.1, 0.05, 0.02, 0.02, 5.0, 0.5])
    return best.x, 2 * best.cost


def measure(img):
    mu, cov, tail = pair_statistics(img)
    order = np.argsort(mu)
    mu, cov, tail = mu[order], cov[order], tail[order]
    lost = np.nonzero(tail < 0.9)[0]
    top = lost[0] if len(lost) else len(mu)
    p, chi2 = fit(mu[:top], cov[:top], img.shape[-1] * img.shape[-2])
    return p, chi2, top, len(mu)


def reduce(folder):
    card, frames, zeros, log = load(folder)
    ncol = card["imaging_columns"]
    fixed = remove_offsets(zeros, ncol).mean(axis=0)
    img = remove_offsets(frames, ncol) - fixed
    img = linearise(img, log, ncol)
    img = undo_trailing(img, ncol)
    p, chi2, top, taken = measure(img)
    held = row_start_profile(img)
    rows_undone = False
    if held:
        # Dim first rows are either charge held back in the parallel transfer (which also correlates
        # neighbouring rows) or an edge that collects less light (which does not). The covariances decide.
        q, chi2_q, top_q, _ = measure(undo_row_trailing(img, *held))
        if chi2_q < chi2:
            p, chi2, top, rows_undone = q, chi2_q, top_q, True
    return {"conversion_gain": float(p[0]), "area_coeff_a00": float(p[1] * 1e-6),
            "levels_used": int(top), "levels_taken": int(taken), "row_trailing_undone": rows_undone}


if __name__ == "__main__":
    result = reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result, indent=1) + "\n")
