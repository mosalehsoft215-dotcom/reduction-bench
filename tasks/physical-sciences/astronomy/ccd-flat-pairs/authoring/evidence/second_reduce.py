"""Independent second reduction, written apart from solution/bf_reduce.py, for the allowances.

    python3 second_reduce.py <sensor folder> <result file>

Where it differs from the reference:
- row offsets: clipped mean of overscan reads 6-15; the column pattern from the median zero frame;
- linearity: log(s / m) = b_range + log(1 - d2 x^2 - d3 x^3), fitted by scipy least squares, then the
  response is inverted pixel by pixel with Newton steps;
- deferred charge: the trail is a geometric fit to the first four overscan reads of flats at 15-45 % of the top, checked
  against the deficit of the first imaging columns; it is undone with the inverse
  IIR filter in scipy.signal.lfilter;
- the same check down the columns: a geometric fit to the deficit of the first four rows gives an inverse filter
  along the rows, which is kept only if the series model then fits the covariances better (a dim edge of the die
  dims the first rows too, without holding charge);
- pair statistics: 64 x 64 tiles, each with its own plane removed, covariances averaged over tiles;
- fit range: the run stops before the first level whose flat has a negative pixel skewness below -0.08;
- model: the covariance series C = mu [delta + a mu + (2/3) a*a mu^2 + (1/3) a*a*a mu^3 + ...] (to sixth
  order, convolutions in real space), the kernel summing to zero, with white noise at lag 0 and a row-offset
  constant on the lags within a row, fitted to six lags.
"""
import csv
import json
import sys
from math import factorial
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import convolve2d, lfilter

TILE = 64
LAGS6 = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
HALF = 10


def clipped_row_offsets(stack, ncol):
    over = stack[..., ncol + 6:]
    med = np.median(over, axis=-1, keepdims=True)
    spread = 1.4826 * np.median(np.abs(over - med), axis=-1, keepdims=True) + 1e-9
    keep = np.abs(over - med) < 4 * spread
    return stack - (np.where(keep, over, 0).sum(-1, keepdims=True) / keep.sum(-1, keepdims=True))


def linear_response(img, log, ncol):
    sig = img[..., :ncol].mean(axis=(-1, -2))
    s = np.array([sig[int(r["pair"]), int(r["frame"])] for r in log])
    m = np.array([float(r["monitor_nC"]) for r in log])
    names = sorted({r["monitor_range"] for r in log})
    which = np.array([names.index(r["monitor_range"]) for r in log])
    x = s / 6e4

    def resid(p):
        b, d2, d3 = p[: len(names)], p[-2], p[-1]
        return np.log(s / m) - b[which] - np.log(1 - d2 * x * x - d3 * x ** 3)

    p0 = [np.median(np.log(s / m)[which == i]) for i in range(len(names))] + [0.0, 0.0]
    fit = least_squares(resid, p0).x
    d2, d3 = fit[-2], fit[-1]
    if abs(d2) < 1e-5 and abs(d3) < 1e-5:
        return img
    y = img / 6e4                       # measured; find u with u (1 - d2 u^2 - d3 u^3) = y
    u = y.copy()
    for _ in range(6):
        f = u * (1 - d2 * u * u - d3 * u ** 3) - y
        u = u - f / (1 - 3 * d2 * u * u - 4 * d3 * u ** 3)
    return u * 6e4


def geometric(values):
    """Fits v_k = A r^k by a log-linear regression over the leading run of positive values (from k = 0, at
    least two); returns (A, r), or (0, 0) when the profile does not start positive."""
    run = 0
    while run < len(values) and values[run] > 0:
        run += 1
    if run < 2:
        return 0.0, 0.0
    v = values[:run]
    slope, icpt = np.polyfit(np.arange(run), np.log(v), 1, w=v / v[0])    # small, noisy values count less
    return float(np.exp(icpt)), float(np.exp(slope))


def lit(img, ncol):
    sig = img[..., :ncol].mean(axis=(-1, -2)).ravel()
    order = np.argsort(sig)
    band = order[(sig[order] > 0.15 * sig.max()) & (sig[order] < 0.45 * sig.max())]
    return img.reshape(-1, *img.shape[-2:])[band]


def deferred_charge(img, ncol):
    bright = lit(img, ncol)
    level = bright[:, :, 64:ncol - 8].mean()
    trail = bright[:, :, ncol:ncol + 4].mean(axis=(0, 1)) / bright[:, :, ncol - 1].mean()
    lead = 1 - bright[:, :, :4].mean(axis=(0, 1)) / level
    trail_a, trail_r = geometric(trail)
    if trail_a < 3e-4 or not (0.7 < lead[0] / trail[0] < 1.4 and lead[1] > 0):
        return img[..., :ncol]
    eps, r = trail_a, min(max(trail_r, 0.0), 0.9)
    return lfilter([1.0, -r], [1.0 - eps, -(r - eps)], img, axis=-1)[..., :ncol]


def row_deferred_charge(img):
    """Returns the image with charge held in the parallel transfer undone, or None if the rows show no deficit."""
    n = img.shape[-2]
    rows = lit(img, img.shape[-1]).mean(axis=(0, 2))
    lead_a, lead_r = geometric(1 - rows[:4] / rows[64:n - 8].mean())
    if lead_a < 2.5e-3 or not (0 < lead_r < 0.9):
        return None
    return lfilter([1.0, -lead_r], [1.0 - lead_a, -(lead_r - lead_a)], img, axis=-2)


def tile_covariances(d):
    n, m = d.shape
    yy, xx = np.mgrid[0:TILE, 0:TILE] / TILE - 0.5
    basis = np.stack([np.ones_like(xx), xx, yy], -1).reshape(-1, 3)
    acc = np.zeros(len(LAGS6))
    count = 0
    for i in range(0, n, TILE):
        for j in range(0, m, TILE):
            t = d[i:i + TILE, j:j + TILE]
            coef, *_ = np.linalg.lstsq(basis, t.ravel(), rcond=None)
            t = t - (basis @ coef).reshape(TILE, TILE)
            acc += [(t[: TILE - a, : TILE - b] * t[a:, b:]).mean() for a, b in LAGS6]
            count += 1
    return acc / count / 2


def skewness(f):
    z = f - f.mean()
    return float((z ** 3).mean() / (z ** 2).mean() ** 1.5)


def kernel_image(q):
    """q = (a01, a10, a11, a02, a20) in 1e-6; a00 closes the sum."""
    a = np.zeros((2 * HALF + 1, 2 * HALF + 1))
    c = HALF
    for (dy, dx), v in zip([(1, 0), (0, 1), (1, 1), (2, 0), (0, 2)], q):
        for sy in ((1,) if dy == 0 else (1, -1)):
            for sx in ((1,) if dx == 0 else (1, -1)):
                a[c + sy * dy, c + sx * dx] += v * 1e-6
    a[c, c] = -a.sum()
    return a


def series_model(p, mu):
    g, noise, q = p[0], p[1], p[2:7]
    a = kernel_image(q)
    powers = [None, a]
    for _ in range(5):
        powers.append(convolve2d(powers[-1], a, mode="same"))
    c = HALF
    out = []
    for level in mu:
        e = level * g
        cv = []
        for dy, dx in LAGS6:
            v = (1.0 if (dy, dx) == (0, 0) else 0.0)
            for n in range(1, 6):
                v += (2 ** n) / factorial(n + 1) * powers[n][c + dy, c + dx] * e ** n
            cv.append(e * v / g ** 2)
        out.append(cv)
    out = np.array(out)
    out[:, 0] += noise
    out[:, [0, 2, 5]] += p[-1]        # row-offset residual: lags (0,0), (0,1), (0,2)
    return out


def series_fit(img):
    mu, cv, sk = [], [], []
    for f1, f2 in img:
        d = f1 - f2 * (f1.mean() / f2.mean())
        mu.append(0.5 * (f1.mean() + f2.mean()))
        cv.append(tile_covariances(d))
        sk.append(skewness(f1))
    order = np.argsort(mu)
    mu, cv, sk = np.array(mu)[order], np.array(cv)[order], np.array(sk)[order]
    bad = np.nonzero(sk < -0.08)[0]
    top = bad[0] if len(bad) else len(mu)
    mu, cv = mu[:top], cv[:top]
    npix = img.shape[-1] * img.shape[-2]
    sd = np.maximum(cv[:, :1], 1.0) / np.sqrt(npix) * np.array([1.414, 1, 1, 1, 1, 1])
    res = lambda p: ((series_model(p, mu) - cv) / sd).ravel()
    fit = least_squares(res, [1.5, 10.0, 0.6, 0.4, 0.1, 0.03, 0.03, 0.0], x_scale=[0.1, 5, 0.1, 0.1, 0.05, 0.02, 0.02, 0.5])
    return fit.x, fit.cost, top


def reduce(folder):
    card = json.loads((folder / "sensor_card.json").read_text())
    ncol = card["imaging_columns"]
    frames = np.load(folder / "flat_pairs.npz")["frames"].astype(float)
    zeros = np.load(folder / "zero_frames.npz")["frames"].astype(float)
    with open(folder / "exposure_log.csv") as f:
        log = list(csv.DictReader(f))
    pattern = np.median(clipped_row_offsets(zeros, ncol), axis=0)
    img = clipped_row_offsets(frames, ncol) - pattern
    img = linear_response(img, log, ncol)
    img = deferred_charge(img, ncol)
    best = series_fit(img)
    undone = row_deferred_charge(img)
    if undone is not None:
        other = series_fit(undone)       # dim first rows: held charge or a dim edge; the fit decides
        if other[1] < best[1]:
            best = other
    p, cost, top = best
    a00 = kernel_image(p[2:7])[HALF, HALF]
    return {"conversion_gain": float(p[0]), "area_coeff_a00": float(a00), "levels_used": int(top)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(Path(sys.argv[1])), indent=1) + "\n")
