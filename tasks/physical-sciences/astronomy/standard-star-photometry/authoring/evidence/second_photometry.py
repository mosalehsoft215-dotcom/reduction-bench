"""Second, separately written reduction of a B/V epoch of CCD stamps (same call as the reference).

Differences from solution/photometry_reduce.py:
- every Moffat fit is a plain nonlinear least-squares fit of all parameters together (centre, FWHM, flux, a constant
  sky and its two slopes), started from image moments, instead of solving flux and sky linearly;
- the companion offset is shared by all target stamps: it starts from the median over the stamps of the residual
  maximum within 8 pixels after one-star fits, and every target stamp is then fitted together in one nonlinear
  least-squares problem with a sparse Jacobian (per stamp: sky plane, FWHM, centre and both fluxes; shared: the
  offset), tried from four starting offsets (the residual offset, 1.3 times it, and both mirrored); the second star is kept
  when chi-square drops by more than 25 per stamp and the fainter-to-brighter flux ratio is steady within each filter
  (median absolute deviation under a quarter of the median, median above 0.03), and in every stamp the brighter of the two is the target;
- the shutter delta is the 20 % trimmed mean of the short/long pairs' individual delays, not a count-weighted median;
- the transformation is fitted by iteratively reweighted least squares (Tukey biweight) instead of hard clipping;
- in the target stamps, pixels the fitted model (two stars, or one when no companion is kept) cannot explain (cosmic
  rays, trails) are tamed by a Cauchy loss of scale 3 sigma, and for a single star the pixels still beyond 5 sigma
  (with their eight neighbours, the core excepted) are then left out of a final plain fit;
- the target is solved visit by visit (its B and V frames taken together): a visit whose V or B - V lies more than
  0.05 mag from the median of the visits is dropped (the two nearest the median are kept when fewer remain), and V
  and B - V come from the frames of the remaining visits.
"""
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

BETA, NL, SIDE, SS = 3.2, 1.2e-6, 25, 5
_o = (np.arange(SS) + 0.5) / SS - 0.5
_x = (np.arange(SIDE)[:, None] + _o).ravel()
GX, GY = np.meshgrid(_x, _x)
PIX = np.arange(SIDE) - SIDE // 2


def moffat(x0, y0, fwhm):
    al = fwhm / (2 * np.sqrt(2 ** (1 / BETA) - 1))
    p = (BETA - 1) / (np.pi * al * al) * (1 + ((GX - x0) ** 2 + (GY - y0) ** 2) / al ** 2) ** (-BETA)
    return p.reshape(SIDE, SS, SIDE, SS).mean(axis=(1, 3))


def model(q, nstar):
    sky = q[0] + q[1] * PIX[None, :] + q[2] * PIX[:, None]
    fw = q[3]
    out = sky.copy()
    for k in range(nstar):
        x0, y0, fl = q[4 + 3 * k:7 + 3 * k]
        out = out + fl * moffat(x0, y0, fw)
    return out


def grow(bad):
    out = bad.copy()
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            out |= np.roll(np.roll(bad, dy, 0), dx, 1)
    return out


def fit(img, nstar, start, clean=False):
    """Least squares of the whole model.  With clean=True: a Cauchy-loss fit (scale 3 sigma) first, then the pixels
    still more than 5 sigma off it (with their eight neighbours, the star's 1.5 px core excepted) are left out of a
    final plain least-squares fit."""
    sig = np.sqrt(np.clip(img, 0, None) + 36.0)
    if not clean:
        sol = least_squares(lambda q: ((model(q, nstar) - img) / sig).ravel(), start, x_scale="jac")
        return sol.x, np.sum(sol.fun ** 2), np.ones(img.shape, bool)
    sol = least_squares(lambda q: ((model(q, nstar) - img) / sig).ravel(), start, x_scale="jac", loss="cauchy",
                        f_scale=3.0)
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    core = np.hypot(xx - sol.x[4], yy - sol.x[5]) < 1.5
    good = ~(grow(np.abs(model(sol.x, nstar) - img) / sig > 5.0) & ~core)
    sol = least_squares(lambda q: (good * (model(q, nstar) - img) / sig).ravel(), sol.x, x_scale="jac")
    return sol.x, np.sum(sol.fun ** 2), good


def one_star(img, clean=False):
    sky0 = np.median(np.r_[img[0], img[-1], img[:, 0], img[:, -1]])
    s = np.clip(img - sky0, 0, None)
    tot = s.sum()
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    x0, y0 = (s * xx).sum() / tot, (s * yy).sum() / tot
    return fit(img, 1, [sky0, 0.0, 0.0, 5.0, x0, y0, tot], clean)


def target_fluxes(imgs, filters):
    """All target stamps together: one companion offset shared by every stamp, two stars of one FWHM per stamp."""
    singles = [one_star(im) for im in imgs]
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    offs = []
    for (q, _, gd), im in zip(singles, imgs):
        resid = np.where(gd, im - model(q, 1), 0.0)
        zone = (np.hypot(xx - q[4], yy - q[5]) > 2.0) & (np.hypot(xx - q[4], yy - q[5]) < 8.0)
        py, px = np.unravel_index(np.argmax(np.where(zone, resid, -np.inf)), resid.shape)
        offs.append((px - q[4], py - q[5]))
    dx0, dy0 = np.median(np.array(offs), axis=0)
    n, k = len(imgs), 8
    sig = [np.sqrt(np.clip(im, 0, None) + 36.0) for im in imgs]
    goods = [np.ones((SIDE, SIDE), bool) for _ in imgs]

    def res(p):
        dx, dy = p[-2:]
        out = []
        for i, im in enumerate(imgs):
            sky, gx, gy, fw, x0, y0, f1, f2 = p[k * i:k * i + k]
            mdl = (sky + gx * PIX[None, :] + gy * PIX[:, None] + f1 * moffat(x0, y0, fw)
                   + f2 * moffat(x0 + dx, y0 + dy, fw))
            out.append((goods[i] * (mdl - im) / sig[i]).ravel())
        return np.concatenate(out)
    sparse = np.zeros((n * SIDE * SIDE, k * n + 2), int)
    for i in range(n):
        sparse[i * SIDE * SIDE:(i + 1) * SIDE * SIDE, k * i:k * i + k] = 1
    sparse[:, -2:] = 1
    best = None
    for scale in (1.0, 1.3, -1.0, -1.3):                      # both sides: the blend is nearly symmetric
        ox, oy = scale * dx0, scale * dy0
        p0 = []
        for q, _, _ in singles:
            p0 += [q[0], q[1], q[2], q[3], q[4] - 0.2 * ox, q[5] - 0.2 * oy, 0.8 * q[6], 0.25 * q[6]]
        sol = least_squares(res, np.r_[p0, ox, oy], jac_sparsity=sparse, x_scale="jac")
        if best is None or sol.cost < best.cost:
            best = sol
    first, other = best.x[6::k][:n], best.x[7::k][:n]
    ratio = np.minimum(first, other) / np.maximum(first, other)
    def steady_in(sel):
        r = ratio[sel]
        med = np.median(r)
        return med > 0.03 and np.median(np.abs(r - med)) < 0.25 * med
    steady = all(steady_in(filters == f) for f in set(filters))
    # a real companion has a steady flux ratio to the target within each filter; hits and trails do not
    if steady and sum(c for _, c, _ in singles) - 2 * best.cost > 25.0 * n:
        # pixels the two-star model cannot explain (cosmic rays, trails) lose their pull in one more joint fit with a
        # Cauchy loss (scale 3 sigma)
        best = least_squares(res, best.x, jac_sparsity=sparse, x_scale="jac", loss="cauchy", f_scale=3.0)
        first, other = best.x[6::k][:n], best.x[7::k][:n]
        return list(np.maximum(first, other))                 # in every stamp the target is the brighter star
    return [one_star(im, clean=True)[0][6] for im in imgs]


def reduce_night(folder):
    folder = Path(folder)
    head = json.loads((folder / "epoch.json").read_text())
    with open(folder / "frames.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    raw = np.load(folder / "stamps.npy").astype(float) - head["bias_pedestal_adu"]
    raw = raw * (1 + NL * np.clip(raw, 0, None))
    cat = {s["star"]: (s["V"], s["B_V"]) for s in json.loads((folder / "standards.json").read_text())}
    tgt = head["target"]
    flux = np.array([0.0 if r["star"] == tgt else one_star(raw[i])[0][6] for i, r in enumerate(rows)])
    tidx = [i for i, r in enumerate(rows) if r["star"] == tgt]
    flux[tidx] = target_fluxes([raw[i] for i in tidx], np.array([rows[i]["filter"] for i in tidx]))
    t_cmd = np.array([float(r["exposure_s"]) for r in rows])
    ut = np.array([float(r["ut_hours"]) for r in rows])
    xm = np.array([float(r["airmass"]) for r in rows])
    # shutter: for each short/long pair, F_l (t_s + d) = F_s (t_l + d) -> d (F_l - F_s) = F_s t_l - F_l t_s
    lhs, rhs = [], []
    for i in range(len(rows) - 2):
        for j in (i + 1, i + 2):
            if rows[i]["star"] == rows[j]["star"] and rows[i]["filter"] == rows[j]["filter"] and \
                    rows[i]["star"] in cat and abs(t_cmd[j] / t_cmd[i] - 10) < 0.5 and ut[j] - ut[i] < 0.05:
                lhs.append(flux[j] - flux[i])
                rhs.append(flux[i] * t_cmd[j] - flux[j] * t_cmd[i])
    each = np.sort(np.array(rhs) / np.array(lhs))       # one delay per pair
    cut = len(each) // 5
    d = float(np.mean(each[cut:len(each) - cut]))         # 20 % trimmed mean of the pairs' delays
    mag = -2.5 * np.log10(flux / (t_cmd + d))
    t0 = ut.min()
    coefs = {}
    for f in ("B", "V"):
        idx = [i for i, r in enumerate(rows) if r["filter"] == f and r["star"] in cat]
        a = np.array([[1.0, xm[i], (ut[i] - t0) * xm[i], cat[rows[i]["star"]][1] * xm[i], cat[rows[i]["star"]][1]]
                      for i in idx])
        y = np.array([mag[i] - (cat[rows[i]["star"]][0] + (cat[rows[i]["star"]][1] if f == "B" else 0.0)) for i in idx])
        w = np.ones(len(idx))
        for _ in range(10):
            c = np.linalg.lstsq(a * w[:, None], y * w, rcond=None)[0]
            res = y - a @ c
            s = 1.4826 * np.median(np.abs(res)) + 1e-4
            u = res / (4.685 * s)
            w = np.sqrt(np.clip(1 - u * u, 0, None) ** 2)
        coefs[f] = c
    rows_t = sorted((i for i, r in enumerate(rows) if r["star"] == tgt), key=lambda i: ut[i])
    a, b = [], []
    for i in rows_t:
        z, k1, kr, k2, eps = coefs[rows[i]["filter"]]
        a.append([1.0, (1.0 if rows[i]["filter"] == "B" else 0.0) + k2 * xm[i] + eps])
        b.append(mag[i] - z - k1 * xm[i] - kr * (ut[i] - t0) * xm[i])
    a, b = np.array(a), np.array(b)
    # visits: consecutive B and V frames of the target (within 0.1 h); each gives its own V and B - V
    visits, cur = [], [0]
    for j in range(1, len(rows_t)):
        if ut[rows_t[j]] - ut[rows_t[cur[0]]] < 0.1:
            cur.append(j)
        else:
            visits.append(cur)
            cur = [j]
    visits.append(cur)
    per = np.array([np.linalg.lstsq(a[g], b[g], rcond=None)[0] for g in visits])
    off = np.abs(per - np.median(per, axis=0)).max(axis=1)
    good = [i for i, o in enumerate(off) if o <= 0.05]
    if len(good) < 2:                                         # no consensus: the two visits nearest the median
        good = list(np.argsort(off)[:2])
    keep = [j for i in good for j in visits[i]]
    v, c = np.linalg.lstsq(a[keep], b[keep], rcond=None)[0]
    return {"target_v_mag": float(v), "target_b_minus_v": float(c)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_night(sys.argv[1])))
