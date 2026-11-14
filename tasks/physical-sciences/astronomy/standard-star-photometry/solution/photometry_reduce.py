"""Reduces one B/V epoch of CCD stamps to the target's V magnitude and B - V colour.

    python3 photometry_reduce.py <epoch directory> <result.json>

0. Pixels: bias pedestal removed, then the CCD response law t = m (1 + 1.2e-6 m) applied pixel by pixel.
1. PSF photometry: every stamp is fitted with the Moffat PSF of the notes (beta 3.2) integrated over the pixels,
   position and FWHM free, flux and a sky plane solved linearly; pixels more than five noise units off the fitted
   model (cosmic-ray hits, satellite trails), widened by two pixels but never within 1.5 px of a star's centre,
   are left out and the fit repeated (twice).
   The target's stamps are fitted together for a possible companion: the median residual of single-star fits,
   aligned on the target, gives its offset; all target stamps are refitted with two stars at a common offset
   (fluxes free per stamp), started from that offset and from its mirror image, the second star is kept when it lowers chi-square by more than 25 per stamp, and in every stamp
   the brighter of the two is the target, whose flux alone is used; the target stamps are then cleaned of
   unexplained pixels the same way and refitted.
2. Shutter: every standard is taken twice in a row, short and ten times longer; delta from
   F_short / (t_short + delta) = F_long / (t_long + delta), count-weighted median over the epoch.
3. Transformation, per filter, by linear least squares over the standard frames:
   m - M = zeta + k1 X + r (t - t_first) X + k2 c X + eps c, m = -2.5 log10(F / (t + delta)).
4. Target: its B and V frames give a linear system in (V, c), solved by least squares; while the worst frame lies
   more than 0.05 mag off the solution (a frame through the wrong filter, or under cloud) it is dropped and the
   system solved again, at least two frames being kept in each filter.
"""
import csv
import json
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.ndimage import binary_dilation  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

BETA_PSF, NONLIN, N = 3.2, 1.2e-6, 25
SUB = 5
_OFF = (np.arange(SUB) + 0.5) / SUB - 0.5
_XS = (np.arange(N)[:, None] + _OFF[None, :]).ravel()
_XX, _YY = np.meshgrid(_XS, _XS)
_G = np.arange(N) - N // 2


def psf(x0, y0, fwhm):
    alpha = fwhm / (2 * np.sqrt(2 ** (1 / BETA_PSF) - 1))
    r2 = (_XX - x0) ** 2 + (_YY - y0) ** 2
    prof = (BETA_PSF - 1) / (np.pi * alpha ** 2) * (1 + r2 / alpha ** 2) ** (-BETA_PSF)
    return prof.reshape(N, SUB, N, SUB).mean(axis=(1, 3))


def linear_solve(cols, img, w):
    a = np.column_stack([c.ravel() for c in cols] + [np.ones(N * N), np.repeat(_G, N) * 0 + np.tile(_G, N),
                                                     np.repeat(_G, N)])
    coef, *_ = np.linalg.lstsq(a * w[:, None], img.ravel() * w, rcond=None)
    return coef, a @ coef


_YI, _XI = np.mgrid[0:N, 0:N]


def hit_mask(img, model, centres):
    """Pixels the model cannot explain (cosmic-ray hits, satellite trails): more than five noise units off the model,
    widened by two pixels all round, except within 1.5 px of a fitted star centre (a star's core is always kept).
    Returns the weight multiplier (0 on those pixels)."""
    z = np.abs(img.ravel() - model) / np.sqrt(np.clip(model, 0, None) + 36.0)
    bad = binary_dilation((z > 5.0).reshape(N, N), np.ones((5, 5), bool))
    for cx, cy in centres:
        bad &= np.hypot(_XI - cx, _YI - cy) >= 1.5
    return (~bad).ravel().astype(float)


def fit_single(img, model_weights=True, keep=None):
    keep = np.ones(N * N) if keep is None else keep
    w = keep / np.sqrt(np.clip(img, 0, None).ravel() + 36.0)
    yc, xc = np.unravel_index(np.argmax(img), img.shape)

    def res(q):
        _, model = linear_solve([psf(*q)], img, w)
        return (model - img.ravel()) * w
    sol = least_squares(res, [xc, yc, 5.0], bounds=([0, 0, 1.5], [N - 1, N - 1, 12.0]), x_scale=[0.3, 0.3, 0.3])
    coef, model = linear_solve([psf(*sol.x)], img, w)
    if not model_weights:
        return sol.x, coef[0], (img.ravel() - model).reshape(N, N), w
    w = keep / np.sqrt(np.clip(model, 0, None) + 36.0)       # weights from the model, not the noisy data
    sol = least_squares(res, sol.x, bounds=([0, 0, 1.5], [N - 1, N - 1, 12.0]), x_scale=[0.3, 0.3, 0.3])
    coef, model = linear_solve([psf(*sol.x)], img, w)
    return sol.x, coef[0], (img.ravel() - model).reshape(N, N), w


def fit_clean(img, hits=True):
    """A single-star fit, refitted twice without the pixels it cannot explain."""
    out = fit_single(img)
    keep = np.ones(N * N)
    for _ in range(2 if hits else 0):
        new = hit_mask(img, img.ravel() - out[2].ravel(), [out[0][:2]]) * keep
        if np.array_equal(new, keep):
            break
        keep = new
        out = fit_single(img, keep=keep)
    return out


def fit_target(imgs, keeps=None):
    """Single-star fits, then a companion search on the aligned mean residual and a joint two-star refit; returns the
    target's flux and the model of every stamp."""
    keeps = keeps or [np.ones(N * N) for _ in imgs]
    singles = [fit_single(im, model_weights=False, keep=kp) for im, kp in zip(imgs, keeps)]
    stack = np.full((len(imgs), 2 * N, 2 * N), np.nan)
    for k, ((q, f, r, w), im) in enumerate(zip(singles, imgs)):
        sx, sy = int(round(N - q[0])), int(round(N - q[1]))
        stack[k, sy:sy + N, sx:sx + N] = r / max(f, 1.0)
    with warnings.catch_warnings():                       # cells outside every stamp are all-NaN
        warnings.simplefilter("ignore", RuntimeWarning)
        acc = np.nan_to_num(np.nanmedian(stack, axis=0), nan=-np.inf)  # median: one stamp's hit does not lead it
    yy, xx = np.mgrid[0:2 * N, 0:2 * N]
    ring = (np.hypot(xx - N, yy - N) > 2.0) & (np.hypot(xx - N, yy - N) < 8.0)
    k = np.argmax(np.where(ring, acc, -np.inf))
    py, px = np.unravel_index(k, acc.shape)
    dx0, dy0 = px - N, py - N
    chi_single = sum(np.sum((kp * r.ravel() / np.sqrt(np.clip(im, 0, None).ravel() + 36.0)) ** 2)
                     for (_, _, r, _), im, kp in zip(singles, imgs, keeps))
    ws = [kp / np.sqrt(np.clip(im, 0, None).ravel() + 36.0) for im, kp in zip(imgs, keeps)]

    def res(q):
        dx, dy = q[-2:]
        out = []
        for i, im in enumerate(imgs):
            x0, y0, fw = q[3 * i:3 * i + 3]
            _, model = linear_solve([psf(x0, y0, fw), psf(x0 + dx, y0 + dy, fw)], im, ws[i])
            out.append((model - im.ravel()) * ws[i])
        return np.concatenate(out)
    sol = None
    for sgn in (1.0, -1.0):                                  # a blend is nearly symmetric: try the mirrored offset too
        q0 = np.r_[np.concatenate([s[0] for s in singles]), sgn * dx0, sgn * dy0]
        trial = least_squares(res, q0, x_scale="jac")
        if sol is None or trial.cost < sol.cost:
            sol = trial
    if chi_single - 2 * sol.cost < 25.0 * len(imgs):        # no companion worth a second star
        return ([f for _, f, _, _ in singles], [im.ravel() - r.ravel() for (_, _, r, _), im in zip(singles, imgs)],
                [[q[:2]] for q, _, _, _ in singles])
    q = sol.x
    fluxes, models, centres = [], [], []
    for i, im in enumerate(imgs):
        x0, y0, fw = q[3 * i:3 * i + 3]
        coef, model = linear_solve([psf(x0, y0, fw), psf(x0 + q[-2], y0 + q[-1], fw)], im, ws[i])
        fluxes.append(max(coef[0], coef[1]))                 # the target is the brighter of the two in every stamp
        models.append(model)
        centres.append([(x0, y0), (x0 + q[-2], y0 + q[-1])])
    return fluxes, models, centres


def fit_target_clean(imgs, hits=True):
    """fit_target, then refitted twice without the pixels its models cannot explain."""
    fluxes, models, centres = fit_target(imgs)
    keeps = [np.ones(N * N) for _ in imgs]
    for _ in range(2 if hits else 0):
        new = [hit_mask(im, md, ct) * kp for im, md, ct, kp in zip(imgs, models, centres, keeps)]
        if all(np.array_equal(a, b) for a, b in zip(new, keeps)):
            break
        keeps = new
        fluxes, models, centres = fit_target(imgs, keeps)
    return fluxes


def load(folder, linear=True, companion=True, hits=True):
    with open(folder / "frames.csv", newline="") as fh:
        frames = list(csv.DictReader(fh))
    head = json.loads((folder / "epoch.json").read_text())
    stamps = np.load(folder / "stamps.npy").astype(float) - head["bias_pedestal_adu"]
    if linear:
        stamps = stamps * (1 + NONLIN * np.clip(stamps, 0, None))
    target = head["target"]
    tidx = [i for i, f in enumerate(frames) if f["star"] == target]
    for i, f in enumerate(frames):
        if f["star"] != target:
            f["net_counts"] = str(fit_clean(stamps[i], hits)[1])
    tf = (fit_target_clean([stamps[i] for i in tidx], hits) if companion
          else [fit_clean(stamps[i], hits)[1] for i in tidx])
    for i, fl in zip(tidx, tf):
        frames[i]["net_counts"] = str(fl)
    cat = {s["star"]: (s["V"], s["B_V"]) for s in json.loads((folder / "standards.json").read_text())}
    return frames, cat, target


def shutter_delta(frames, cat):
    by = {}
    for f in frames:
        if f["star"] in cat:
            by.setdefault((f["star"], f["filter"]), []).append(f)
    deltas, weights = [], []
    for group in by.values():
        group.sort(key=lambda f: float(f["ut_hours"]))
        for a, b in zip(group, group[1:]):
            ta, tb = float(a["exposure_s"]), float(b["exposure_s"])
            if abs(float(b["ut_hours"]) - float(a["ut_hours"])) > 0.05 or abs(tb / ta - 10) > 0.5:
                continue
            ca, cb = float(a["net_counts"]), float(b["net_counts"])
            deltas.append((ca * tb - cb * ta) / (cb - ca))
            weights.append(ca)
    order = np.argsort(deltas)
    cw = np.cumsum(np.array(weights)[order])
    return float(np.array(deltas)[order][np.searchsorted(cw, cw[-1] / 2)])


def reduce_night(folder, shutter=True, drift=True, linear=True, logged_shutter=False, companion=True, hits=True,
                 discordant=True):
    frames, cat, target = load(Path(folder), linear, companion, hits)
    if logged_shutter:
        delta = json.loads((Path(folder) / "epoch.json").read_text())["shutter_correction_s"]
    else:
        delta = shutter_delta(frames, cat) if shutter else 0.0
    t_first = min(float(f["ut_hours"]) for f in frames)
    coef = {}
    for filt in ("B", "V"):
        rows, rhs = [], []
        for f in frames:
            if f["filter"] != filt or f["star"] not in cat:
                continue
            v, c = cat[f["star"]]
            X, t = float(f["airmass"]), float(f["ut_hours"]) - t_first
            m = -2.5 * np.log10(float(f["net_counts"]) / (float(f["exposure_s"]) + delta))
            rows.append([1.0, X, t * X if drift else 0.0, c * X, c])
            rhs.append(m - (v + c if filt == "B" else v))
        a_mat, y_vec = np.array(rows), np.array(rhs)
        keep = np.ones(len(y_vec), bool)
        for _ in range(5):                                   # 4-sigma clipping of discordant standard frames
            coef[filt] = np.linalg.lstsq(a_mat[keep], y_vec[keep], rcond=None)[0]
            res = y_vec - a_mat @ coef[filt]
            sig = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
            new = np.abs(res) < 4 * max(sig, 1e-4)
            if np.array_equal(new, keep):
                break
            keep = new
    a, b = [], []
    for f in frames:
        if f["star"] != target:
            continue
        z, k1, r, k2, eps = coef[f["filter"]]
        X, t = float(f["airmass"]), float(f["ut_hours"]) - t_first
        m = -2.5 * np.log10(float(f["net_counts"]) / (float(f["exposure_s"]) + delta))
        known = z + k1 * X + (r * t * X if drift else 0.0)
        # m - known = V + (1 if B else 0) c + (k2 X + eps) c
        a.append([1.0, (1.0 if f["filter"] == "B" else 0.0) + k2 * X + eps])
        b.append(m - known)
    a, b = np.array(a), np.array(b)
    filt = np.array([f["filter"] for f in frames if f["star"] == target])
    use = np.ones(len(b), bool)
    while True:
        v, c = np.linalg.lstsq(a[use], b[use], rcond=None)[0]
        res = np.abs(b - a @ [v, c])
        worst = int(np.argmax(np.where(use, res, -1.0)))
        # a target frame more than 0.05 mag off the others (the wheel left in the other filter, or cloud) is dropped,
        # keeping at least two frames in each filter
        if not discordant or res[worst] < 0.05 or np.sum(use & (filt == filt[worst])) <= 2:
            break
        use[worst] = False
    return {"target_v_mag": float(v), "target_b_minus_v": float(c)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_night(sys.argv[1])))
