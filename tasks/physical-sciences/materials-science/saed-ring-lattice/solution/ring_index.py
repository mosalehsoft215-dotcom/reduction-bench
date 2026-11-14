"""Cubic lattice parameter from an azimuthally averaged SAED ring profile.

Run as: python3 ring_index.py <pattern folder> <result.json>

Coarse powders (crystallites of 10 nm and more) give Bragg rings: pseudo-Voigt rings at K sqrt(s)/a over a smooth
background, with the gold islands' rings fixing the camera constant K of the refocused area (bragg_measure).
Ultrafine powders (a few nm) do not: each crystallite scatters as a whole, and the orientation-averaged intensity of a
sphere cut from the ideal lattice (Debye scattering equation) has maxima away from the Bragg radii.
For those the specimen enters the fit as one Debye column: its scattering-vector scale fixes K/a, the gold rings fix
K, and the crystallite is rebuilt at the fitted lattice parameter before a last refinement.
"""
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, nnls

from scipy.spatial.distance import cdist, pdist

# the electron scattering factors of /app/data/scattering_factors.json (f = sum a_k exp(-b_k s^2), angstrom)
TABLE = {"O": {"a": [0.455, 0.917, 0.472, 0.138], "b": [0.24, 1.79, 6.81, 20.3]}, "Mg": {"a": [0.21, 0.604, 1.267, 1.65], "b": [0.25, 1.95, 7.5, 32.0]}, "Ca": {"a": [0.38, 1.12, 2.25, 3.62], "b": [0.22, 1.88, 8.6, 41.0]}, "Mn": {"a": [0.35, 1.08, 2.07, 1.68], "b": [0.16, 1.39, 5.6, 21.6]}, "Co": {"a": [0.36, 1.12, 2.06, 1.45], "b": [0.15, 1.31, 5.3, 19.8]}, "Ni": {"a": [0.36, 1.13, 2.05, 1.38], "b": [0.15, 1.28, 5.2, 19.1]}, "Fe": {"a": [0.35, 1.1, 2.07, 1.56], "b": [0.15, 1.35, 5.45, 20.8]}, "Zr": {"a": [0.49, 1.52, 2.94, 3.04], "b": [0.17, 1.47, 6.4, 26.5]}, "Y": {"a": [0.48, 1.5, 2.96, 3.24], "b": [0.17, 1.5, 6.6, 28.0]}, "Ce": {"a": [0.55, 1.71, 3.38, 3.98], "b": [0.16, 1.44, 6.5, 28.4]}, "Au": {"a": [0.52, 1.66, 3.27, 3.15], "b": [0.12, 1.15, 5.1, 21.9]}}
SEARCH = {"rocksalt": (3.9, 5.0), "fluorite": (4.9, 5.7), "spinel": (7.9, 8.8)}

GOLD_A = 4.0786
KIND = {"cubic spinel (Fd-3m)": "spinel", "rock salt (Fm-3m)": "rocksalt", "fluorite (Fm-3m)": "fluorite"}


def reflections(kind, top=40):
    out = set()
    for h, k, l in itertools.product(range(0, 9), repeat=3):
        if (h, k, l) == (0, 0, 0) or not (h >= k >= l):
            continue
        s = h * h + k * k + l * l
        if s > top or len({h % 2, k % 2, l % 2}) > 1:
            continue
        if kind == "spinel" and h % 2 == 0 and s % 4 != 0:
            continue
        if kind == "fluorite" and h % 2 == 0 and (h + k + l) % 4 == 2:
            continue
        out.add(s)
    return sorted(out)


def shapes(r, centres, w0, w1, eta):
    cols = []
    for c in centres:
        x = (r - c) / (w0 + w1 * c)
        cols.append((1 - eta) * np.exp(-0.5 * x * x) + eta / (1 + x * x))
    return cols


def placed(c, distortion):
    """Where the projector lenses put a ring meant for radius c: radial distortion grows as the cube of radius."""
    c = np.asarray(c, dtype=float)
    return c * (1 + distortion * (c / 1000.0) ** 2)


def design(r, K, a, w, kind, gold, distortion=0.0):
    w0, w1, eta = w
    cols = shapes(r, placed([K * math.sqrt(s) / a for s in reflections(kind)], distortion), w0, w1, eta)
    if gold:
        cols += shapes(r, placed([K * math.sqrt(s) / GOLD_A for s in reflections("rocksalt")], distortion), w0, w1, eta)
    cols += [np.exp(-r / lam) for lam in (60.0, 120.0, 240.0, 480.0, 960.0)] + [np.ones_like(r)]
    return np.column_stack(cols)


def misfit(p, r, y, kind, gold, k_fixed=None):
    p = list(p)
    distortion = p.pop() if len(p) == (5 if k_fixed is not None else 6) else 0.0
    K, a, w0, w1, eta = (k_fixed, *p) if k_fixed is not None else p
    m = design(r, K, a, (w0, w1, eta), kind, gold, distortion)
    scale = np.sqrt(np.maximum(y, 1.0))
    coef, _ = nnls(m / scale[:, None], y / scale)
    return (m @ coef - y) / scale


def bragg_measure(folder, *, gold_rings=True, distortion=False):
    folder = Path(folder)
    meta = json.loads((folder / "pattern.json").read_text())
    kind = KIND[meta["phase"]]
    r, y = np.loadtxt(folder / "profile.csv", delimiter=",", skiprows=1, unpack=True)
    k0 = meta["camera_constant_A_px_session"]
    # first pass at the session camera constant: the strong sample rings fix K/a
    best = None
    rc, yc = r[::2], y[::2]                     # the coarse grid needs only every other point
    for a0 in np.arange(3.6, 9.2, 0.025):
        cost = float(np.sum(misfit([a0, 1.6, 0.0028, 0.4], rc, yc, kind, False, k0) ** 2))
        if best is None or cost < best[0]:
            best = (cost, a0)
    # NNLS lets unobserved rings take zero amplitude, so a doubled cell fits as well as the true one: refine the
    # best grid cell and its half, and keep the smaller cell unless it fits clearly worse
    refined = []
    for a_try in (best[1] / 2, best[1]):
        if a_try < 3.0:
            continue
        f_try = least_squares(misfit, [a_try, 1.6, 0.0028, 0.4], args=(r, y, kind, False, k0),
                              bounds=([3.0, 0.3, 0.0, 0.0], [10.0, 6.0, 0.02, 1.0]), x_scale=[1e-3, 0.05, 1e-4, 0.05])
        refined.append((2 * f_try.cost, f_try))
    fit = refined[0][1] if refined[0][0] <= refined[-1][0] * 1.01 else refined[-1][1]
    # The projector lenses may bend the outer rings outwards: free that distortion with the sample rings alone.
    bend = 0.03 if distortion else 1e-12
    fit = least_squares(misfit, [*fit.x, 0.0], args=(r, y, kind, False, k0),
                        bounds=([3.0, 0.3, 0.0, 0.0, -bend], [10.0, 6.0, 0.02, 1.0, bend]),
                        x_scale=[1e-3, 0.05, 1e-4, 0.05, 1e-4])
    a_session = fit.x[0]
    if not gold_rings:
        return {"lattice_parameter_A": float(a_session)}
    # The gold islands' rings sit at known spacings in every pattern: with them the camera constant of this very
    # pattern is fitted together with the lattice parameter (the ratio K/a stays pinned by the sample rings).
    ratio = k0 / a_session
    trials = []
    for f in np.arange(0.95, 1.0501, 0.0025):
        K = k0 * f
        p = [K, K / ratio, *fit.x[1:]]          # widths, shape and distortion carried over
        trials.append((float(np.sum(misfit(p, r, y, kind, True) ** 2)), p))
    start = min(trials, key=lambda t: t[0])[1]
    final = least_squares(misfit, start, args=(r, y, kind, True),
                          bounds=([k0 * 0.9, 3.0, 0.3, 0.0, 0.0, -bend], [k0 * 1.1, 10.0, 6.0, 0.02, 1.0, bend]),
                          x_scale=[0.5, 1e-3, 0.05, 1e-4, 0.05, 1e-4])
    return {"lattice_parameter_A": float(final.x[1])}


# Whole-crystallite (Debye) scattering
SITES = {
    "rocksalt": ([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
                 [[.5, 0, 0], [0, .5, 0], [0, 0, .5], [.5, .5, .5]]),
    "fluorite": ([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
                 [[x, y, z] for x in (.25, .75) for y in (.25, .75) for z in (.25, .75)]),
}


def factor(table, composition, s):
    """Composition-weighted electron scattering factor (A) at s (1/A)."""
    s2 = np.asarray(s) ** 2
    out = 0.0
    for element, share in composition.items():
        a, b = table[element]["a"], table[element]["b"]
        out = out + share * sum(ai * np.exp(-bi * s2) for ai, bi in zip(a, b))
    return out


def cluster(kind, a, diameter_A):
    cation, anion = (np.array(x, dtype=float) for x in SITES[kind])
    n = int(diameter_A / (2 * a)) + 2
    cells = np.array(list(itertools.product(range(-n, n + 1), repeat=3)), dtype=float)
    cat = ((cells[:, None, :] + cation[None]).reshape(-1, 3)) * a
    an = ((cells[:, None, :] + anion[None]).reshape(-1, 3)) * a
    r = diameter_A / 2
    return cat[np.linalg.norm(cat, axis=1) <= r + 1e-9], an[np.linalg.norm(an, axis=1) <= r + 1e-9]


def histograms(cat, an, width=0.002):
    out = []
    for kind, d in (("cc", pdist(cat)), ("aa", pdist(an)), ("ca", cdist(cat, an).ravel())):
        h, e = np.histogram(d, bins=np.arange(0.0, d.max() + 2 * width, width))
        keep = h > 0
        out.append((kind, 0.5 * (e[1:] + e[:-1])[keep], h[keep].astype(float)))
    return out


def intensity(s, cat, an, f_cat, f_an, hist=None):
    """Orientation-averaged intensity of the crystallite (f_cat, f_an evaluated at s)."""
    hist = hist or histograms(cat, an)
    total = len(cat) * f_cat ** 2 + len(an) * f_an ** 2
    weights = {"cc": 2 * f_cat * f_cat, "aa": 2 * f_an * f_an, "ca": 2 * f_cat * f_an}
    for kind, r, h in hist:
        x = 2 * np.pi * np.outer(s, r)
        total = total + weights[kind] * (np.sinc(x / np.pi) @ h)
    return total


def debye_table(kind, a, diameter_nm, composition):
    """Crystallite intensity on a fine scattering-vector grid at lattice parameter a (1/A grid, intensity)."""
    s = np.linspace(0.02, 1.3, 5200)
    cat, an = cluster(kind, a, diameter_nm * 10.0)
    f_cat = factor(TABLE, composition, s)
    f_an = factor(TABLE, {"O": 1.0}, s)
    return s, intensity(s, cat, an, f_cat, f_an)


def debye_columns(r, K, a, blur, table, a_table):
    s_ref, inten = table
    # distances scale with a: the crystallite at a scatters at s like the table crystallite at s a / a_table
    col = np.interp(r / K * a / a_table, s_ref, inten)
    step = r[1] - r[0]
    half = int(4 * blur / step) + 1
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) * step / blur) ** 2)
    return np.convolve(col, kernel / kernel.sum(), mode="same")


def debye_design(r, K, a, blur, table, a_table, gold_w=None):
    cols = [debye_columns(r, K, a, blur, table, a_table)]
    if gold_w is not None:
        w0, w1, eta = gold_w
        cols += shapes(r, [K * math.sqrt(s) / GOLD_A for s in reflections("rocksalt")], w0, w1, eta)
    cols += [np.exp(-r / lam) for lam in (60.0, 120.0, 240.0, 480.0, 960.0)] + [np.ones_like(r)]
    return np.column_stack(cols)


def solve(m, y):
    scale = np.sqrt(np.maximum(y, 1.0))
    a, b = m / scale[:, None], y / scale
    norms = np.linalg.norm(a, axis=0)
    try:
        coef, _ = nnls(a / norms, b, maxiter=20 * a.shape[1])
        coef = coef / norms
    except RuntimeError:
        coef = np.linalg.lstsq(a, b, rcond=None)[0]
    return (m @ coef - y) / scale


def debye_measure(folder, meta, kind):
    r, y = np.loadtxt(folder / "profile.csv", delimiter=",", skiprows=1, unpack=True)
    k0 = meta["camera_constant_A_px_session"]
    diameter, composition = meta["crystallite_diameter_nm"], meta["cation_composition"]
    lo, hi = SEARCH[kind]
    a_table = 0.5 * (lo + hi)
    table = debye_table(kind, a_table, diameter, composition)
    # 1. a coarse grid over K and a together (every fourth profile point), gold rings included
    rd, yd = r[::4], y[::4]
    gold_w = (1.6, 0.0028, 0.4)
    cells = []
    # the crystallite (fixed diameter) holds a different set of atoms at each a: tables every 2 % of a
    anchors = np.exp(np.arange(np.log(lo), np.log(hi) + 0.02, 0.02))
    tables = {float(x): debye_table(kind, float(x), diameter, composition) for x in anchors}
    for k in k0 * np.arange(0.95, 1.0501, 0.005):
        for a in np.exp(np.arange(np.log(lo), np.log(hi), 0.002)):
            at = float(anchors[int(np.argmin(np.abs(anchors - a)))])
            cells.append((float(np.sum(solve(debye_design(rd, k, a, 1.5, tables[at], at, gold_w), yd) ** 2)), k, a))
    cells.sort()
    picks, seen = [], []
    for c, k, a in cells:                      # three distinct starting cells
        if all(abs(k / k2 - 1) > 0.01 or abs(a / a2 - 1) > 0.01 for k2, a2 in seen):
            picks.append((k, a))
            seen.append((k, a))
        if len(picks) == 3:
            break
    best = None
    for K, a in picks:
        table = debye_table(kind, a, diameter, composition)
        a_table = a

        def resid(p, tab=table, at=a_table):
            k, aa, blur, w0, w1, eta = p
            return solve(debye_design(r, k, aa, blur, tab, at, (w0, w1, eta)), y)
        fit = least_squares(resid, [K, a, 1.5, *gold_w], x_scale=[0.2, 1e-4, 0.05, 0.05, 1e-4, 0.05],
                            bounds=([k0 * 0.9, a * 0.97, 0.3, 0.3, 0.0, 0.0], [k0 * 1.1, a * 1.03, 5.0, 6.0, 0.02, 1.0]))
        if best is None or fit.cost < best.cost:
            best = fit
    K, a, blur = best.x[:3]
    gold_w = tuple(best.x[3:])
    # 3. rebuild the crystallite at the fitted lattice parameter and refine again
    for _ in range(2):
        table = debye_table(kind, a, diameter, composition)
        a_table = a

        def resid(p, tab=table, at=a_table):
            k, aa, blur, w0, w1, eta = p
            return solve(debye_design(r, k, aa, blur, tab, at, (w0, w1, eta)), y)
        fit = least_squares(resid, [K, a, blur, *gold_w], x_scale=[0.2, 1e-4, 0.05, 0.05, 1e-4, 0.05],
                            bounds=([k0 * 0.9, a * 0.99, 0.3, 0.3, 0.0, 0.0], [k0 * 1.1, a * 1.01, 5.0, 6.0, 0.02, 1.0]))
        K, a, blur = fit.x[:3]
        gold_w = tuple(fit.x[3:])
    return {"lattice_parameter_A": float(a)}


def measure(folder, *, crystallites=True, gold_rings=True):
    folder = Path(folder)
    meta = json.loads((folder / "pattern.json").read_text())
    kind = KIND[meta["phase"]]
    if crystallites and meta.get("crystallite_diameter_nm", 100.0) < 10.0:
        return debye_measure(folder, meta, kind)
    return bragg_measure(folder, gold_rings=gold_rings)


if __name__ == "__main__":
    out = measure(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(out, allow_nan=False) + "\n", encoding="utf-8")
