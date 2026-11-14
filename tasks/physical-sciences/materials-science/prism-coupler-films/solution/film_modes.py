"""Film index and thickness from the TE and TM angle scans of a prism coupler.

Run as: python3 film_modes.py <film folder> <result.json>

The scans are reflectance of the pressed stack: prism, the air gap left under the head, film and substrate.  Light
reaches the film by tunnelling across that gap, so each dip sits at the resonance of the coupled stack, not at the
guided mode of the free air / film / substrate slab; the thinner the gap, the wider the dips and the further they
move (TE and TM in opposite directions).  The reduction therefore
1. finds the substrate cut-off (the stage zero) and the dip centres, and numbers the dips from the top;
2. solves the free-slab mode conditions for a first index and thickness;
3. fits index, thickness, stage zero and gap to the dip centres through the complex resonances of the coupled stack;
4. refines everything on the whole measured curves with the exact stack reflectance, blurred by the beam divergence.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.special import erf

PRISM_N, PRISM_A, LAM_UM = 1.9648, 50.0, 0.6328
K = 2 * math.pi / LAM_UM
PRIOR = 0.02          # width (deg) of the cut-off prior on the stage zero


def effective(theta_deg):
    return PRISM_N * np.sin(np.radians(PRISM_A) + np.arcsin(np.sin(np.radians(theta_deg)) / PRISM_N))


def angle_of(n_eff):
    return np.degrees(np.arcsin(PRISM_N * np.sin(np.arcsin(n_eff / PRISM_N) - np.radians(PRISM_A))))


def load(path):
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1]


def knee_angle(theta, level, n_sub):
    """Stage angle of the substrate cut-off.  Above it the stack reflects totally (apart from the mode dips); below it
    light leaks into the substrate and the level starts to fall as sqrt(n_s^2 - N^2), with broad leaky resonances
    further down.  Scanning down from high angles, the cut-off is the first place where the upper envelope stays
    clearly below the plateau above it; a local square-root onset then refines it."""
    expect = angle_of(n_sub)
    step = theta[1] - theta[0]
    half = int(round(0.03 / step))
    sel = np.flatnonzero(np.abs(theta - expect) < 0.9)
    env = np.array([np.percentile(level[max(0, i - half):i + half + 1], 90) for i in sel])
    t = theta[sel]
    noise = 1.4826 * np.median(np.abs(np.diff(level))) / math.sqrt(2)
    first = None
    for k in range(len(t) - 1, -1, -1):
        if t[k] > expect + 0.6:
            continue
        above = env[(t > t[k] + 0.02) & (t < t[k] + 0.25)]
        below = env[(t >= t[k] - 0.12) & (t <= t[k])]
        if above.size < 10 or below.size < 12:
            continue
        plateau = np.median(above)
        if np.all(below < plateau * (1 - max(0.012, 4 * noise / plateau))):
            first = t[k]
            break
    if first is None:
        return expect, 0.0
    near = (t > first - 0.12) & (t < first + 0.15)
    from scipy.ndimage import percentile_filter

    def model(p):
        # the onset passed through the same running 90th percentile as the data, so the envelope's own smearing
        # of the kink does not move the fitted cut-off
        zero, plateau, c = p
        n = effective(t - zero)
        curve = plateau * (1 - c * np.sqrt(np.clip(n_sub ** 2 - n * n, 0, None)))
        return percentile_filter(curve, 90, size=2 * half + 1, mode="nearest")[near]
    guess = first + 0.01 - expect
    fit = least_squares(lambda p: model(p) - env[near], [guess, float(np.median(env[(t > first + 0.02) & (t < first + 0.25)])), 1.0],
                        bounds=([guess - 0.08, 0, 0], [guess + 0.08, 3, 80]), x_scale=[0.005, 0.01, 0.1])
    return expect + fit.x[0], fit.x[2]


def dips(theta, level, start):
    """Centres of the dips above the cut-off: prominent minima, each refined by a Lorentzian on a straight local
    baseline over a window that stops short of its neighbours."""
    from scipy.signal import find_peaks
    keep = theta > start + 0.03
    t, y = theta[keep], level[keep]
    noise = 1.4826 * np.median(np.abs(np.diff(y))) / math.sqrt(2)
    smooth = np.convolve(y, np.ones(3) / 3, mode="same")
    peaks, props = find_peaks(-smooth, prominence=max(8 * noise, 0.015), distance=4, wlen=161)
    centres = []
    for n, i in enumerate(peaks):
        left = (peaks[n - 1] + i) // 2 if n else max(0, i - 30)
        right = (peaks[n + 1] + i) // 2 if n + 1 < len(peaks) else min(y.size - 1, i + 30)
        left, right = max(left, i - 30), min(right, i + 30)
        tt, yy = t[left:right + 1], y[left:right + 1]
        lor = lambda p: p[3] + p[4] * (tt - t[i]) - p[1] * (p[2] / 2) ** 2 / ((tt - p[0]) ** 2 + (p[2] / 2) ** 2)
        try:
            fit = least_squares(lambda p: lor(p) - yy, [t[i], props["prominences"][n], 0.02, yy.max(), 0.0],
                                bounds=([tt[0], 0, 0.002, 0, -5], [tt[-1], 2, 0.3, 3, 5]))
            var = fit.fun @ fit.fun / max(len(yy) - 5, 1)
            try:
                sigma = math.sqrt(max(np.linalg.inv(fit.jac.T @ fit.jac)[0, 0] * var, 0))
            except np.linalg.LinAlgError:
                sigma = 0.01
            centres.append((fit.x[0], sigma))
        except ValueError:
            centres.append((t[i], 0.01))
    return np.array(centres).reshape(-1, 2)


def film_kz(N, n_x, n_z, pol, n_y=None):
    """Normal wavenumber in the film (units of K) and the divisor of its admittance.  Obliquely sputtered films are
    biaxial with principal axes x (in plane, along the plane of incidence), y (in plane, across it) and z (normal):
    TE light is polarised along y and sees n_y; TM light has kz = (n_x / n_z) sqrt(n_z^2 - N^2) and admittance
    kz / n_x^2."""
    if pol == "TE":
        n_y = n_x if n_y is None else n_y
        return np.sqrt(n_y * n_y - N * N + 0j), 1.0
    return n_x / n_z * np.sqrt(n_z * n_z - N * N + 0j), n_x * n_x


def mode_phase(N, m, n_f, d, n_s, pol, n_e=None, n_y=None):
    n_e = n_f if n_e is None else n_e
    g, eps = film_kz(N, n_f, n_e, pol, n_y)
    g = np.maximum(np.real(g), 1e-6)
    rs, rc = ((eps / n_s ** 2), eps) if pol == "TM" else (1.0, 1.0)
    return (K * d * g - m * math.pi - np.arctan(rs * np.sqrt(np.maximum(N ** 2 - n_s ** 2, 0)) / g)
            - np.arctan(rc * np.sqrt(N ** 2 - 1) / g))


def free_modes(n_f, d, n_s, pol, n_e=None, n_y=None):
    n_e = n_f if n_e is None else n_e
    n_y = n_f if n_y is None else n_y
    top = n_y if pol == "TE" else n_e
    out, m = [], 0
    while mode_phase(n_s + 1e-9, m, n_f, d, n_s, pol, n_e, n_y) > 0:
        out.append(brentq(lambda N: mode_phase(N, m, n_f, d, n_s, pol, n_e, n_y), n_s + 1e-9, top - 1e-9,
                          xtol=1e-14))
        m += 1
    return np.array(out)


def kz(n, N):
    z = np.sqrt(n * n - N * N + 0j)
    return K * np.where(np.imag(z) < 0, -z, z)


def stack_terms(N, n_f, kappa, d, gap, n_s, pol, n_e=None, n_y=None):
    n_e = n_f if n_e is None else n_e
    n_y = n_f if n_y is None else n_y
    layers = [PRISM_N, 1.0, n_f + 1j * kappa, n_s]
    z = [kz(n, N) for n in layers]
    q = [zj / (n * n) if pol == "TM" else zj for zj, n in zip(z, layers)]
    g, eps = film_kz(N, n_f + 1j * kappa, n_e + 1j * kappa, pol, n_y + 1j * kappa)
    z[2] = K * g
    q[2] = z[2] / eps
    r = lambda i, j: (q[i] - q[j]) / (q[i] + q[j])
    phase_f = np.exp(2j * z[2] * d)
    lower = (r(1, 2) + r(2, 3) * phase_f) / (1 + r(1, 2) * r(2, 3) * phase_f)
    coupling = lower * np.exp(2j * z[1] * gap)
    return r(0, 1), coupling


def reflectance(N, n_f, kappa, d, gap, n_s, pol, n_e=None, n_y=None):
    top, coupling = stack_terms(N, n_f, kappa, d, gap, n_s, pol, n_e, n_y)
    return np.abs((top + coupling) / (1 + top * coupling)) ** 2


def pole(guess, n_f, d, gap, n_s, pol, n_e=None, n_y=None):
    """Complex resonance of the pressed lossless stack next to a free-slab mode.  With the lower reflection written
    r = A / B (A = r12 + r23 e^{2i kz_f d}, B = 1 + r12 r23 e^{2i kz_f d}), the resonance is B + r01 e^{2i kz_g g} A = 0,
    which has no pole at the free mode (where B = 0)."""
    n_e = n_f if n_e is None else n_e
    n_y = n_f if n_y is None else n_y

    def f(N):
        layers = [PRISM_N, 1.0, n_f, n_s]
        z = [K * np.sqrt(n * n - N * N + 0j) for n in layers]
        z[0] = z[0] if z[0].real >= 0 else -z[0]            # light leaves into the prism
        z[3] = z[3] if z[3].imag >= 0 else -z[3]            # and decays into the substrate
        q = [zj / (n * n) if pol == "TM" else zj for zj, n in zip(z, layers)]
        g, eps = film_kz(N, n_f, n_e, pol, n_y)
        z[2] = K * g
        q[2] = z[2] / eps
        r = lambda i, j: (q[i] - q[j]) / (q[i] + q[j])
        pf = np.exp(2j * z[2] * d)
        return (1 + r(1, 2) * r(2, 3) * pf) + r(0, 1) * np.exp(2j * z[1] * gap) * (r(1, 2) + r(2, 3) * pf)
    a, b = complex(guess), complex(guess) * (1 + 1e-7) + 1e-7j
    fa, fb = f(a), f(b)
    for _ in range(60):
        c = b - fb * (b - a) / (fb - fa)
        a, fa, b, fb = b, fb, c, f(c)
        if abs(b - a) < 1e-15:
            break
    return complex(b)


def first_guess(lines, n_sub):
    """Free-slab index and thickness without numbering the dips: on a grid, every observed dip must sit next to a
    mode and every mode well above the cut-off next to a dip; then a least-squares polish with that numbering."""
    top = max(v.max() for v in lines.values())
    best = (np.inf, None)
    d_grid = np.geomspace(0.4, 25.0, 2400)
    for n in np.arange(top + 1e-4, top + 0.05, 1e-4):
        cost = np.zeros_like(d_grid)
        for pol, obs in lines.items():
            N = np.linspace(n_sub + 1e-7, n - 1e-7, 700)
            rs, rc = ((n / n_sub) ** 2, n ** 2) if pol == "TM" else (1.0, 1.0)
            g = np.sqrt(n * n - N * N)
            phi = np.arctan(rs * np.sqrt(N * N - n_sub ** 2) / g) + np.arctan(rc * np.sqrt(N * N - 1) / g)
            pred = []
            for m in range(60):
                dm = (m * math.pi + phi) / (K * g)          # thickness that puts mode m at N (rising with N)
                if dm[0] > d_grid[-1]:
                    break
                pred.append(np.where(d_grid >= dm[0], np.interp(d_grid, dm, N), np.nan))
            pred = np.array(pred)
            gap_obs = np.nanmin(np.abs(pred[:, :, None] - obs[None, None, :]), axis=0)        # (d, obs)
            cost += np.minimum(np.nan_to_num(gap_obs, nan=1.0), 2e-3).sum(axis=1)
            seen = pred > n_sub + 3e-3
            gap_pred = np.min(np.abs(pred[:, :, None] - obs[None, None, :]), axis=2)
            cost += np.where(seen, np.minimum(gap_pred, 2e-3), 0).sum(axis=0)
        k = int(np.argmin(cost))
        if cost[k] < best[0]:
            best = (cost[k], (n, d_grid[k]))
    n0, d0 = best[1]
    order = {pol: [int(np.argmin(np.abs(free_modes(n0, d0, n_sub, pol) - N))) for N in obs] for pol, obs in lines.items()}

    def residual(p):
        return np.concatenate([mode_phase(N, np.array(order[pol]), p[0], p[1], n_sub, pol) for pol, N in lines.items()])
    fit = least_squares(residual, [n0, d0], bounds=([top + 1e-6, 0.2], [top + 0.2, 40.0]), x_scale=[1e-4, 1e-2])
    return fit.x, order


def numbered_fit(angles, offset, n_sub, first=None, model="biaxial"):
    """Free-slab fit with the dips above the cut-off numbered consecutively from the top; first[pol] is the mode
    number of the top dip.  model: biaxial (n_x, n_y, n_z), uniaxial (n_y = n_x) or isotropic.
    Returns (cost, n_x, n_y, n_z, d)."""
    lines = {}
    for pol, a in angles.items():
        N = np.sort(effective(a - offset))[::-1]
        lines[pol] = N[N > n_sub + 2e-4]
        if lines[pol].size < 2:
            return np.inf, None, None, None, None
    first = first or {pol: 0 for pol in lines}

    def unpack(p):
        if model == "isotropic":
            return p[0], p[0], p[0], p[1]
        if model == "uniaxial":
            return p[0], p[0], p[1], p[2]
        return p[0], p[1], p[2], p[3]

    def residual(p):
        n_x, n_y, n_z, d = unpack(p)
        return np.concatenate([mode_phase(N, first[pol] + np.arange(N.size), n_x, d, n_sub, pol, n_z, n_y)
                               for pol, N in lines.items()])
    te, tm = lines["TE"][0], lines["TM"][0]
    best = None
    for d0 in (0.8, 1.6, 3.0, 5.0, 9.0):
        if model == "isotropic":
            top = max(te, tm)
            x0, lo, hi = [top + 0.003, d0], [top + 1e-6, 0.2], [top + 0.25, 40.0]
        elif model == "uniaxial":
            x0, lo, hi = [te + 0.003, tm + 0.003, d0], [te + 1e-6, tm + 1e-6, 0.2], [te + 0.25, tm + 0.25, 40.0]
        else:
            x0 = [tm + 0.001, te + 0.003, tm + 0.006, d0]
            lo, hi = [n_sub + 1e-3, te + 1e-6, tm + 1e-6, 0.2], [tm + 0.25, te + 0.25, tm + 0.25, 40.0]
        try:
            fit = least_squares(residual, np.clip(x0, lo, hi), bounds=(lo, hi), x_scale=[1e-3] * (len(x0) - 1) + [0.1])
        except ValueError:
            continue
        if best is None or fit.cost < best.cost:
            best = fit
    if best is None:
        return np.inf, None, None, None, None
    return (best.cost, *unpack(best.x))


def staged_fit(angles, offset, n_sub):
    """Biaxial free-slab fit in stages: the TE ladder alone gives n_y and d (TE sees n_y only), the TM ladder with
    that thickness gives n_x and n_z, each with its own mode numbering; a joint polish follows.
    Returns (cost, n_x, n_y, n_z, d)."""
    lines = {}
    for pol, a in angles.items():
        N = np.sort(effective(a - offset))[::-1]
        lines[pol] = N[N > n_sub + 2e-4]
        if lines[pol].size < 2:
            return np.inf, None, None, None, None
    te, tm = lines["TE"], lines["TM"]
    best_te = None
    for a in range(9):
        for d0 in (0.8, 1.6, 3.0, 5.0, 9.0):
            res = lambda p: mode_phase(te, a + np.arange(te.size), p[0], p[1], n_sub, "TE")
            try:
                fit = least_squares(res, [te[0] + 0.003, d0], bounds=([te[0] + 1e-6, 0.2], [te[0] + 0.25, 40.0]),
                                    x_scale=[1e-3, 0.1])
            except ValueError:
                continue
            if best_te is None or fit.cost < best_te[0].cost:
                best_te = (fit, a)
    (n_y, d), a = best_te[0].x, best_te[1]
    best_tm = None
    for b in range(9):
        for dz in (0.002, 0.01, 0.03):
            res = lambda p: mode_phase(tm, b + np.arange(tm.size), p[0], d, n_sub, "TM", p[1])
            x0 = [tm[0] + dz, tm[0] + dz]
            lo, hi = [n_sub + 1e-3, tm[0] + 1e-6], [tm[0] + 0.3, tm[0] + 0.25]
            try:
                fit = least_squares(res, x0, bounds=(lo, hi), x_scale=[1e-3, 1e-3])
            except ValueError:
                continue
            if best_tm is None or fit.cost < best_tm[0].cost:
                best_tm = (fit, b)
    (n_x, n_z), b = best_tm[0].x, best_tm[1]

    def residual(p):
        return np.concatenate([mode_phase(te, a + np.arange(te.size), p[0], p[3], n_sub, "TE", p[2], p[1]),
                               mode_phase(tm, b + np.arange(tm.size), p[0], p[3], n_sub, "TM", p[2], p[1])])
    x0 = [n_x, n_y, n_z, d]
    lo = [n_sub + 1e-3, te[0] + 1e-6, tm[0] + 1e-6, 0.2]
    hi = [n_x + 0.05, te[0] + 0.25, tm[0] + 0.25, 40.0]
    fit = least_squares(residual, np.clip(x0, lo, hi), bounds=(lo, hi), x_scale=[1e-3, 1e-3, 1e-3, 0.1])
    return (fit.cost, *fit.x)


def coupled_fit(found, offset, start, n_sub, model="biaxial"):
    """Principal indices, thickness, stage zero and gap from the dip centres through the coupled-stack resonances."""
    n_x0, n_y0, n_z0, d0 = start
    keep, order = {}, {}
    for pol, rows in found.items():
        N = effective(rows[:, 0] - offset)
        above = N > n_sub + 2e-4
        keep[pol] = rows[above][np.argsort(N[above])[::-1]]
        free = free_modes(n_x0, d0, n_sub, pol, n_z0, n_y0)
        order[pol] = [int(np.argmin(np.abs(free - v))) for v in np.sort(N[above])[::-1]]
    count = sum(len(v) for v in keep.values())
    width = {"isotropic": 1, "uniaxial": 2, "biaxial": 3}[model]

    def unpack(p):
        idx, (t, zero, gap) = p[:width], p[width:]
        if model == "isotropic":
            return idx[0], idx[0], idx[0], t, zero, gap
        if model == "uniaxial":
            return idx[0], idx[0], idx[1], t, zero, gap
        return idx[0], idx[1], idx[2], t, zero, gap

    def residual(p):
        n_x, n_y, n_z, t, zero, gap = unpack(p)
        out = []
        for pol, rows in keep.items():
            free = free_modes(n_x, t, n_sub, pol, n_z, n_y)
            if free.size <= max(order[pol]):
                return np.full(count + 1, 1e3)
            res = np.array([pole(free[m], n_x, t, gap, n_sub, pol, n_z, n_y).real for m in order[pol]])
            sigma = np.hypot(rows[:, 1], 0.0008)
            out.append((angle_of(np.clip(res, n_sub, PRISM_N * 0.999)) + zero - rows[:, 0]) / sigma)
        out.append([(zero - offset) / PRIOR])                  # the cut-off: a fair but rounded edge
        return np.concatenate(out)

    idx0 = {"isotropic": [n_x0], "uniaxial": [n_x0, n_z0], "biaxial": [n_x0, n_y0, n_z0]}[model]
    best = None
    for gap0 in (0.05, 0.08, 0.13, 0.22):
        x0 = idx0 + [d0, offset, gap0]
        lo = [n_sub + 1e-4] * width + [0.2, offset - 0.15, 0.02]
        hi = [v + 0.03 for v in idx0] + [40, offset + 0.15, 0.6]
        try:
            fit = least_squares(residual, x0, x_scale=[1e-4] * width + [1e-3, 1e-3, 1e-2], diff_step=1e-7,
                                bounds=(lo, hi))
        except ValueError:
            continue
        if best is None or fit.cost < best.cost:
            best = fit
    return unpack(best.x)


def measure(folder, *, coupled=True, numbering=True, use_knee=True, model="biaxial"):
    folder = Path(folder)
    n_sub = json.loads((folder / "film.json").read_text())["substrate_index"]
    found, knees = {}, []
    for pol in ("TE", "TM"):
        theta, level = load(folder / f"scan_{pol}.csv")
        found[pol] = dips(theta, level, theta[0] - 1)
        knees.append(knee_angle(theta, level, n_sub)[::-1])
    # One mounting for both scans: the clearer substrate cut-off of the two fixes the stage zero.
    offset = max(knees)[1] - angle_of(n_sub) if use_knee else 0.0
    angles = {pol: rows[:, 0] for pol, rows in found.items()}
    if model == "biaxial" and numbering:
        _, n_x, n_y, n_z, d = staged_fit(angles, offset, n_sub)
    elif numbering:
        # Mode numbers: every pair of top-dip orders for TE and TM; only the right pair lets one thickness carry
        # both polarisations (TE through n_y, TM through n_x and n_z).
        trials = [numbered_fit(angles, offset, n_sub, {"TE": a, "TM": b}, model) for a in range(9) for b in range(9)]
        _, n_x, n_y, n_z, d = min(trials, key=lambda x: x[0])
    else:
        _, n_x, n_y, n_z, d = numbered_fit(angles, offset, n_sub, None, model)
    if not coupled:
        return {"film_index": float(n_x), "film_thickness_um": float(d)}
    n_x, n_y, n_z, d, _, _ = coupled_fit(found, offset, (n_x, n_y, n_z, d), n_sub, model)
    return {"film_index": float(n_x), "film_thickness_um": float(d)}


if __name__ == "__main__":
    result = measure(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result, allow_nan=False) + "\n", encoding="utf-8")
