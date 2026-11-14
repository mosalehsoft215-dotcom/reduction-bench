"""Independent reduction of one film, written apart from solution/film_modes.py.

Cut-off and dip centres: copies of the reference's routines (an envelope-consistent square-root onset; Lorentzian
dip fits), so that this program stands alone.  Mode numbers: every pair of top-dip orders for TE and TM is tried with a forward mode
solver and Nelder-Mead.  Coupling across the air gap is then taken into account: for trial film, gap and stage zero,
the resonance next to each mode is the pole of a rational fit to the stack's complex reflection coefficient sampled
along real N, and those resonances are fitted to the measured dip centres.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import percentile_filter
from scipy.optimize import brentq, least_squares, minimize, minimize_scalar
from scipy.signal import find_peaks

NP, A_DEG, WAVELENGTH = 1.9648, 50.0, 0.6328
EDGE = 0.02
NARROW = 2.0
STEP = 3e-5


def to_index(theta):
    return NP * np.sin(math.radians(A_DEG) + np.arcsin(np.sin(np.radians(theta)) / NP))


def to_angle(n):
    return math.degrees(math.asin(NP * math.sin(math.asin(n / NP) - math.radians(A_DEG))))


def forward(n_f, d, n_s, pol, count, n_e=None):
    """Free-slab modes; TM light in a uniaxial film (optic axis along the normal) has kz = (n_o / n_e) sqrt(n_e^2 - N^2)
    and admittance kz / n_o^2."""
    n_e = n_f if n_e is None else n_e
    k = 2 * math.pi / WAVELENGTH
    out = []
    for m in range(count):
        def f(N):
            if pol == "TE":
                g, rs, rc = math.sqrt(max(n_f ** 2 - N ** 2, 1e-14)), 1.0, 1.0
            else:
                g = n_f / n_e * math.sqrt(max(n_e ** 2 - N ** 2, 1e-14))
                rs, rc = (n_f / n_s) ** 2, n_f ** 2
            return (k * d * g - m * math.pi - math.atan(rs * math.sqrt(N ** 2 - n_s ** 2) / g)
                    - math.atan(rc * math.sqrt(N ** 2 - 1) / g))
        lo, hi = n_s + 1e-10, (n_f if pol == "TE" else n_e) - 1e-10
        out.append(brentq(f, lo, hi, xtol=1e-12) if f(lo) > 0 else n_s)
    return np.array(out)


# The reference's cut-off fit, copied so that this program stands alone.
def knee_angle(theta, level, n_sub):
    """Stage angle of the substrate cut-off.  Above it the stack reflects totally (apart from the mode dips); below it
    light leaks into the substrate and the level starts to fall as sqrt(n_s^2 - N^2), with broad leaky resonances
    further down.  Scanning down from high angles, the cut-off is the first place where the upper envelope stays
    clearly below the plateau above it; a local square-root onset then refines it."""
    expect = to_angle(n_sub)
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
        n = to_index(t - zero)
        curve = plateau * (1 - c * np.sqrt(np.clip(n_sub ** 2 - n * n, 0, None)))
        return percentile_filter(curve, 90, size=2 * half + 1, mode="nearest")[near]
    guess = first + 0.01 - expect
    fit = least_squares(lambda p: model(p) - env[near], [guess, float(np.median(env[(t > first + 0.02) & (t < first + 0.25)])), 1.0],
                        bounds=([guess - 0.08, 0, 0], [guess + 0.08, 3, 80]), x_scale=[0.005, 0.01, 0.1])
    return expect + fit.x[0], fit.x[2]


def cut_off(path, n_s):
    t, y = np.loadtxt(path, delimiter=",", skiprows=1, unpack=True)
    return knee_angle(t, y, n_s)


# The reference's Lorentzian dip fit, copied so that this program stands alone.
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


def centroids(path, start):
    t, y = np.loadtxt(path, delimiter=",", skiprows=1, unpack=True)
    rows = dips(t, y, start)
    return rows[np.argsort(rows[:, 0])[::-1]]


def stack_amplitude(N, n_f, d, gap, n_s, pol, n_e=None):
    """Complex reflection coefficient of the lossless pressed stack (plane-wave Airy sums), uniaxial film."""
    n_e = n_f if n_e is None else n_e
    k = 2 * math.pi / WAVELENGTH
    layers = [NP, 1.0, n_f, n_s]
    z = []
    for n in layers:
        w = np.sqrt(n * n - N * N + 0j)
        z.append(k * np.where(np.imag(w) >= 0, w, -w))
    q = [zj / (n * n) if pol == "TM" else zj for zj, n in zip(z, layers)]
    if pol == "TM":
        w = np.sqrt(n_e * n_e - N * N + 0j)
        z[2] = k * n_f / n_e * np.where(np.imag(w) >= 0, w, -w)
        q[2] = z[2] / (n_f * n_f)
    r = lambda a, b: (q[a] - q[b]) / (q[a] + q[b])
    phase = np.exp(2j * z[2] * d)
    inner = (r(1, 2) + r(2, 3) * phase) / (1 + r(1, 2) * r(2, 3) * phase)
    outer = inner * np.exp(2j * z[1] * gap)
    return (r(0, 1) + outer) / (1 + r(0, 1) * outer)


def resonance(N0, n_f, d, gap, n_s, pol, n_e=None):
    """Pole of r(N) next to a free mode, from a rational fit r ~ a + b / (N - p) to samples along real N
    (r N = p r + a N + c, linear in p, a and c), re-centred once on the fitted pole."""
    centre, width = N0, 3e-4
    for _ in range(5):
        N = centre + width * np.linspace(-2, 2, 9)
        r = stack_amplitude(N, n_f, d, gap, n_s, pol, n_e)
        design = np.column_stack((r, N, np.ones_like(N)))
        p, a, c = np.linalg.lstsq(design.astype(complex), r * N, rcond=None)[0]
        centre, width = p.real, max(NARROW * abs(p.imag), 1e-5)
    return p.real


def reduce(folder):
    folder = Path(folder)
    n_s = json.loads((folder / "film.json").read_text())["substrate_index"]
    cuts = {pol: cut_off(folder / f"scan_{pol}.csv", n_s) for pol in ("TE", "TM")}
    edge = max(cuts.values(), key=lambda c: c[1])[0]       # the clearer cut-off of the one mounting
    offset = edge - to_angle(n_s)
    rows = {pol: centroids(folder / f"scan_{pol}.csv", edge) for pol in ("TE", "TM")}
    seen = {pol: r[:, 0] for pol, r in rows.items()}
    spread = {pol: np.hypot(r[:, 1], 0.001) for pol, r in rows.items()}
    top = max(to_index(v[0] - offset) for v in seen.values())

    def free_cost(p, first):
        n_f, n_e, d = p
        if n_f <= top_te or n_e <= top_tm or d <= 0:
            return 1e6
        total = 0.0
        for pol, a in seen.items():
            N = np.sort(to_index(a - offset))[::-1]
            model = forward(n_f, d, n_s, pol, first[pol] + N.size, n_e)[first[pol]:]
            total += float(((model - N) ** 2).sum())
        return total

    top_te = to_index(seen["TE"][0] - offset)
    top_tm = to_index(seen["TM"][0] - offset)
    best = None
    for a in range(9):
        for b in range(max(0, a - 2), min(9, a + 3)):
            first = {"TE": a, "TM": b}
            fit = minimize(free_cost, [top_te + 0.01, top_tm + 0.01, 3.0], args=(first,), method="Nelder-Mead",
                           options={"xatol": 1e-7, "fatol": 1e-15, "maxiter": 2500})
            if best is None or fit.fun < best[0]:
                best = (fit.fun, fit.x, first)
    _, (n_f, n_e, d), first = best

    # Gap coupling: measured centroids against the reflectance minima of the pressed stack.
    def minima(p):
        n, ne, t, gap, zero = p
        out = []
        for pol, a in seen.items():
            modes = forward(n, t, n_s, pol, first[pol] + a.size, ne)[first[pol]:]
            for N0, obs, s in zip(modes, a, spread[pol]):
                out.append((to_angle(resonance(N0, n, t, gap, n_s, pol, ne)) + zero - obs) / s)
        out.append((zero - offset) / EDGE)       # the edge keeps the zero near; the dips fix it
        return np.array(out)

    fit = least_squares(minima, [n_f, n_e, d, 0.08, offset], x_scale=[1e-4, 1e-4, 1e-3, 1e-2, 1e-3], diff_step=STEP,
                        bounds=([n_s + 1e-3, n_s + 1e-3, 0.3, 0.02, offset - 0.15],
                                [n_f + 0.03, n_e + 0.03, 40, 0.5, offset + 0.15]))
    return {"film_index": float(fit.x[0]), "film_thickness_um": float(fit.x[2])}


if __name__ == "__main__":
    out = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(out) + chr(10), encoding="utf-8")
