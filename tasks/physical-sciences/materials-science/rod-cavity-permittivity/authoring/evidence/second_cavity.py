"""Independent reduction of one rod folder, written apart from solution/rod_eps.py.

Resonances: 1/|S21|^2 is exactly quadratic in frequency for a Lorentzian, so each sweep's upper part is fitted by
weighted linear least squares in (1, f, f^2), giving f0, the loaded Q and the peak in closed form.  The unloaded Q
follows from the peak of two equal probes.

Eigenvalues by point matching, without any addition theorem: around the rod centre the field is
sum_m [A_m J_m(k rho) + B_m Y_m(k rho)] cos(m phi) outside and C_m J_m(k sqrt(eps) rho) cos(m phi) inside; the two
rod-surface conditions fix B_m / A_m, and E_z = 0 is imposed at as many points of the wall's upper half as there are
A_m.  The resonance is a complex root of that determinant with eps'(1 - i tan delta), so its imaginary part is the
dielectric loss directly.  The wall loss of the loaded field comes from a boundary-integral (Rellich) form of the
stored energy, with derivatives on the wall taken by finite differences of the matched field.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import jv, jvp, yv, yvp, jn_zeros

ROOTS = jn_zeros(0, 2)
C_MM = 299792458.0e3
ETA = 376.730313668
HEIGHT = 10.0
TERMS = 24


def peak_parameters(path):
    rows = np.genfromtxt(path, delimiter=",", names=True)
    f = rows["freq_GHz"]
    power = 10 ** (rows["s21_dB"] / 10)
    keep = power > 0.3 * power.max()
    weight = power ** 3                       # amplitude noise makes var(1/P) grow as 1/P^3; weight it back
    for _ in range(3):                        # second and third passes: wider band, weights from the fitted curve
        x = f[keep] - f[keep].mean()
        design = np.column_stack((np.ones_like(x), x, x * x))
        root = np.sqrt(weight[keep])
        a, b, c = np.linalg.lstsq(design * root[:, None], (1 / power[keep]) * root, rcond=None)[0]
        centre = -b / (2 * c) + f[keep].mean()
        model = 1 / (a + b * (f - f[keep].mean()) + c * (f - f[keep].mean()) ** 2)
        keep, weight = model > 0.12 * model.max(), model ** 3
    floor = a - b * b / (4 * c)               # 1/|S21(f0)|^2
    f0 = centre
    q_loaded = f0 * np.sqrt(c / floor) / 2
    peak = 1 / np.sqrt(floor)
    return f0, q_loaded / (1 - peak)


def ratio(k, eps, b, m):
    """Scattered-to-regular ratio B_m / A_m of the rod; orders above 14 scatter nothing at these sizes."""
    m = np.asarray(m)
    out = np.zeros(m.shape, dtype=complex if np.iscomplexobj(k) or np.iscomplexobj(eps) else float)
    low = m <= 14
    out[low] = _ratio(k, eps, b, m[low])
    return out


def _ratio(k, eps, b, m):
    k1 = k * np.sqrt(eps)
    top = k * jvp(m, k * b) * jv(m, k1 * b) - k1 * jv(m, k * b) * jvp(m, k1 * b)
    return top / (k1 * yv(m, k * b) * jvp(m, k1 * b) - k * yvp(m, k * b) * jv(m, k1 * b))


def wall_points(a, d):
    theta = np.pi * (np.arange(TERMS) + 0.5) / TERMS
    x, y = a * np.cos(theta) - d, a * np.sin(theta)
    return np.hypot(x, y), np.arctan2(y, x)


def matching(k, eps, b, d, a):
    m = np.arange(TERMS)
    rho, phi = wall_points(a, d)
    t = ratio(k, eps, b, m)
    regular, scattered = jv(m[None, :], k * rho[:, None]), t * yv(m[None, :], k * rho[:, None])
    cols = (regular + scattered) * np.cos(m[None, :] * phi[:, None])
    scale = np.max(np.abs(regular) + np.abs(scattered), axis=0)      # never zero at a root, unlike max |column|
    return cols / scale, scale, t


def root(guess, eps, b, d, a):
    f = lambda k: np.linalg.det(matching(k, eps, b, d, a)[0])
    k0, k1 = complex(guess), complex(guess) * (1 - 5e-6)
    f0, f1 = f(k0), f(k1)
    for _ in range(60):
        k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
        k0, f0, k1, f1 = k1, f1, k2, f(k2)
        if abs(k1 - k0) < 1e-13 * abs(k1):
            break
    return k1


def field(k, eps, b, d, a):
    """Callable E_z(x, y) outside the rod and on its surface, normalised arbitrarily."""
    cols, scale, t = matching(k, eps, b, d, a)
    coeff = np.linalg.svd(cols)[2][-1] / scale
    m = np.arange(TERMS)

    def e(x, y):
        rho, phi = np.hypot(x - d, y), np.arctan2(y, x - d)
        return ((jv(m, k * rho[..., None]) + t * yv(m, k * rho[..., None])) * np.cos(m * phi[..., None])) @ coeff
    return e


def wall_q(k, eps, b, d, a, rs):
    e = field(k, eps, b, d, a)
    theta = np.linspace(0, 2 * np.pi, 800, endpoint=False)
    h = 1e-4 * a
    inner = e((a - h) * np.cos(theta), (a - h) * np.sin(theta))
    inner2 = e((a - 2 * h) * np.cos(theta), (a - 2 * h) * np.sin(theta))
    slope = (-4 * inner + inner2) / (2 * h) * -1           # one-sided derivative toward the wall (E = 0 there)
    w_side = np.mean(slope ** 2) * 2 * np.pi * a
    phi = np.linspace(0, 2 * np.pi, 800, endpoint=False)
    on_rod = e(d + b * np.cos(phi), b * np.sin(phi))
    s = a / (2 * k * k) * w_side + (eps - 1) / 2 * np.mean((d * np.cos(phi) + b) * on_rod ** 2) * 2 * np.pi * b
    return k ** 3 * ETA * HEIGHT * s / (rs * (HEIGHT * w_side + 2 * k * k * s))


def reduce(folder):
    folder = Path(folder)
    b = json.loads((folder / "rod.json").read_text())["rod_diameter_mm"] / 2
    seen = {}
    with open(folder / "sweeps.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            seen[(row["insertion"], int(row["mode"][3]) - 1)] = peak_parameters(folder / "sweeps" / row["sweep"])
    labels = sorted({k[0] for k in seen if k[0] != "empty"})
    radius = [ROOTS[n] * C_MM / (2 * math.pi * seen[("empty", n)][0] * 1e9) for n in range(2)]
    target = {(lab, n): ROOTS[n] / radius[n] * seen[(lab, n)][0] / seen[("empty", n)][0] for lab in labels for n in range(2)}

    def residual(p):
        return [root(target[(lab, n)], p[0], b, p[1 + j], radius[n]).real / target[(lab, n)] - 1
                for j, lab in enumerate(labels) for n in range(2)]
    # Start from first-order perturbation (shift = (eps-1)/2 (b/a)^2 J0^2(x d/a)/J1^2(x)), solved on a coarse grid.
    from scipy.special import j0, j1
    drop = np.array([[1 - target[(lab, n)] / (ROOTS[n] / radius[n]) for n in range(2)] for lab in labels])
    grid_d = np.linspace(0, 6, 241)
    ratio_model = (j0(ROOTS[1] * grid_d / 20) / j0(ROOTS[0] * grid_d / 20)) ** 2 * j1(ROOTS[0]) ** 2 / j1(ROOTS[1]) ** 2
    start_d = [grid_d[int(np.argmin(abs(ratio_model - drop[j, 1] / drop[j, 0])))] for j in range(len(labels))]
    start_e = 1 + np.mean([2 * drop[j, 0] * j1(ROOTS[0]) ** 2 / ((b / 20) ** 2 * j0(ROOTS[0] * start_d[j] / 20) ** 2)
                           for j in range(len(labels))])
    best = least_squares(residual, [start_e] + start_d, bounds=([1.01] + [0] * len(labels), [80] + [8] * len(labels)),
                         diff_step=1e-5, max_nfev=40)
    fit = least_squares(residual, best.x, bounds=([1.01] + [0] * len(labels), [80] + [8] * len(labels)),
                        diff_step=1e-7, xtol=1e-14, ftol=1e-15)
    eps = float(fit.x[0])
    # Loss: dielectric part from the complex root at a trial tan delta (it is linear in tan delta), walls apart.
    trial = 1e-3
    tangents = []
    for j, lab in enumerate(labels):
        d = fit.x[1 + j]
        for n in range(2):
            a = radius[n]
            k_empty = ROOTS[n] / a
            q_empty = seen[("empty", n)][1]
            rs = k_empty * ETA * a * HEIGHT / (2 * q_empty * (HEIGHT + a))
            k = root(target[(lab, n)], eps, b, d, a).real
            lossy = root(k, eps * (1 - 1j * trial), b, d, a)
            per_unit = 2 * lossy.imag / lossy.real / trial
            q_wall = wall_q(k, eps, b, d, a, rs * math.sqrt(k / k_empty))
            tangents.append((1 / seen[(lab, n)][1] - 1 / q_wall) / per_unit)
    return {"eps_real": eps, "loss_tangent": float(np.median(tangents))}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + chr(10), encoding="utf-8")
