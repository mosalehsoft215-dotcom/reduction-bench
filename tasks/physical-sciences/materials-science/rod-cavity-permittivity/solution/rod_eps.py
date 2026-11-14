"""Permittivity and loss tangent of a dielectric rod from TM010 and TM020 sweeps of a cylindrical cavity.

Run as: python3 rod_eps.py <rod folder> <answer.json>

The rod fills the full height, so the loaded fields do not vary with z and each resonance is an exact 2-D eigenvalue:
cavity Bessel waves J_n(kr) cos(n theta) plus the rod's scattered waves Y_m(k rho) cos(m phi), tied together by
Graf's addition theorem, with the rod's Mie coefficients.  The rod's offset from the axis is fitted per insertion,
its permittivity is shared.  The dielectric filling of each resonance follows from the derivative of the exact
eigenvalue with eps (Hellmann-Feynman), and the copper wall loss of the loaded cavity from the same field.
"""
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import jn_zeros, jv, jvp, yv, yvp

C_MM = 299792458.0e3
ETA = 376.730313668
RADIUS_MM, HEIGHT_MM = 20.0, 10.0
ZEROS = jn_zeros(0, 2)
ORDERS = 9


def resonance(path):
    """Lorentzian magnitude fit: centre frequency (Hz), loaded Q and the peak transmission |S21(f0)|."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    f, mag = data[:, 0] * 1e9, 10 ** (data[:, 1] / 20)
    i = int(np.argmax(mag))
    half = f[mag >= mag[i] / np.sqrt(2)]
    q0 = f[i] / max(half[-1] - half[0], (f[1] - f[0]) * 2)
    scale = f[i]

    def model(p):
        fc, q, peak = p
        return peak / np.sqrt(1 + (2 * q * (f / (fc * scale) - 1)) ** 2)

    fit = least_squares(lambda p: model(p) - mag, [1.0, q0, mag[i]], x_scale=[1 / q0, q0 * 0.05, mag[i] * 0.05])
    fc, q, peak = fit.x
    return fc * scale, q, peak


class Cavity:
    """Even (cos-type) field family of the cavity of radius a with a rod of radius b at distance d."""

    def __init__(self, a, b):
        self.a, self.b = a, b
        self.n = np.arange(-ORDERS, ORDERS + 1)
        self.fold = np.zeros((2 * ORDERS + 1, ORDERS + 1))
        for j in range(ORDERS + 1):
            self.fold[ORDERS + j, j] = 1
            self.fold[ORDERS - j, j] = (-1) ** j     # J_-n = (-1)^n J_n: mirror-even needs alpha_-n = (-1)^n alpha_n

    def mie(self, k, eps):
        n, b, k1 = self.n, self.b, k * np.sqrt(eps)
        top = k * jvp(n, k * b) * jv(n, k1 * b) - k1 * jv(n, k * b) * jvp(n, k1 * b)
        return top / (k1 * yv(n, k * b) * jvp(n, k1 * b) - k * yvp(n, k * b) * jv(n, k1 * b))

    def system(self, k, eps, d):
        n = self.n
        shift = jv(n[None, :] - n[:, None], k * d)                  # shift[m, n] = J_{n-m}(k d)
        t = self.mie(k, eps)
        jn, yn = jv(n, k * self.a), yv(n, k * self.a)
        norm = 1 / np.hypot(np.abs(jn), np.abs(yn))
        full = np.diag(norm * jn) + (norm * yn)[:, None] * (shift.T @ (t[:, None] * shift))
        return full[ORDERS:] @ self.fold, shift, t

    def det(self, k, eps, d):
        return np.linalg.det(self.system(k, eps, d)[0])

    def root(self, guess, eps, d):
        k0, k1 = complex(guess), complex(guess) * (1 - 3e-6)
        f0, f1 = self.det(k0, eps, d), self.det(k1, eps, d)
        for _ in range(50):
            k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
            k0, f0, k1, f1 = k1, f1, k2, self.det(k2, eps, d)
            if abs(k1 - k0) < 1e-14 * abs(k1):
                break
        return k1 if np.iscomplexobj(eps) and np.imag(eps) else k1.real

    def centred_root(self, eps, x_zero):
        """Root of the axially symmetric family (n = 0) for a rod on the axis, scanning down from the empty mode."""
        from scipy.optimize import brentq
        f = lambda k: jv(0, k * self.a) + yv(0, k * self.a) * self.mie(k, eps)[ORDERS]
        top = x_zero / self.a * (1 - 1e-10)
        grid = top * (1 - np.linspace(0, 0.6, 400))
        values = [f(k) for k in grid]
        i = next(i for i in range(len(grid) - 1) if np.sign(values[i]) != np.sign(values[i + 1]))
        return brentq(f, grid[i + 1], grid[i], xtol=1e-15, rtol=1e-15)

    def wall_loss(self, k, eps, d, rs):
        """Q of the copper walls (side wall and both end plates) for the loaded field at real eps."""
        mtx, shift, t = self.system(k, eps, d)
        vh = np.linalg.svd(mtx)[2]
        alpha = self.fold @ vh[-1]
        incident = shift @ alpha
        beta = t * incident
        k1 = k * np.sqrt(eps)
        inner = (incident * jv(self.n, k * self.b) + beta * yv(self.n, k * self.b)) / jv(self.n, k1 * self.b)
        outer = shift.T @ beta
        side = 2 * np.pi * self.a * k * k * np.sum((alpha * jvp(self.n, k * self.a) + outer * yvp(self.n, k * self.a)) ** 2)
        x = k1 * self.b
        rod = 2 * np.pi * np.sum(inner ** 2 * self.b ** 2 / 2 * (jvp(self.n, x) ** 2 + (1 - self.n ** 2 / x ** 2) * jv(self.n, x) ** 2))
        h = 1e-6 * eps
        slope = (math.log(self.root(k, eps + h, d)) - math.log(self.root(k, eps - h, d))) / (2 * h)
        energy = rod / (-2 * slope)                                  # integral of eps |E|^2 over the disk
        filling = -2 * eps * slope                                   # electric energy share in the rod
        q_wall = k ** 3 * ETA * HEIGHT_MM * energy / (rs * (HEIGHT_MM * side + 2 * k * k * energy))
        return q_wall, filling


def table(folder, unloaded=True):
    out = {}
    with open(folder / "sweeps.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            f0, ql, peak = resonance(folder / "sweeps" / row["sweep"])
            # Each probe pair loads the resonator: Q0 = QL / (1 - |S21(f0)|) for two equal probes.
            out[(row["insertion"], int(row["mode"][3]) - 1)] = (f0, ql / (1 - peak) if unloaded else ql)
    return out


def measure(folder, *, model="exact", walls=True, unloaded=True):
    """model: "exact" (the reference); evidence routes "first" (first-order perturbation with the offset from the two
    modes), "hybrid" (exact centred-rod eigenvalues times the first-order offset factor) and "centred" (exact, rod
    taken on the axis)."""
    from scipy.special import j0, j1
    folder = Path(folder)
    b = json.loads((folder / "rod.json").read_text())["rod_diameter_mm"] / 2
    sweeps = table(folder, unloaded)
    labels = sorted({k[0] for k in sweeps if k[0] != "empty"})
    a = [ZEROS[n] * C_MM / (2 * math.pi * sweeps[("empty", n)][0]) for n in range(2)]
    cavity = [Cavity(a[n], b) for n in range(2)]
    k_empty = ZEROS / np.asarray(a)
    observed = np.array([[k_empty[n] * sweeps[(lab, n)][0] / sweeps[("empty", n)][0] for n in range(2)]
                         for lab in labels])
    drop = 1 - observed / k_empty
    fill = (b / RADIUS_MM) ** 2 / j1(ZEROS) ** 2
    count = len(labels)
    low, high = [1.01] + [0] * count, [200] + [9] * count

    # Start: first-order estimate of eps and the offsets.
    best = None
    for start_eps in (2.5, 5.0, 10.0, 20.0, 40.0):
        for start_d in (0.3, 1.5, 2.5, 3.5):
            g = least_squares(lambda p: np.concatenate([(p[0] - 1) / 2 * fill * j0(ZEROS * p[1 + j] / RADIUS_MM) ** 2
                                                         / drop[j] - 1 for j in range(count)]),
                              [start_eps] + [start_d] * count, bounds=(low, high))
            if best is None or g.cost < best.cost:
                best = g
    params = best.x
    if model in ("exact", "centred", "hybrid"):
        def centred(e):
            return np.array([cavity[n].centred_root(e, ZEROS[n]) for n in range(2)])

        def residual(p):
            rows = []
            for j in range(count):
                if model == "exact":
                    k = [cavity[n].root(observed[j, n], p[0], p[1 + j]) for n in range(2)]
                elif model == "centred":
                    k = centred(p[0])
                else:
                    k = k_empty * (1 - (1 - centred(p[0]) / k_empty) * j0(ZEROS * p[1 + j] / RADIUS_MM) ** 2)
                rows.append(np.asarray(k) / observed[j] - 1)
            return np.concatenate(rows)
        if model == "centred":
            params[1:] = 0.0
            high = [200] + [1e-9] * count
        trials = []
        started = time.monotonic()
        for scale in (1.0, 0.93, 0.85):
            if trials and time.monotonic() - started > 15:      # a slow host keeps the best fit so far
                break
            x0 = np.clip(np.r_[1 + (params[0] - 1) * scale, params[1:]], low, high)
            trials.append(least_squares(residual, x0, bounds=(low, high), x_scale=np.r_[0.01 * params[0], [0.05] * count],
                                        diff_step=1e-6, xtol=1e-13, ftol=1e-15))
        params = min(trials, key=lambda r: r.cost).x
    eps = float(params[0])

    tangents = []
    for j, lab in enumerate(labels):
        for n in range(2):
            q_empty = sweeps[("empty", n)][1]
            rs = k_empty[n] * ETA * a[n] * HEIGHT_MM / (2 * q_empty * (HEIGHT_MM + a[n]))
            d = params[1 + j]
            if model == "exact":
                k = cavity[n].root(observed[j, n], eps, d)
                q_wall, filling = cavity[n].wall_loss(k, eps, d, rs * math.sqrt(k / k_empty[n]))
            elif model == "first":
                q_wall, filling = q_empty, eps * fill[n] * j0(ZEROS[n] * d / RADIUS_MM) ** 2
            else:
                k = cavity[n].centred_root(eps, ZEROS[n])
                q_wall, filling = cavity[n].wall_loss(k, eps, 0.0, rs * math.sqrt(k / k_empty[n]))
                filling *= j0(ZEROS[n] * d / RADIUS_MM) ** 2
            if not walls:
                q_wall = q_empty
            tangents.append((1 / sweeps[(lab, n)][1] - 1 / q_wall) / filling)
    return {"eps_real": eps, "loss_tangent": float(np.mean(tangents))}


if __name__ == "__main__":
    answer = measure(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(answer, allow_nan=False) + "\n", encoding="utf-8")
