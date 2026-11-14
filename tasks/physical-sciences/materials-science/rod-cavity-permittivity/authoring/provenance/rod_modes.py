"""Exact TM0n0-family resonances of a cylindrical cavity holding a full-height dielectric rod off its axis.

The rod runs through the whole height, so the fields do not vary with z and E_z obeys the 2-D Helmholtz equation
with k^2 eps(x) in a disk of radius a, E_z = 0 on the wall.  Outside the rod
    E = sum_n alpha_n J_n(k r) e^{i n theta} + sum_m beta_m Y_m(k rho) e^{i m phi}
(cavity coordinates r, theta; rod coordinates rho, phi about the rod centre at distance d on the x axis), inside
    E = sum_m gamma_m J_m(k sqrt(eps) rho) e^{i m phi}.
Graf's addition theorem moves the cavity waves to the rod (J_n e^{in theta} = sum_m J_{n-m}(kd) J_m e^{im phi}) and the
rod's waves to the wall (Y_m e^{im phi} = sum_n J_{n-m}(kd) Y_n(kr) e^{in theta}, r > d); the rod's Mie coefficient
gives beta_m = T_m a_m, and the wall condition leaves M(k) alpha = 0.

The stored-energy integral S = int eps |E|^2 dA comes from the Rellich identity (boundary integrals only):
    S = a/(2k^2) * oint_wall (dE/dr)^2 ds + (eps - 1)/2 * oint_rod (d cos phi + R) E^2 R dphi,
and the wall loss of the loaded cavity from P ~ h oint_wall (dE/dr)^2 ds + 2 k^2 S (side wall plus both end plates),
so that Q_c = k^3 eta h S / (Rs (h W + 2 k^2 S)).
"""
import numpy as np
from scipy.special import jv, jvp, yv, yvp

ORDERS = 12


def mie(k, eps, radius, m):
    k1 = k * np.sqrt(eps)
    num = k * jvp(m, k * radius) * jv(m, k1 * radius) - k1 * jv(m, k * radius) * jvp(m, k1 * radius)
    den = k1 * yv(m, k * radius) * jvp(m, k1 * radius) - k * yvp(m, k * radius) * jv(m, k1 * radius)
    return num / den


def matrix(k, eps, radius, d, a, orders=ORDERS):
    n = np.arange(-orders, orders + 1)
    g = jv(n[None, :] - n[:, None], k * d)                     # g[m, n] = J_{n-m}(k d)
    t = mie(k, eps, radius, n)
    jn, yn = jv(n, k * a), yv(n, k * a)
    scale = 1 / np.sqrt(np.abs(jn) ** 2 + np.abs(yn) ** 2)       # row scale without zeros or poles
    return (scale * jn)[:, None] * np.eye(len(n)) + (scale * yn)[:, None] * (g.T @ (t[:, None] * g)), g, t, n


def even(mtx, orders=ORDERS):
    """The rod sits on the x axis, so cos-type fields (alpha_n = alpha_-n) form their own family."""
    fold = np.zeros((2 * orders + 1, orders + 1))
    for j in range(orders + 1):
        fold[orders + j, j] = 1
        fold[orders - j, j] = (-1) ** j          # J_-n = (-1)^n J_n: mirror-even fields need alpha_-n = (-1)^n alpha_n
    return mtx[orders:] @ fold, fold


def det(k, eps, radius, d, a):
    return np.linalg.det(even(matrix(k, eps, radius, d, a)[0])[0])


def secant(k_guess, eps, radius, d, a, step=-2e-6):
    k0, k1 = complex(k_guess), complex(k_guess * (1 + step))
    f0, f1 = det(k0, eps, radius, d, a), det(k1, eps, radius, d, a)
    for _ in range(60):
        k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
        k0, f0, k1, f1 = k1, f1, k2, det(k2, eps, radius, d, a)
        if abs(k1 - k0) < 1e-14 * abs(k1):
            break
    return k1


def follow(x_zero, eps, radius, d, a, reach=(0.2, 0.12)):
    """Real root of the cos-type mode that is TM0n0 in the empty cavity: the only sign change of the determinant
    below the empty value (scanning down from it): rods only lower the frequencies, and the first root met is the
    one connected to the empty TM0n0 mode even where it mixes with TM210."""
    from scipy.optimize import brentq
    top = x_zero / a * (1 - 1e-9)
    grid = top * (1 - np.linspace(0, reach[0] if x_zero < 3 else reach[1], 241))
    values = [det(k, eps.real, radius, d, a) for k in grid]
    for i in range(len(grid) - 1):
        if np.sign(values[i]) != np.sign(values[i + 1]):
            root = brentq(lambda k: det(k, eps.real, radius, d, a), grid[i + 1], grid[i], xtol=1e-15, rtol=1e-15)
            return root if np.isreal(eps) else secant(root, eps, radius, d, a)
    raise ValueError("no TM0n0 root in range")


def coefficients(k, eps, radius, d, a):
    mtx, g, t, n = matrix(k, eps, radius, d, a)
    reduced, fold = even(mtx)
    u, s, vh = np.linalg.svd(reduced)
    alpha = fold @ vh[-1].conj()
    incident = g @ alpha                                         # a_m on the rod
    beta = t * incident
    k1 = k * np.sqrt(eps)
    gamma = (incident * jv(n, k * radius) + beta * yv(n, k * radius)) / jv(n, k1 * radius)
    wall = g.T @ beta                                            # rod waves expressed on the wall
    return n, alpha, beta, gamma, wall


def wall_q(k, eps, radius, d, a, height, rs, eta=376.730313668):
    """Wall-loss Q of the loaded cavity (real eps), with Rs the copper surface resistance at this frequency."""
    n, alpha, beta, gamma, wall = coefficients(k, eps, radius, d, a)
    theta = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    ephase = np.exp(1j * np.outer(theta, n))
    radial = ephase @ (k * (alpha * jvp(n, k * a) + wall * yvp(n, k * a)))
    w_side = np.mean(np.abs(radial) ** 2) * 2 * np.pi * a
    on_rod = ephase @ (gamma * jv(n, k * np.sqrt(eps) * radius))
    interface = np.mean((d * np.cos(theta) + radius) * np.abs(on_rod) ** 2) * 2 * np.pi * radius
    s = a / (2 * k * k) * w_side + (eps - 1) / 2 * interface
    return k ** 3 * eta * height * s / (rs * (height * w_side + 2 * k * k * s)), s, w_side
