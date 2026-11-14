"""Linear optics of the quadrupole-scan line, shared by the generator and the reference.

The scanned quadrupole (effective length LQ, gradient G) is followed by a drift DRIFT to the screen.
An electron with relative momentum deviation delta sees k = 0.299792458 G[T/m] / p0[GeV/c] / (1 + delta),
focusing in the scan plane. At the screen its rms size is
    sigma^2 = eps (M11^2 beta - 2 M11 M12 alpha + M12^2 (1 + alpha^2) / beta)
for the Twiss parameters (beta, alpha) and geometric emittance eps at the quadrupole entrance,
which are taken to be the same for every momentum slice.
"""
import numpy as np

ME = 0.51099895                       # electron rest energy, MeV
C_GEV = 0.299792458                   # k = C_GEV * G[T/m] / p[GeV/c]
LQ, DRIFT = 0.200, 1.500              # m


def momentum_gev(t_mev):
    return np.sqrt(t_mev ** 2 + 2 * t_mev * ME) / 1000.0


def gamma_rel(t_mev):
    return 1.0 + t_mev / ME


def m11_m12(k):
    """thick focusing quadrupole (k > 0) followed by the drift; returns (M11, M12)"""
    k = np.asarray(k, float)
    sk = np.sqrt(k)
    c, s, sp = np.cos(sk * LQ), np.sin(sk * LQ) / sk, -sk * np.sin(sk * LQ)
    return c + DRIFT * sp, s + DRIFT * c


def size2(k, eps, beta, alpha):
    m11, m12 = m11_m12(k)
    return eps * (m11 ** 2 * beta - 2 * m11 * m12 * alpha + m12 ** 2 * (1 + alpha ** 2) / beta)


def hermite(n=31):
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


def profile(xc, k0, eps, beta, alpha, spread, res, n=31):
    """expected line-out shape (unit area, per metre) at positions xc for nominal k0: a mixture over a
    Gaussian momentum distribution of rms `spread` of Gaussians blurred by the screen resolution"""
    xq, wq = hermite(n)
    out = np.zeros_like(np.asarray(xc, float))
    for d, w in zip(xq * spread, wq):
        s2 = size2(k0 / (1 + d), eps, beta, alpha) + res ** 2
        out += w * np.exp(-xc ** 2 / (2 * s2)) / np.sqrt(2 * np.pi * s2)
    return out
