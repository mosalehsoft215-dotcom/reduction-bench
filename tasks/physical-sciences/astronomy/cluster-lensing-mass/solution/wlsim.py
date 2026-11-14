"""Weak-lensing model and synthetic catalogues for a stack of identical NFW clusters.

Cosmology: flat LCDM, H0 = 70, Om = 0.3. Distances in Mpc (physical), masses in Msun.
NFW shear: Wright & Brainerd (2000). Source ellipticity: e = (e_int + g) / (1 + g* e_int).
"""
import math
import numpy as np

C_KMS = 299792.458
G_MPC = 4.300917270e-9           # Mpc (km/s)^2 / Msun
H0, OM = 70.0, 0.3


def E(z):
    return np.sqrt(OM * (1 + z) ** 3 + 1 - OM)


_zg = np.linspace(0, 6, 6001)
_dc = np.concatenate([[0], np.cumsum(0.5 * (1 / E(_zg[1:]) + 1 / E(_zg[:-1])) * np.diff(_zg))]) * C_KMS / H0


def Dc(z):
    return np.interp(z, _zg, _dc)


def Da(z):
    return Dc(z) / (1 + z)


def Da12(z1, z2):
    return np.maximum(Dc(z2) - Dc(z1), 0) / (1 + z2)


def sigma_crit_inv(zl, zs):
    """1 / Sigma_crit in Mpc^2 / Msun (physical); 0 for zs <= zl."""
    return np.where(zs > zl, 4 * math.pi * G_MPC / C_KMS ** 2 * Da(zl) * Da12(zl, zs) / np.maximum(Da(zs), 1e-9), 0.0)


def rho_crit(z):
    H = H0 * E(z)
    return 3 * H * H / (8 * math.pi * G_MPC)          # Msun / Mpc^3


def nfw_params(M200, c, z):
    r200 = (3 * M200 / (4 * math.pi * 200 * rho_crit(z))) ** (1 / 3)
    rs = r200 / c
    dc = 200 / 3 * c ** 3 / (math.log(1 + c) - c / (1 + c))
    return rs, dc * rho_crit(z)


def _g_less(x):
    return (8 * np.arctanh(np.sqrt((1 - x) / (1 + x))) / (x * x * np.sqrt(1 - x * x)) + 4 / (x * x) * np.log(x / 2)
            - 2 / (x * x - 1) + 4 * np.arctanh(np.sqrt((1 - x) / (1 + x))) / ((x * x - 1) * np.sqrt(1 - x * x)))


def _g_more(x):
    return (8 * np.arctan(np.sqrt((x - 1) / (1 + x))) / (x * x * np.sqrt(x * x - 1)) + 4 / (x * x) * np.log(x / 2)
            - 2 / (x * x - 1) + 4 * np.arctan(np.sqrt((x - 1) / (1 + x))) / ((x * x - 1) ** 1.5))


def _f_sigma(x):
    out = np.empty_like(x)
    lo, hi = x < 1 - 1e-6, x > 1 + 1e-6
    out[lo] = 2 / (x[lo] ** 2 - 1) * (1 - 2 / np.sqrt(1 - x[lo] ** 2) * np.arctanh(np.sqrt((1 - x[lo]) / (1 + x[lo]))))
    out[hi] = 2 / (x[hi] ** 2 - 1) * (1 - 2 / np.sqrt(x[hi] ** 2 - 1) * np.arctan(np.sqrt((x[hi] - 1) / (1 + x[hi]))))
    mid = ~(lo | hi)
    out[mid] = 2 / 3
    return out


def nfw_sigma_dsigma(R, M200, c, z):
    """Surface density Sigma(R) and Delta Sigma(R), Msun / Mpc^2."""
    rs, rhos = nfw_params(M200, c, z)
    x = np.asarray(R, float) / rs
    sig = rs * rhos * _f_sigma(x)
    g = np.empty_like(x)
    lo, hi = x < 1 - 1e-6, x > 1 + 1e-6
    g[lo] = _g_less(x[lo]); g[hi] = _g_more(x[hi])
    mid = ~(lo | hi)
    g[mid] = 10 / 3 + 4 * math.log(0.5)
    return sig, rs * rhos * g
