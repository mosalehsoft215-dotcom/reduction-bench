"""Physics behind every record of the nanofluid-carboy bench (mode-mismatched dual-beam thermal lens).

Probe centre intensity after the excitation light reaches the cell at t0 (Shen, Lowe & Snook 1992):
    I(t)/I0 = [1 - (theta/2) atan(2 m V / (((1 + 2m)^2 + V^2) tc / (2 (t - t0)) + 1 + 2m + V^2))]^2
    theta = - P_abs (dn/dT) / (k lambda_p),   tc = w_e^2 / (4 D),   m = (w_p(z1) / w_e)^2,   V = z1 / zc
    P_abs = eta P_meter (1 - exp(-alpha l))       (Beer-Lambert along the cell; nothing fluoresces)
k and D belong to the liquid in the cell. For a reference liquid they are those of its base fluid. For a nanofluid
carboy dn/dT is that of the base fluid, while its own conductivity is k_carboy = rho c_p D_carboy with rho and c_p from the
carboy sheet; the carboy's conductivity is k_carboy = gain * k_base(T), so D_carboy = gain * k_base(T) / (rho c_p).
"""
import math

import numpy as np

PROBE_NM, PUMP_NM = 632.8, 532.0

# base fluids against temperature in deg C: value at 20 C and slope per K
BASE = {
    "water":           {"k": (0.5980, 1.6e-3),  "dndt": (-0.91e-4, -3.9e-6), "D": (1.430e-7, 4.6e-10)},
    "ethanol":         {"k": (0.1690, -2.3e-4), "dndt": (-4.01e-4, -1.6e-6), "D": (8.92e-8, -2.4e-10)},
    "ethylene glycol": {"k": (0.2520, 1.0e-4),  "dndt": (-2.60e-4, -1.0e-6), "D": (9.36e-8, 1.0e-10)},
}


def base(fluid, temp_c):
    return {q: a + b * (temp_c - 20.0) for q, (a, b) in BASE[fluid].items()}


def mode_terms(w0, z_waist, z_cell, w_e):
    zc = math.pi * w0 * w0 / (PROBE_NM * 1e-9)
    z1 = z_cell - z_waist
    return (w0 * w0 * (1.0 + (z1 / zc) ** 2)) / (w_e * w_e), z1 / zc


def probe_ratio(t, t0, theta, tc, m, V):
    out = np.ones_like(t)
    lit = t > t0
    dt = t[lit] - t0
    q = 2 * m * V / (((1 + 2 * m) ** 2 + V * V) * tc / (2 * dt) + 1 + 2 * m + V * V)
    out[lit] = (1.0 - 0.5 * theta * np.arctan(q)) ** 2
    return out


def lens_strength(p_meter_w, eta, alpha_cm, path_cm, k, dndt):
    p_abs = eta * p_meter_w * (1.0 - math.exp(-alpha_cm * path_cm))
    return -p_abs * dndt / (k * PROBE_NM * 1e-9)
