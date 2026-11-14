"""Thin-film resonance on a coplanar waveguide, field swept at a fixed microwave frequency (model for the generator).

Easy-plane film (effective magnetisation Meff, tesla), static field in the x-z plane at theta_B from the normal z.
Equilibrium polar angle theta of m:   B sin(theta - theta_B) = (Meff / 2) sin(2 theta)
Smit-Beljers stiffness fields:         B1 = B cos(theta - theta_B) - Meff cos(2 theta)
                                       B2 = B cos(theta - theta_B) - Meff cos^2(theta)
Resonance:                             (omega / gamma)^2 = B1 B2,  gamma = g mu_B / hbar
Gilbert half width in field:           alpha (omega / gamma) (B1 + B2) / |d(B1 B2)/dB|
Inhomogeneous half width:              |dB_res / dMeff| dM   (Lorentzian spread dM of Meff)
Two-magnon half width (in plane only): G2 asin sqrt((x - f0/2)/(x + f0/2)),  x = sqrt(f^2 + f0^2/4),  f0 = gamma Meff / 2 pi
Line: complex Lorentzian in field, 1 / (B_res - B - i HW); the field modulation of peak amplitude a turns it into the
first harmonic (2/a)(u / sqrt(u^2 - a^2) - 1), u = B_res - B - i HW (-> a / u^2 for small a); the detector reads
Re(e^{i phi} . ) on a sloping background, then the output filter (n first-order sections of time constant tau) acts
along the sweep in time.  Gaussmeter: reading = s(T) B_true + o,  s(T) = s25 (1 + beta (T - 25 C)).
"""
import math

import numpy as np
from scipy.signal import lfilter

MU_B = 9.2740100783e-24
HBAR = 1.054571817e-34
G_MARKER = 2.0036


def gamma_of(g):
    return g * MU_B / HBAR


def w_field(freq, g):
    """omega / gamma in tesla."""
    return 2 * math.pi * freq / gamma_of(g)


def angle_m(B, th_b, meff):
    B = np.asarray(B, float)
    lo = np.full(B.shape, float(th_b))
    hi = np.full(B.shape, math.pi / 2)
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        neg = B * np.sin(mid - th_b) - 0.5 * meff * np.sin(2 * mid) < 0
        lo = np.where(neg, mid, lo)
        hi = np.where(neg, hi, mid)
    return 0.5 * (lo + hi)


def stiff(B, th_b, meff):
    th = angle_m(B, th_b, meff)
    c = np.cos(th - th_b)
    return B * c - meff * np.cos(2 * th), B * c - meff * np.cos(th) ** 2


def res_field(freq, g, meff, th_b, b_max=3.0):
    """Upper (saturated-branch) resonance field in tesla."""
    w2 = w_field(freq, g) ** 2
    grid = np.linspace(1e-4, b_max, 3000)
    b1, b2 = stiff(grid, th_b, meff)
    k = int(np.argmin(b1 * b2))
    lo, hi = grid[k], b_max
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        b1, b2 = stiff(np.array([mid]), th_b, meff)
        if b1[0] * b2[0] < w2:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def widths(freq, film, th_b, in_plane):
    """(B_res, Gilbert part, inhomogeneous part, two-magnon part), tesla."""
    g, meff = film["g"], film["meff"]
    br = res_field(freq, g, meff, th_b)
    h = 1e-6
    b1, b2 = stiff(np.array([br - h, br, br + h]), th_b, meff)
    F = b1 * b2
    dF = (F[2] - F[0]) / (2 * h)
    hom = film["alpha"] * w_field(freq, g) * (b1[1] + b2[1]) / abs(dF)
    dm = 1e-5
    inh = abs(res_field(freq, g, meff + dm, th_b) - res_field(freq, g, meff - dm, th_b)) / (2 * dm) * film["dM"]
    tm = 0.0
    if in_plane and film.get("G2", 0) > 0:
        f0 = gamma_of(g) * meff / (2 * math.pi)
        x = math.sqrt(freq ** 2 + f0 ** 2 / 4)
        tm = film["G2"] * math.asin(math.sqrt((x - f0 / 2) / (x + f0 / 2)))
    return br, hom, inh, tm


def first_harmonic(B, br, hw, a):
    u = br - B - 1j * hw
    if a <= 0:
        return 1.0 / u ** 2
    s = np.sqrt(u - a) * np.sqrt(u + a)
    return (2.0 / a) * (u / s - 1.0)


def output_filter(y, step_s, tau_s, poles):
    """n first-order sections along acquisition order, settled on the first value."""
    if tau_s <= 0 or poles <= 0:
        return np.asarray(y, float)
    p = math.exp(-step_s / tau_s)
    out = np.asarray(y, float)
    for _ in range(poles):
        out = lfilter([1 - p], [1, -p], out, zi=[p * out[0]])[0]
    return out
