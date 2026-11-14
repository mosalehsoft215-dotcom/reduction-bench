"""Physical model of one radiotracer serial-sectioning anneal, and the records it leaves.

Diffusion (type-B kinetics of a polycrystal, instantaneous thin-film source of total activity M at the surface):
    c(z) = M / sqrt(pi D t) exp(-z^2 / (4 D t))  +  G exp(-K z^(6/5))
the first term is lattice diffusion, the second the grain-boundary tail, with Le Claire's relation
    s delta D_gb = 1.322 sqrt(D / t) K^(-5/3).
t is the time the furnace thermocouple spent within 5 K of the set point (heating and cooling are fast).

Sectioning: the face (diameter d) is ground off in sections; section i spans z_{i-1}..z_i, with thickness
dm_i / (rho A) from the balance readings before and after it. Its activity at the deposition time is
A_face times the integral of c over the section.

Counting: a well-type NaI detector and a scaler preset on real time. The true rate of a section is
eff_det x activity x 2^(-(t_count - t_deposit) / T_half) + background(t); the scaler loses counts by a
non-paralysable dead time tau, m = n / (1 + n tau); the recorded count is Poisson(m x real time).
Background runs (no sample) and reference-source runs (the tracer isotope, known activity at its own
reference date) are counted on the same scaler and recorded beside the sections.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.special import erf

DAY = 86400.0


def boundary_k(d, t, p):
    """Tail constant K (m^-6/5) from Le Claire's relation."""
    return (1.322 * math.sqrt(d / t) / p) ** 0.6


def _ierf(u):
    """Antiderivative of erf: u erf(u) + exp(-u^2) / sqrt(pi)."""
    return u * erf(u) + math.exp(-u * u) / math.sqrt(math.pi)


def lattice_integral(z0, z1, d, t, m_area, h):
    """Integral over z0..z1 (z0 >= 0) of the lattice term per unit area. h = 0: thin film,
    M / sqrt(pi D t) exp(-z^2 / 4Dt); h > 0: a deposit of thickness h at a reflecting surface (Crank),
    M / (2 h) [erf((h + z) / s) + erf((h - z) / s)], s = 2 sqrt(D t)."""
    s = 2.0 * math.sqrt(d * t)
    if h <= 0.0:
        return m_area * (erf(z1 / s) - erf(z0 / s))
    plus = _ierf((h + z1) / s) - _ierf((h + z0) / s)
    minus = -(_ierf((h - z1) / s) - _ierf((h - z0) / s))
    return m_area / (2.0 * h) * s * (plus + minus)


def section_activity(z0, z1, d, t, m_area, g, k, area, h=0.0, wedge=0.0, holdup=0.0, holdup_depth=1.5e-7):
    """Activity of the slice z0..z1 (Bq). A disc ground at a tilt (height difference wedge across the circular
    face) cuts a plane whose depth over the face follows the chord law: the share of the face at height u from the
    mean goes as sqrt(1 - (2 u / wedge)^2). holdup is the share of the deposit held in a surface film
    0..holdup_depth that never diffused."""
    offsets = [0.0] if wedge <= 0 else list(tilt_offsets(wedge))
    total = 0.0
    for u in offsets:
        a, b = max(z0 + u, 0.0), max(z1 + u, 0.0)
        if b <= a:
            continue
        part = (1.0 - holdup) * lattice_integral(a, b, d, t, m_area, h)
        part += quad(lambda z: g * math.exp(-k * z ** 1.2), a, b, limit=200)[0]
        if holdup:
            part += holdup * m_area * max(0.0, min(b, holdup_depth) - a) / holdup_depth
        total += part
    return area * total / len(offsets)


def tilt_offsets(wedge, n=21):
    """n equal-share heights of a tilted circular face (quantile midpoints of the chord law)."""
    u = np.linspace(-1.0, 1.0, 4001)
    cdf = np.cumsum(np.sqrt(1.0 - u * u))
    cdf /= cdf[-1]
    return 0.5 * wedge * np.interp((np.arange(n) + 0.5) / n, cdf, u)


def effective_time(times_s, temps_c, setpoint_c):
    """Time within 5 K of the set point, from a furnace record (trapezoid on the in-band intervals)."""
    t = np.asarray(times_s, float)
    inband = np.abs(np.asarray(temps_c, float) - setpoint_c) <= 5.0
    both = inband[1:] & inband[:-1]
    return float(np.sum(np.diff(t)[both]))
