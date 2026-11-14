"""Uniform-field drift tube for ion mobility, with the physics the acquisitions are drawn from.

    K_low = (3 z e / (16 N)) * sqrt(2 pi / (mu k T)) / Omega            Mason-Schamp, N = P / (k T)
    1 / K_mix = sum_g x_g / K_g                                          Blanc, for a gas mixture
    v_d = K_low E / sqrt(1 + M_gas v_d^2 / (3 k T))                     Wannier heating, Omega held constant
    t_A = t0 + L / v_d                                                   arrival after the trap gate opens
    sigma_t^2 = gate^2 + 2 D (L / v_d) / v_d^2,  D = K k T / (z e)      longitudinal diffusion
Pressures in Torr, temperatures in K, lengths in cm, times in s, cross sections in A^2.
"""
import math

import numpy as np

E = 1.602176634e-19
KB = 1.380649e-23
AMU = 1.66053906660e-27
TORR = 133.322368
GAS_MASS = {"He": 4.002602, "N2": 28.0134}


def number_density(p_torr, t_k):
    return p_torr * TORR / (KB * t_k)                     # m^-3


def mobility_low(omega_a2, z, m_ion, gas, p_torr, t_k):
    """Low-field mobility (m^2 / V s) in a pure gas."""
    mu = m_ion * GAS_MASS[gas] / (m_ion + GAS_MASS[gas]) * AMU
    n = number_density(p_torr, t_k)
    return 3 * z * E / (16 * n) * math.sqrt(2 * math.pi / (mu * KB * t_k)) / (omega_a2 * 1e-20)


def mobility_mix(ion, mix, p_torr, t_k):
    """Low-field mobility in a mixture {gas: mole fraction}, by Blanc's law."""
    inv = sum(x / mobility_low(ion["ccs"][g], ion["z"], ion["mass"], g, p_torr, t_k) for g, x in mix.items() if x > 0)
    return 1.0 / inv


def drift_velocity(k_low, field_v_m, mix, t_k):
    """Wannier heating with the mean gas mass of the mixture."""
    m_gas = sum(x * GAS_MASS[g] for g, x in mix.items()) * AMU
    v = k_low * field_v_m
    for _ in range(60):
        v = k_low * field_v_m / math.sqrt(1.0 + m_gas * v * v / (3 * KB * t_k))
    return v


def arrival(ion, mix, p_torr, t_k, volts, length_cm, t0_s, gate_s):
    """Centre and width (s) of the arrival-time distribution."""
    L = length_cm / 100.0
    k_low = mobility_mix(ion, mix, p_torr, t_k)
    v = drift_velocity(k_low, volts / L, mix, t_k)
    td = L / v
    k_eff = v / (volts / L)
    D = k_eff * KB * t_k / (ion["z"] * E)
    sigma = math.sqrt(gate_s ** 2 + 2 * D * td / v ** 2)
    return t0_s + td, sigma


def atd(rng, centre, sigma, counts, bin_s=20e-6, half_width=6.0, baseline=3.0):
    lo = centre - half_width * sigma
    edges = np.arange(np.floor(lo / bin_s) * bin_s, centre + half_width * sigma + bin_s, bin_s)
    mids = 0.5 * (edges[1:] + edges[:-1])
    from scipy.special import erf
    cdf = 0.5 * (1 + erf((edges - centre) / (sigma * math.sqrt(2))))
    expected = counts * np.diff(cdf) + baseline
    return mids, rng.poisson(expected)
