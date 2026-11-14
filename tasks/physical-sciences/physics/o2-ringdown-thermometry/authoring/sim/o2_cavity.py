"""Physics of the simulated O2 A-band pulsed cavity ring-down campaigns.

Molecule: 16O2, b1Sigma_g+ (v'=0) <- X3Sigma_g- (v''=0), magnetic-dipole A band near 763 nm.
  * ground-state rotational levels: Hund's case (b) triplet formulas with B, D, lambda, gamma
    (odd N'' only for 16O2), energies relative to the lowest level N''=1, J''=0;
  * upper state: 1Sigma, F(J) = B J (J+1) - D J^2 (J+1)^2;
  * line intensities at 296 K from a Boltzmann population, a simple branch strength and the
    stimulated-emission factor, scaled to a band maximum of 9.4e-24 cm/molecule (HITRAN order);
  * partition sum: direct sum over the ground-state levels times the harmonic vibrational factor.
Absorption: Voigt lines (Doppler from the gas temperature, air broadening with a temperature exponent),
number density from the ideal gas.
Cavity: ring-down rate k(nu) = c (1 - R(nu)) / L + c * fill * alpha(nu). A pulse of the laser excites all
frequencies of its spectral profile at once; the detector sees the sum of their exponentials, each with its
own rate, through a single-pole response with time constant tau_d.
"""
import numpy as np
from scipy.special import voigt_profile

C_CM = 2.99792458e10                 # speed of light, cm/s
C2 = 1.438776877                     # second radiation constant, cm K
KB = 1.380649e-23
AMU = 1.66053907e-27
M_O2 = 31.98983                      # 16O2, u
TORR_TO_CM3 = 133.322368 / KB * 1e-6  # number density per Torr per (1/K), cm^-3

B0, D0, LAMBDA, GAMMA = 1.4376766, 4.8418e-6, 1.984751, -0.008446
B1, D1, NU_BAND = 1.3912, 5.37e-6, 13122.0057
OMEGA_E = 1556.4                     # vibrational wavenumber of the ground state, cm-1
S_MAX = 9.4e-24


def f_lower(n, j):
    """X3Sigma_g- term value (cm-1) of rotational level N, J"""
    rot = B0 * n * (n + 1) - D0 * (n * (n + 1)) ** 2
    if j == n:
        return rot
    if j == n + 1:
        return rot + (2 * n + 3) * B0 - LAMBDA - np.sqrt((2 * n + 3) ** 2 * B0 ** 2 + LAMBDA ** 2 - 2 * LAMBDA * B0) \
            + GAMMA * (n + 1)
    return rot - (2 * n - 1) * B0 - LAMBDA + np.sqrt((2 * n - 1) ** 2 * B0 ** 2 + LAMBDA ** 2 - 2 * LAMBDA * B0) \
        - GAMMA * n


def f_upper(j):
    return B1 * j * (j + 1) - D1 * (j * (j + 1)) ** 2


E_ZERO = min(f_lower(1, j) for j in (0, 1, 2))


def ground_levels(n_max=161):
    out = []
    for n in range(1, n_max, 2):
        for j in (n - 1, n, n + 1):
            if j >= 0:
                out.append((f_lower(n, j) - E_ZERO, 2 * j + 1))
    return np.array(out)


_LEV = ground_levels()


def partition_sum(t):
    t = np.asarray(t, float)
    rot = np.sum(_LEV[None, :, 1] * np.exp(-C2 * _LEV[None, :, 0] / t.reshape(-1, 1)), axis=1)
    return (rot / (1 - np.exp(-C2 * OMEGA_E / t.reshape(-1)))).reshape(t.shape)


def band_lines(n_max=35):
    """rows: branch, N'', J'', wavenumber (cm-1), S(296 K) (cm/molecule), E'' (cm-1), gamma_air (cm-1/atm), n_air"""
    rows = []
    q296 = float(partition_sum(296.0))
    for n in range(1, n_max + 1, 2):
        for branch, dn, jpp in (("RR", 1, n), ("RQ", 1, n + 1), ("PP", -1, n), ("PQ", -1, n - 1)):
            jp = n + dn
            if jp < 0 or jpp < 0 or (branch == "PQ" and n == 1):
                continue
            e = f_lower(n, jpp) - E_ZERO
            nu = NU_BAND + f_upper(jp) - (f_lower(n, jpp) - E_ZERO) - E_ZERO
            strength = {"RR": (n + 2) / 2, "RQ": (n + 1) / 2, "PP": (n - 1 + 1e-9) / 2, "PQ": n / 2}[branch]
            s = nu * strength * np.exp(-C2 * e / 296.0) * (1 - np.exp(-C2 * nu / 296.0)) / q296
            g_air = 0.0482 - 0.00041 * n + 0.0000035 * n * n
            rows.append([branch, n, jpp, nu, s, e, g_air, 0.71 - 0.002 * n])
    smax = max(r[4] for r in rows)
    for r in rows:
        r[4] *= S_MAX / smax
    return rows


def intensity_at(s296, e_low, nu, t):
    """HITRAN-convention temperature conversion of line intensities"""
    return s296 * partition_sum(296.0) / partition_sum(t) * np.exp(-C2 * e_low * (1 / t - 1 / 296.0)) \
        * (1 - np.exp(-C2 * nu / t)) / (1 - np.exp(-C2 * nu / 296.0))


def absorption(nu, lines, t, p_torr, x_o2):
    """alpha(nu) in cm-1 of the gas (lines: array of nu, S296, E'', gamma_air, n_air)"""
    nu = np.asarray(nu, float)
    dens = TORR_TO_CM3 * p_torr / t * x_o2
    s_t = intensity_at(lines[:, 1], lines[:, 2], lines[:, 0], t)
    sig_d = lines[:, 0] * np.sqrt(KB * t / (M_O2 * AMU)) / 2.99792458e8
    hw_l = lines[:, 3] * (p_torr / 760.0) * (296.0 / t) ** lines[:, 4]
    out = np.zeros_like(nu)
    near = np.abs(nu.reshape(-1, 1) - lines[:, 0]).min(axis=0) < 2.0
    for i in np.flatnonzero(near):
        out += s_t[i] * voigt_profile(nu - lines[i, 0], sig_d[i], hw_l[i])
    return dens * out


def laser_offsets(shape, fwhm, n=321):
    """spectral profile of the laser: offsets (cm-1) and normalised weights on a uniform grid"""
    if shape == "gaussian":
        sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
        d = np.linspace(-6 * sig, 6 * sig, n)
        w = np.exp(-0.5 * (d / sig) ** 2)
    elif shape == "sech2":
        a = fwhm / (2 * np.arccosh(np.sqrt(2)))
        d = np.linspace(-9 * a, 9 * a, n)
        w = 1 / np.cosh(d / a) ** 2
    else:
        raise ValueError(shape)
    return d, w / w.sum()


def detected_decay(t, rates, weights, tau_d):
    """sum of exponentials seen through a single-pole detector, for an impulse at t = 0"""
    r = rates[:, None]
    return np.sum(weights[:, None] * (np.exp(-r * t[None, :]) - np.exp(-t[None, :] / tau_d)) / (1 - r * tau_d), axis=0)
