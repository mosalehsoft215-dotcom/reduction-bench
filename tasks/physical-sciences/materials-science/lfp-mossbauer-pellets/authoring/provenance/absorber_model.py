"""57Fe transmission Mossbauer model used to write the pellet runs.

Velocity v (mm/s) Doppler-shifts the source line. Absorber resonant thickness of one line j (dimensionless):
    T(E) = sum_j t_j (G0/2)(Ga_j/2) / ((E - v_j)^2 + (Ga_j/2)^2),        G0 = 0.097 mm/s natural width (FWHM)
so a line of natural width peaks at t_j. A site of 57Fe areal density n and recoil-free fraction f has total
effective thickness t = sigma0 n f, shared among its lines by their relative strengths.
Counts per channel:  C(v) = B [ 1 - r + r * (L_s conv exp(-T))(v) ],  L_s the source Lorentzian (FWHM Gs),
r the recoilless 14.4 keV fraction of the counting window.
Recoil-free fraction (Debye):  f = exp(-(6 E_R / k theta) [1/4 + (T/theta)^2 int_0^{theta/T} x / (e^x - 1) dx]).
"""
import math

import numpy as np
from scipy.integrate import quad

SIGMA0_CM2 = 2.56e-18
ABUND_57 = 0.02119
G0 = 0.097
E_R_EV = 1.956e-3
K_B_EV = 8.617333262e-5
N_A = 6.02214076e23
M_FE = 55.845
FOIL_LINES = np.array([-5.3123, -3.0760, -0.8397, 0.8397, 3.0760, 5.3123])
FOIL_STRENGTH = np.array([3, 4, 1, 1, 4, 3]) / 16.0
GRID = np.arange(-16.0, 16.0, 0.004)          # energy grid, mm/s


def debye_f(temp_k, theta_k):
    y = theta_k / temp_k
    integral = quad(lambda x: x / math.expm1(x) if x > 0 else 1.0, 0, y)[0]
    return math.exp(-(6 * E_R_EV / (K_B_EV * theta_k)) * (0.25 + integral / y ** 2))


def sod_mm_s(temp_k, theta_k):
    """Second-order Doppler shift of the centre (Debye lattice), mm/s:
    -(9 k theta)/(16 M c) - (9 k T)/(2 M c) (T/theta)^3 int_0^{theta/T} x^3/(e^x - 1) dx,  M = mass of 57Fe."""
    k, m, c = 1.380649e-23, 56.9354 * 1.66053907e-27, 299792458.0
    y = theta_k / temp_k
    integral = quad(lambda x: x ** 3 / math.expm1(x) if x > 0 else 0.0, 0, y)[0]
    return -1e3 * (9 * k * theta_k / (16 * m * c) + 9 * k * temp_k / (2 * m * c) * integral / y ** 3)


def n57_per_cm2(fe_mg_per_cm2):
    return fe_mg_per_cm2 * 1e-3 / M_FE * N_A * ABUND_57


def thickness_profile(lines):
    """lines: iterable of (position, effective thickness, absorber FWHM)."""
    T = np.zeros_like(GRID)
    for pos, t, ga in lines:
        T += t * (G0 / 2) * (ga / 2) / ((GRID - pos) ** 2 + (ga / 2) ** 2)
    return T


def transmission(velocities, lines, gs):
    """(L_s conv exp(-T)) evaluated at the given velocities."""
    absorbed = 1.0 - np.exp(-thickness_profile(lines))
    step = GRID[1] - GRID[0]
    k = np.arange(-len(GRID) // 2, len(GRID) // 2) * step
    src = (gs / 2) / math.pi / (k ** 2 + (gs / 2) ** 2) * step
    conv = np.fft.irfft(np.fft.rfft(absorbed) * np.fft.rfft(np.fft.ifftshift(src)), n=len(GRID))
    return 1.0 - np.interp(velocities, GRID, conv)


def doublet(is_mm_s, qs_mm_s, t_site, ga, weight_low=0.5):
    """Two lines of a quadrupole doublet; weight_low is the share of the lower-velocity line (texture)."""
    return [(is_mm_s - qs_mm_s / 2, t_site * weight_low, ga), (is_mm_s + qs_mm_s / 2, t_site * (1 - weight_low), ga)]


def foil_lines(t_foil, ga=G0):
    return [(p, t_foil * s, ga) for p, s in zip(FOIL_LINES, FOIL_STRENGTH)]
