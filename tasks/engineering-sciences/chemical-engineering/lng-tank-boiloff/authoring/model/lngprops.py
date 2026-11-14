"""LNG property model used by the tank simulation and stated in the record guide.

Liquid: ideal solution at its bubble point. Vapour: ideal gas, equilibrium composition by Raoult's law.
Units: T in K, pressures in kPa, molar quantities per kmol (kJ/kmol, kJ/(kmol K), m3/kmol).
"""
import numpy as np

COMP = ("N2", "CH4", "C2H6", "C3H8", "nC4H10")
M = np.array([28.0134, 16.0425, 30.069, 44.0956, 58.1222])
# Antoine, log10(p / bar) = A - B / (T + C)
ANT_A = np.array([3.7362, 3.9895, 4.50706, 4.53678, 4.35576])
ANT_B = np.array([264.651, 443.028, 791.3, 1149.36, 1175.581])
ANT_C = np.array([-6.788, -0.49, -6.422, 24.906, -2.071])
# Watson: dHvap(T) = dH_nb ((Tc - T) / (Tc - Tnb))^0.38
DH_NB = np.array([5570.0, 8190.0, 14690.0, 19040.0, 22440.0])  # kJ/kmol
T_NB = np.array([77.35, 111.67, 184.57, 231.04, 272.66])
T_C = np.array([126.19, 190.56, 305.32, 369.83, 425.12])
CP_LIQ = np.array([60.0, 54.8, 68.2, 84.6, 114.0])  # kJ/(kmol K)
V_LIQ = np.array([0.0460, 0.0379, 0.0465, 0.0610, 0.0775])  # m3/kmol
R = 8.314462  # kPa m3 / (kmol K)


def psat(T):
    """kPa, one value per component (T scalar)."""
    return 100.0 * 10.0 ** (ANT_A - ANT_B / (T + ANT_C))


def dpsat_dT(T):
    p = psat(T)
    return p * np.log(10.0) * ANT_B / (T + ANT_C) ** 2


def bubble_T(P, x, T0=112.0):
    T = T0
    for _ in range(40):
        f = np.dot(x, psat(T)) - P
        T -= f / np.dot(x, dpsat_dT(T))
        if abs(f) < 1e-10 * P:
            break
    return T


def latent(T):
    return DH_NB * np.clip((T_C - T) / (T_C - T_NB), 0.0, None) ** 0.38


def vapour_of(P, x, T):
    return x * psat(T) / P


def bubble_derivs(P, x, T):
    """dT/dP and dT/dx_i (each x_i taken as free) at the bubble point."""
    s = np.dot(x, dpsat_dT(T))
    return 1.0 / s, -psat(T) / s
