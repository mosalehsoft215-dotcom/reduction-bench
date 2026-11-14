"""Gas thermodynamics and stoichiometry of a sulfur-burning contact plant.

Molar flows in kmol/h, temperatures in K, enthalpies in kJ/mol. Dry air is O2 20.946, N2 78.084, Ar 0.934 and
CO2 0.036 mol-%. The furnace burns sulfur to SO2 and a fraction f of it to SO3; the converter oxidises SO2 with
half a mole of O2; the absorbers take out all SO3.
"""
import numpy as np

AIR = {"O2": 0.20946, "N2": 0.78084, "Ar": 0.00934, "CO2": 0.00036}
M_S, M_SO2, M_H2SO4, M_AIR = 32.06, 64.066, 98.079, 28.9647
V_N = 22.414  # m3 per kmol at 0 C and 101.325 kPa

# NIST Shomate coefficients A..F; H_abs(T) = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F (kJ/mol), t = T/1000
SHOMATE = {
    "N2": [(500.0, (28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914)),
           (2000.0, (19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202))],
    "O2": [(700.0, (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471)),
           (2000.0, (30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468))],
    "Ar": [(6000.0, (20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.197350))],
    "CO2": [(1200.0, (24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075))],
    "SO2": [(1200.0, (21.43049, 74.35094, -57.75217, 16.35534, 0.086731, -305.7688))],
    "SO3": [(1200.0, (24.02503, 119.4607, -94.38686, 26.96237, -0.117517, -407.8526))],
}
SPECIES = ("N2", "O2", "Ar", "CO2", "SO2", "SO3")


def h_abs(sp, T):
    T = np.asarray(T, dtype=float)
    out = np.zeros_like(T)
    lo = 0.0
    for hi, (A, B, C, D, E, F) in SHOMATE[sp]:
        m = (T >= lo) & (T < hi) if hi < 6000.0 else (T >= lo)
        t = T[m] / 1000.0
        out[m] = A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F
        lo = hi
    return out


def stream(F_air, n_s, x_so3, absorbed=0.0):
    """Molar flows of each species once a share x_so3 of the sulfur is SO3 and a share `absorbed` of it has been
    taken out as SO3 by the absorbers upstream."""
    F_air = np.asarray(F_air, float)
    n_s = np.asarray(n_s, float)
    o2 = AIR["O2"] * F_air - n_s * (1.0 + 0.5 * x_so3)
    return {"N2": AIR["N2"] * F_air, "Ar": AIR["Ar"] * F_air, "CO2": AIR["CO2"] * F_air, "O2": o2,
            "SO2": n_s * (1.0 - x_so3), "SO3": n_s * (x_so3 - absorbed)}


def total(s):
    return sum(s[k] for k in SPECIES)


def enthalpy(s, T):
    return sum(s[k] * h_abs(k, T) for k in SPECIES)


def bed_outlet_T(s_in, s_out, T_in, loss_share):
    """Adiabatic outlet temperature of a catalyst bed; the shell loses `loss_share` of the reaction heat."""
    H_in = enthalpy(s_in, T_in)
    q_rx = H_in - enthalpy(s_out, T_in)
    target = H_in - loss_share * q_rx
    T = np.asarray(T_in, float) + 0.0
    for _ in range(30):
        f = enthalpy(s_out, T) - target
        dT = 0.05
        g = (enthalpy(s_out, T + dT) - enthalpy(s_out, T - dT)) / (2 * dT)
        T = T - f / g
    return T


def log10_kp(T):
    """SO2 + 1/2 O2 = SO3, Kp in atm^-1/2."""
    return 4956.0 / T - 4.678


def equilibrium_ratio_needed(s, P_atm):
    y = {k: s[k] / total(s) for k in SPECIES}
    return y["SO3"] / (y["SO2"] * np.sqrt(y["O2"] * P_atm))
