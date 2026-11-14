"""Physics of the synthetic Knudsen effusion records (generator side).

Vapour pressure of a crystal form with a constant heat-capacity difference dCp = Cp(gas) - Cp(crystal):
    ln p(T) = ln p298 - dH298 / R (1/T - 1/T0) + dCp / R [ (T0/T - 1) + ln(T/T0) ],   T0 = 298.15 K
A second crystal form b reached through a transition at Tt (enthalpy h taken up on heating, heat capacity Cp_b):
    ln p_b(T) = ln p_a(T) - [ (h - c Tt)(1/Tt - 1/T) + c ln(T/Tt) ] / R,   c = Cp_b - Cp_a
(equal pressures at Tt; the sublimation enthalpies differ by the transition enthalpy at every temperature).
Mass leaving the cell per second through an orifice of effusion conductance C (transmission probability times
orifice area, m^2): dm/dt = p C sqrt(M / (2 pi R T)), T the crystal temperature.
A liner (a closed cup whose lid is the orifice foil) is integrated in series.py: it takes heat from the block
through its conductance and loses the sublimation heat and the radiation of its orifice to the 77 K shroud.
"""
import numpy as np

R = 8.314462618
T0 = 298.15


def ln_p(T, dH298, p298, dCp):
    T = np.asarray(T, dtype=float)
    return np.log(p298) - dH298 / R * (1 / T - 1 / T0) + dCp / R * ((T0 / T - 1) + np.log(T / T0))


def ln_p_form_b(T, lnp_a, Tt_K, h, c):
    T = np.asarray(T, dtype=float)
    return lnp_a - ((h - c * Tt_K) * (1 / Tt_K - 1 / T) + c * np.log(T / Tt_K)) / R


def program(T_start_C, T_set_C, hold_s, up_K_min, down_K_min, T_end_C, dt=1.0):
    """Block temperature of one run's program on a fine grid (true values, K)."""
    t_up = (T_set_C - T_start_C) / up_K_min * 60.0
    t_down = (T_set_C - T_end_C) / down_K_min * 60.0
    total = t_up + hold_s + t_down
    t = np.linspace(0.0, total, int(round(total / dt)) + 1)
    T = np.where(t < t_up, T_start_C + up_K_min * t / 60.0,
                 np.where(t < t_up + hold_s, T_set_C, T_set_C - down_K_min * (t - t_up - hold_s) / 60.0))
    return t, T + 273.15


def rate(T_K, lnp, C, M_kg):
    """Mass leaving per second (kg/s) at temperature T_K for ln p values lnp."""
    return np.exp(lnp) * C * np.sqrt(M_kg / (2 * np.pi * R * T_K))


def clausing_W(L_over_r):
    """Transmission probability of a short cylindrical channel (Berman's closed form of Clausing's integral)."""
    x = float(L_over_r)
    if x == 0:
        return 1.0
    s = np.sqrt(x * x + 4)
    num = (8 - x * x) * s + x ** 3 - 16
    den = 72 * x * s - 288 * np.log(x / 2 + s / 2)
    return 1 + x * x / 4 - x / 4 * s - num ** 2 / den


def dH_sub(T, dH298, dCp):
    """Sublimation enthalpy (J/mol) at T for a constant dCp = Cp(gas) - Cp(crystal)."""
    return dH298 + dCp * (np.asarray(T, dtype=float) - T0)


def crystal_temperature(T_block, lnp_fun, dH_fun, C, M_kg, G):
    """Temperature of a crystal held in a liner that exchanges heat with the block at G (W/K): the heat the
    sublimation takes, rate * dH / M, flows in through G, so T = T_block - rate(T) dH(T) / (M G).
    Without a liner (G None) the crystal is at the block temperature."""
    T_block = np.asarray(T_block, dtype=float)
    if G is None:
        return T_block
    T = T_block.copy()
    for _ in range(40):
        q = rate(T, lnp_fun(T), C, M_kg) * dH_fun(T) / M_kg
        T_new = T_block - q / G
        if np.max(np.abs(T_new - T)) < 1e-10:
            return T_new
        T = T_new
    return T
