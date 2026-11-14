"""Shaft circuit of generator unit GT-ST 11 (520 MVA, static excitation), per harmonic, as the plant study states it.

Nodes: T (shaft at the turbine-end bearing), X (shaft at the exciter-end bearing), ground.
  T-ground : grounding brush (contact resistance r_brush) in series with the grounding module resistor,
             in parallel with the turbine-end bearing oil film (r_film || c_film).
  X-ground : exciter-end bearing insulation (r_ins || c_ins).
  T -> X   : the shaft, carrying the axial EMF e_ax (raises X above T) through r_shaft.
Injections: the electrostatic charging current from the LP blading enters at T (DC only); the static
exciter's ripple couples through the field-to-body capacitance and enters at X (harmonics 3 and 6).
Monitors: DE voltage brush (V_T), grounding-brush current (I_gb), NDE voltage brush (V_X, high-impedance input).
"""
import math

import numpy as np

HARMONICS = (0, 1, 3, 6)
NOMINAL = dict(
    r_shaft=0.004,
    r_module=3.0,
    r_brush=1.2,
    r_film=8000.0,
    c_film=30e-9,
    r_ins=250e3,
    c_ins=8e-9,
    c_field=0.30e-6,
    axial_v_per_ka=0.045,
    third_v_per_ka_field=0.12,
    charge_a_per_kg_s=1.2e-5,
    monitor_r_in=1.0e6,
    monitor_c_in=2.0e-9,
)


def par(z1, z2):
    return z1 * z2 / (z1 + z2)


def cap(c, w):
    return 1.0 / (1j * w * c) if w > 0 else 1e18


def solve(h, p, e_ax, i_x, i_t, contact=0.0):
    """Phasors (V_T, I_gb, V_X as the NDE monitor reads it) for harmonic h with element values p."""
    w = 2 * math.pi * 50.0 * h
    z_gb = p["r_brush"] + p["r_module"]
    z_t = par(z_gb, par(p["r_film"], cap(p["c_film"], w)))
    z_mon = par(p["monitor_r_in"], cap(p["monitor_c_in"], w))
    z_x = par(par(p["r_ins"], cap(p["c_ins"], w)), contact + z_mon)
    # unknowns V_T, V_X, I_s (shaft current T->X); V_X - V_T = e_ax - I_s r_shaft
    a = np.array([[1 / z_t, 0, 1], [0, 1 / z_x, -1], [-1, 1, p["r_shaft"]]], dtype=complex)
    b = np.array([i_t, i_x, e_ax], dtype=complex)
    v_t, v_x, _ = np.linalg.solve(a, b)
    i_gb = v_t / z_gb
    v_read = v_x * z_mon / (contact + z_mon)
    return v_t, i_gb, v_read


def sources(h, p, stator_ka, field_ka, ripple_v, lp_steam_kg_s):
    """Axial EMF, injection at X and injection at T for harmonic h."""
    w = 2 * math.pi * 50.0 * h
    if h == 0:
        return 0.0, 0.0, p["charge_a_per_kg_s"] * lp_steam_kg_s
    if h == 1:
        return p["axial_v_per_ka"] * stator_ka, 0.0, 0.0
    if h == 3:
        return p["third_v_per_ka_field"] * field_ka, 0.0, 0.0
    return 0.0, 1j * w * p["c_field"] * ripple_v, 0.0
