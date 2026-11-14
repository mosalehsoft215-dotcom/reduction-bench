"""Thermal history of a power transformer through a short-circuit temperature-rise test.

Lumped model of the IEC 60076-7 kind, driven by the losses and current actually applied:
  top-oil rise over the ambient the oil sees    tau_o dD/dt = D_ss - D,   D_ss = D_or (P/P_N)^x
  top-to-bottom oil difference                  tau_B dB/dt = B_ss - B,   B_ss = B_r  (P/P_N)^x
  winding-to-average-oil gradient (each)        tau_w dg/dt = g_ss - g,   g_ss = g_r  (I/I_N)^y
The ambient the oil sees is the air temperature passed through a first-order lag of tau_o; the IEC oil cups around
the tank are built with the same lag. tau_B = tau_o while losses flow and tau_o/2 after switch-off.
top = seen + D, average oil = top - B/2, bottom = top - B, winding = average oil + g.
Times in minutes, temperatures in degC. Everything runs on a 5 s grid.
"""
import math

import numpy as np

DT = 1.0 / 12.0  # minutes (5 s)


def lag(u, tau, y0, dt=DT):
    out = np.empty_like(u)
    y = y0
    for i, ui in enumerate(u):
        a = math.exp(-dt / tau[i] if np.ndim(tau) else -dt / tau)
        y = ui + (y - ui) * a
        out[i] = y
    return out


def simulate(p, power_kw, current_pu, air_c):
    """Arrays on the 10 s grid: applied losses (kW), current (per unit of rated), air temperature (degC)."""
    rel = np.clip(power_kw / p["P_N"], 0.0, None)
    seen = lag(air_c, p["tau_o"], air_c[0])
    D = lag(p["D_or"] * rel ** p["x"], p["tau_o"], 0.0)
    tau_B = np.where(power_kw > 0, p["tau_o"], p["tau_o"] / 2.0)
    B = lag(p["B_r"] * rel ** p["x"], tau_B, 0.0)
    top = seen + D
    avg = top - B / 2.0
    out = dict(seen=seen, D=D, B=B, top=top, avg=avg, bottom=top - B)
    for w in ("HV", "LV"):
        g = lag(p["g_r"][w] * np.clip(current_pu, 0.0, None) ** p["y"], p["tau_w"][w], 0.0)
        out["g_" + w] = g
        out["wind_" + w] = avg + g
    return out
