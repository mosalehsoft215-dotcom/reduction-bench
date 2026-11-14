"""The true HP circuits (circuits_truth.json) blown at a given drum pressure: drum steam flow and the smallest
cleaning force ratio along SH1, SH2 and the main steam line.

Steady continuous blow:
  * saturated steam leaves the drum at the drum pressure;
  * exhaust gas (flow m_g and specific heat cp_g, both fixed at the blowing load) meets SH2 first, then SH1;
  * each superheater is a counter-flow bank with a fixed overall conductance UA (kW/K); its duty solves
    UA = integral of m dh / (T_gas - T_steam) over the steam enthalpy rise, the gas temperature rising by
    m (h - h_in) / (m_g cp_g) from the gas outlet at the steam-inlet end, the steam pressure falling linearly
    with enthalpy;
  * every section loses P_in - P_out = k m^2 v (Pa), v at the section's mean pressure and mean enthalpy;
  * the attemperator mixes in the spray that passes its shut valve (feedwater at the pump discharge);
  * the exit at the end of the temporary pipe is a nozzle discharging to the barometric pressure p_b:
    m = C sqrt(P / v) psi(p_b / P), psi = 1 when choked (p_b / P at or below the critical ratio of an ideal
    gas with isentropic exponent 1.3) and the isentropic flow function scaled to 1 at the critical ratio above it.
Steam properties are IAPWS-IF97 (region 2 by its backward T(P, h) equations; regions 1, 2 and 4 when wet).
Cleaning force ratio K = (m / m_MCR)^2 v / v_MCR on 41 points per section, both states taken linearly in pressure
and enthalpy between their section ends.
"""
import math
from functools import lru_cache

import numpy as np
from iapws import IAPWS97
from iapws.iapws97 import Const, _Region1, _Region2
from scipy.optimize import brentq

POINTS = 41
EXIT_EXPONENT = 1.3
NODES, WEIGHTS = np.polynomial.legendre.leggauss(24)
_N4 = (0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5, -0.32325550322333e7,
       0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849, 0.65017534844798e3)
_SUB = {"a": (Const.Backward2a_T_Ph_n, Const.Backward2a_T_Ph_Li, Const.Backward2a_T_Ph_Lj, 0.0, 2.1),
        "b": (Const.Backward2b_T_Ph_n, Const.Backward2b_T_Ph_Li, Const.Backward2b_T_Ph_Lj, -2.0, 2.6),
        "c": (Const.Backward2c_T_Ph_n, Const.Backward2c_T_Ph_Li, Const.Backward2c_T_Ph_Lj, 25.0, 1.8)}


def sat_temperature(P_bar):
    """IF97 region 4 saturation temperature (C); works on arrays."""
    n = _N4
    beta = (np.asarray(P_bar, dtype=float) / 10.0) ** 0.25
    E = beta ** 2 + n[2] * beta + n[5]
    F = n[0] * beta ** 2 + n[3] * beta + n[6]
    G = n[1] * beta ** 2 + n[4] * beta + n[7]
    D = 2 * G / (-F - np.sqrt(F ** 2 - 4 * E * G))
    return (n[9] + D - np.sqrt((n[9] + D) ** 2 - 4 * (n[8] + n[9] * D))) / 2 - 273.15


def _backward_T(P, h):
    """Region-2 backward T(P, h) in C for superheated states (arrays of bar and kJ/kg)."""
    P = P / 10.0
    hbc = 0.26526571908428e4 + np.sqrt(np.clip(P - 0.45257578905948e1, 0.0, None) / 1.2809002730136e-4)
    region = np.where(P <= 4.0, "a", np.where((P <= 6.546699678) | (h >= hbc), "b", "c"))
    T = np.empty(h.shape)
    for key, (n, Li, Lj, shift, nu0) in _SUB.items():
        sel = region == key
        if sel.any():
            T[sel] = (n * (P[sel, None] + shift) ** Li * (h[sel, None] / 2000.0 - nu0) ** Lj).sum(axis=1)
    return np.maximum(T - 273.15, sat_temperature(P * 10.0))


@lru_cache(maxsize=None)
def _saturation(P_key):
    Ts = float(sat_temperature(P_key)) + 273.15
    g, f = _Region2(Ts, P_key / 10.0), _Region1(Ts, P_key / 10.0)
    return {"T": Ts - 273.15, "hg": g["h"], "hf": f["h"], "vg": g["v"], "vf": f["v"]}


def saturation(P_bar):
    return _saturation(round(float(P_bar), 9))


def T_ph(P_bar, h):
    """Steam temperature (C): IF97 region-2 backward equations when superheated, saturation when wet."""
    s = saturation(P_bar)
    if h <= s["hg"]:
        return s["T"]
    return float(_backward_T(np.array([float(P_bar)]), np.array([float(h)]))[0])


def v_ph(P_bar, h):
    """Specific volume (m3/kg): region 2 above saturation, the two-phase mixture below it."""
    s = saturation(P_bar)
    if h > s["hg"]:
        return _Region2(T_ph(P_bar, h) + 273.15, P_bar / 10.0)["v"]
    return s["vf"] + (h - s["hf"]) / (s["hg"] - s["hf"]) * (s["vg"] - s["vf"])


def h_pT(P_bar, T_C):
    return IAPWS97(P=P_bar / 10.0, T=T_C + 273.15).h


def exit_factor(ratio, k=EXIT_EXPONENT):
    """Isentropic nozzle mass flux at back-pressure ratio p_b/P over the choked flux (1 when choked)."""
    critical = (2.0 / (k + 1.0)) ** (k / (k - 1.0))
    if ratio <= critical:
        return 1.0
    if ratio >= 1.0:
        return 0.0
    flux = 2.0 * k / (k - 1.0) * (ratio ** (2.0 / k) - ratio ** ((k + 1.0) / k))
    choked = k * (2.0 / (k + 1.0)) ** ((k + 1.0) / (k - 1.0))
    return math.sqrt(flux / choked)


def outlet_pressure(P_in, h_in, h_out, flow, k):
    """P_out = P_in - k m^2 v(mean P, mean h); None when the pressure runs out."""
    P_out = P_in
    for _ in range(300):
        new = P_in - k * flow * flow * v_ph(0.5 * (P_in + P_out), 0.5 * (h_in + h_out)) * 1e-5
        if new < 0.5:
            return None
        if abs(new - P_out) < 1e-10:
            return new
        P_out = new
    return None


def _conductance_needed(P_in, P_out, h_in, h_out, flow, gas_in, gas_capacity):
    """UA (kW/K) a counter-flow bank needs to lift the steam from h_in to h_out; huge if the gas is not hotter."""
    gas_out = gas_in - flow * (h_out - h_in) / gas_capacity
    hg = saturation(P_in)["hg"]
    edges = [h_in, hg, h_out] if h_in < hg < h_out else [h_in, h_out]
    total = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        h = 0.5 * (a + b) + 0.5 * (b - a) * NODES
        P = P_in + (h - h_in) / (h_out - h_in) * (P_out - P_in)
        steam = sat_temperature(P) if b <= hg else _backward_T(P, h)
        gap = gas_out + flow * (h - h_in) / gas_capacity - steam
        if np.any(gap <= 0.0):
            return 1e12
        total += 0.5 * (b - a) * float(np.dot(WEIGHTS, flow / gap))
    return total


def bank(P_in, h_in, flow, UA, k, gas_in, gas_capacity):
    """Steam outlet (P, h) and gas leaving temperature of a counter-flow superheater; None if pressure runs out."""
    P_out = P_in - k * flow * flow * v_ph(P_in, h_in) * 1e-5
    if P_out < 0.5:
        return None
    for _ in range(60):
        def short(h_out):
            return _conductance_needed(P_in, P_out, h_in, h_out, flow, gas_in, gas_capacity) - UA
        hi = h_in + 20.0
        while short(hi) < 0.0:
            hi = h_in + 1.6 * (hi - h_in)
        h_out = brentq(short, h_in + 1e-6, hi, xtol=1e-8, rtol=1e-13)
        new = outlet_pressure(P_in, h_in, h_out, flow, k)
        if new is None:
            return None
        if abs(new - P_out) < 1e-9:
            P_out = new
            break
        P_out = new
    return P_out, h_out, gas_in - flow * (h_out - h_in) / gas_capacity


def path(c, Pd, flow):
    """States along the path for a drum steam flow; 'excess' > 0 when the choked exit cannot pass that flow."""
    gas_capacity = c["gas_flow_kg_s"] * c["gas_cp_kJ_kgK"]
    drum = saturation(Pd)
    leak = c["passing_spray_kg_s"]
    flow2 = flow + leak
    water = h_pT(Pd + c["feed_over_drum_bar"], c["feed_temperature_C"])

    def through(gas_mid):
        sh1 = bank(Pd, drum["hg"], flow, c["UA1"], c["k1"], gas_mid, gas_capacity)
        if sh1 is None:
            return None
        hm = (flow * sh1[1] + leak * water) / flow2
        sh2 = bank(sh1[0], hm, flow2, c["UA2"], c["k2"], c["gas_in_C"], gas_capacity)
        return None if sh2 is None else (sh1, hm, sh2)

    # the gas meets SH2 first but SH2 takes its steam from SH1: secant on the gas temperature between the banks
    x0, x1 = c["gas_in_C"] - 30.0, c["gas_in_C"] - 60.0
    r0 = through(x0)
    if r0 is None:
        return None
    f0 = r0[2][2] - x0
    for _ in range(40):
        r1 = through(x1)
        if r1 is None:
            return None
        f1 = r1[2][2] - x1
        if abs(f1) < 1e-7:
            break
        x0, x1, f0 = x1, x1 - f1 * (x1 - x0) / (f1 - f0), f1
    (P1, h1, _), hm, (P2, h2, _) = r1
    P3 = outlet_pressure(P2, h2, h2, flow2, c["k3"])
    P4 = P3 and outlet_pressure(P3, h2, h2, flow2, c["k4"])
    if not P4:
        return None
    choked = c["C"] * math.sqrt(P4 * 1e5 / v_ph(P4, h2)) * exit_factor(c["ambient_bar"] / P4)
    return {"excess": flow2 - choked, "flow": flow, "flow2": flow2, "Pd": Pd, "P1": P1, "h0": drum["hg"],
            "h1": h1, "hm": hm, "P2": P2, "h2": h2, "P3": P3, "P4": P4, "gas_mid": x1,
            "sections": [("SH1", Pd, P1, drum["hg"], h1, flow), ("SH2", P1, P2, hm, h2, flow2),
                         ("MS", P2, P3, h2, h2, flow2)]}


@lru_cache(maxsize=None)
def _mcr_volumes(P_a, T_a, P_b, T_b, points):
    ha = saturation(P_a)["hg"] if T_a is None else h_pT(P_a, T_a)
    hb = h_pT(P_b, T_b)
    f = np.linspace(0.0, 1.0, points)
    return tuple(IAPWS97(P=(P_a + t * (P_b - P_a)) / 10.0, h=ha + t * (hb - ha)).v for t in f)


def solve(c, Pd):
    """Path states of circuit c at drum pressure Pd (bar abs), at the flow the choked exit passes."""
    def excess(flow):
        state = path(c, Pd, flow)
        return 1e3 if state is None else state["excess"]
    lo = 4.0
    while excess(lo * 1.15) < 0.0:
        lo *= 1.15
    return path(c, Pd, brentq(excess, lo, lo * 1.15, xtol=1e-9))


def smallest_ratio(c, state, points=None):
    points = points or POINTS
    worst = math.inf
    for name, Pa, Pb, ha, hb, m in state["sections"]:
        sec = c["mcr"][name]
        v_mcr = _mcr_volumes(sec["inlet"][0], sec["inlet"][1], sec["outlet"][0], sec["outlet"][1], points)
        for i in range(points):
            f = i / (points - 1)
            ratio = (m / sec["flow_kg_s"]) ** 2 * v_ph(Pa + f * (Pb - Pa), ha + f * (hb - ha)) / v_mcr[i]
            worst = min(worst, ratio)
    return worst


def blow(c, Pd):
    """(drum steam flow kg/s, smallest cleaning force ratio) of circuit c at drum pressure Pd (bar abs)."""
    state = solve(c, Pd)
    return state["flow"], smallest_ratio(c, state)
