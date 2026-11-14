"""Forward model of a power-compensation titration calorimeter, used to write the run folders.

Physics, in the order it is applied:

* Binding: one class of independent sites, one site per active macromolecule; only the ligand monomer binds.
* Ligand self-association, where the compound has a dimer entry: monomer and dimer in equilibrium in every solution
  (Kdim = m^2 / d, dissociation enthalpy dHdis per mole of dimer). Free monomer from the mass balance
  L = m + 2 m^2/Kdim + S m/(Kd + m), solved by bracketing.
* Cell bookkeeping: the cell is overfilled, so its sensed volume V stays fixed. An injection of v pushes the same
  volume of the cell's pre-injection solution out of the sensed region before it mixes: X -> X (1 - v/V) for the
  cell species, Y -> Y (1 - v/V) + Y_s v/V for the titrant. The heat of injection i is
  dH V (C_i - C_{i-1} k) - dHdis (V (D_i - D_{i-1} k) - v D_s) + h_mol Y_s v, with k = 1 - v/V, D_s the dimer
  concentration of the syringe solution and h_mol the titrant's heat of dilution per mole injected.
* Every fill draws its titrant from its own stock (its own Y_s and D_s).
* Before the first injection of each fill, titrant diffuses out of the needle into the cell during the
  equilibration (a leak of lam microlitres); that heat is released before the record starts and the first
  injection then delivers v_1 - lam.
* The delivered volume of every injection scatters around its nominal value (relative sigma vol_sd); where an
  operator remark gives delivered volumes (a bubble), those are the volumes delivered.
* Heat release: the titrant arrives at a constant rate over the injection's duration and mixes with a first-order
  lag tau_mix; the instrument turns the heat flow into its reading through a first-order lag tau_inst that belongs
  to the record's feedback mode. An electrical calibration pulse is Joule heat in the cell wall and sees tau_inst only.
* Reading: dp = baseline - scale * (heat flow released in the cell), the mean over each filter period, plus white
  noise; scale is the instrument's power scale on that day, the same in every record of the day.
"""
import numpy as np
from scipy.optimize import brentq

CAL_J = 4.184


def complex_conc(sites, lig, kd):
    b = sites + lig + kd
    return 0.5 * (b - np.sqrt(np.maximum(b * b - 4.0 * sites * lig, 0.0)))


def equilibrium(sites, lig, kd, kdim):
    """(complex, dimer) concentrations for total sites and total ligand (monomer units)."""
    if kdim is None:
        return complex_conc(sites, lig, kd), 0.0
    if lig <= 0.0:
        return 0.0, 0.0
    f = lambda m: m + 2.0 * m * m / kdim + sites * m / (kd + m) - lig
    m = brentq(f, 0.0, lig, xtol=1e-15, rtol=1e-14)
    return sites * m / (kd + m), m * m / kdim


def heats_per_injection(run, fills):
    """Heat (microcal, negative when released) of every injection, fill by fill in order.

    run: V_ul, forward (bool), cell_um (cell species total; 0 for a blank), active_fraction, kd_um, dh_kj,
    dimer (None or (kdim_um, dh_dissociation_kj)), dil_ucal_per_nmol.
    fills: list of dicts with titrant_um, leak_ul and delivered (array of delivered volumes).
    """
    V = run["V_ul"]
    f = run["active_fraction"]
    kdim, dhdis = run["dimer"] if run["dimer"] else (None, 0.0)
    X, Y = run["cell_um"], 0.0

    def state(x, y):
        sites, lig = (f * x, y) if run["forward"] else (f * y, x)
        return equilibrium(sites, lig, run["kd_um"], kdim if run["forward"] else None)

    out = []
    for fill in fills:
        Ys = fill["titrant_um"]
        Ds = equilibrium(0.0, Ys, run["kd_um"], kdim)[1] if (run["forward"] and kdim) else 0.0
        lam = fill["leak_ul"]
        keep = 1.0 - lam / V
        X, Y = X * keep, Y * keep + Ys * lam / V
        C, D = state(X, Y)
        q_fill = []
        for k, v in enumerate(fill["delivered"]):
            v_in = v - lam if k == 0 else v
            keep = 1.0 - v_in / V
            Xn, Yn = X * keep, Y * keep + Ys * v_in / V
            Cn, Dn = state(Xn, Yn)
            q = (run["dh_kj"] * V * (Cn - C * keep) - dhdis * (V * (Dn - D * keep) - v_in * Ds)) * 1e-3 / CAL_J
            q += run["dil_ucal_per_nmol"] * Ys * v_in * 1e-3
            q_fill.append(q)
            X, Y, C, D = Xn, Yn, Cn, Dn
        out.append(np.array(q_fill))
    return out


def boxcar_response(t, start, dur, tau_a, tau_b):
    """Unit-area response at times t to a constant release over [start, start+dur] through two first-order lags
    (tau_b may be 0 for a single lag)."""
    x = np.asarray(t, float) - start

    def ramp_out(u):
        u = np.maximum(u, 0.0)
        if tau_b <= 0:
            return 1.0 - np.exp(-u / tau_a)
        if abs(tau_a - tau_b) < 1e-6:
            return 1.0 - (1.0 + u / tau_a) * np.exp(-u / tau_a)
        return 1.0 - (tau_a * np.exp(-u / tau_a) - tau_b * np.exp(-u / tau_b)) / (tau_a - tau_b)
    return (ramp_out(x) - ramp_out(x - dur)) / dur


def record(times_end, period, events, baseline, scale, sd, rng, sub=10):
    """Readings at the end of each filter period: baseline - scale * mean heat flow over the period + noise.

    events: list of (heat_ucal, start_s, duration_s, tau_a, tau_b), heat as the enthalpy change of the cell
    contents: negative when heat is released (binding with dH < 0, or the Joule heat of the calibration pulse).
    """
    acc = np.zeros_like(times_end)
    offs = (np.arange(sub) + 0.5) / sub * period
    for heat, start, dur, ta, tb in events:
        lo, hi = start - 1.0, start + dur + 40.0 * max(ta, tb) + 50.0
        m = (times_end > lo) & (times_end - period < hi)
        if not m.any():
            continue
        tt = times_end[m]
        r = np.zeros_like(tt)
        for o in offs:
            r += boxcar_response(tt - period + o, start, dur, ta, tb)
        acc[m] += heat * r / sub
    # heat released in the cell lets the feedback take power off the sample heater, so dp falls
    return baseline + scale * acc + rng.normal(0.0, sd, len(times_end))


def drift_curve(t, p):
    span = max(t[-1], 1.0)
    u = t / span
    return (p["b0"] + p["b1"] * u + p["b2"] * u * u
            + p["w1"] * np.sin(2 * np.pi * t / p["P1"] + p["ph1"])
            + p["w2"] * np.sin(2 * np.pi * t / p["P2"] + p["ph2"]))
