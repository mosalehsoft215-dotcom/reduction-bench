"""Forward model of a power-compensation titration calorimeter, used to write the run folders.

Physics, in the order it is applied:

* Binding: one class of independent sites, one site per active macromolecule. The complex concentration follows
  from the exact quadratic of the mass-action law.
* Cell bookkeeping: the cell is overfilled, so its sensed volume V stays fixed. An injection of v pushes the same
  volume of the cell's pre-injection solution up the access tube, out of the sensed region, before it mixes. Totals
  therefore go as X -> X (1 - v/V) for the cell species and Y -> Y (1 - v/V) + Y_syringe v/V for the titrant, and
  the heat of injection i is dH V (C_i - C_{i-1} (1 - v/V)) plus the heat of dilution and injection of the titrant.
* Before the first injection of each syringe fill, titrant has diffused out of the needle into the cell during the
  equilibration period (a leak of lam microlitres). That heat is released before the record starts; the first
  injection then delivers v_1 - lam, so the total reached after it equals the nominal one.
* The delivered volume of every injection scatters around its nominal value (relative sigma vol_sd).
* Heat release: the titrant arrives at a constant rate over the injection's duration and mixes with a first-order
  lag tau_mix. The instrument turns the heat flow into its reading through a first-order lag tau_inst.
  The electrical calibration pulse is Joule heat in the cell wall, so it sees tau_inst only.
* Reading: dp = baseline - scale * (heat flow released in the cell), recorded as the mean over each filter period,
  plus white noise. The baseline is a slow polynomial with a slow wander.
"""
import numpy as np

CAL_J = 4.184


def complex_conc(sites, lig, kd):
    b = sites + lig + kd
    return 0.5 * (b - np.sqrt(np.maximum(b * b - 4.0 * sites * lig, 0.0)))


def injection_heats(run, delivered):
    """Heat (microcal, negative when exothermic) of every injection, fill by fill in order.

    run: dict with V_ul, cell_species ('macromolecule' or 'ligand'), conc_cell_um, conc_syringe_um, active_fraction,
    kd_um, dh_kj, dil_ucal_per_ul, leaks_ul (one per fill), fills (list of lists of nominal volumes).
    delivered: list (per fill) of arrays of delivered volumes.
    """
    V = run["V_ul"]
    f = run["active_fraction"]
    forward = run["cell_species"] == "macromolecule"
    X = run["conc_cell_um"]          # cell species total
    Y = 0.0                          # titrant total in the cell
    Ys = run["conc_syringe_um"]

    def comp(x, y):
        sites, lig = (f * x, y) if forward else (f * y, x)
        return complex_conc(sites, lig, run["kd_um"])

    C = comp(X, Y)
    out = []
    for fill, vols in enumerate(delivered):
        lam = run["leaks_ul"][fill]
        keep = 1.0 - lam / V
        X, Y = X * keep, Y * keep + Ys * lam / V
        C = comp(X, Y)               # leak equilibrated before the record
        q_fill = []
        for k, v in enumerate(vols):
            v_in = v - lam if k == 0 else v
            keep = 1.0 - v_in / V
            Xn, Yn = X * keep, Y * keep + Ys * v_in / V
            Cn = comp(Xn, Yn)
            q = run["dh_kj"] * V * (Cn - C * keep) * 1e-3 / CAL_J + run["dil_ucal_per_ul"] * v_in
            q_fill.append(q)
            X, Y, C = Xn, Yn, Cn
        out.append(np.array(q_fill))
    return out


def lagged_box(t, start, dur, tau_a, tau_b):
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
            r += lagged_box(tt - period + o, start, dur, ta, tb)
        acc[m] += heat * r / sub
    # heat released in the cell lets the feedback take power off the sample heater, so dp falls
    return baseline + scale * acc + rng.normal(0.0, sd, len(times_end))


def baseline_curve(t, p):
    span = max(t[-1], 1.0)
    u = t / span
    return (p["b0"] + p["b1"] * u + p["b2"] * u * u
            + p["w1"] * np.sin(2 * np.pi * t / p["P1"] + p["ph1"])
            + p["w2"] * np.sin(2 * np.pi * t / p["P2"] + p["ph2"]))
