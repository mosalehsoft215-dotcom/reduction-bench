"""Recomputes a pipeline's polarized potential profile for a set of rectifier currents, from the true line
parameters in line_truth.json (the parameters the survey records were generated with).

Steady direct-current model on a 10 m node grid with insulated ends:
  * pipe: longitudinal resistance r = rho_steel / (pi (D - t) t) per metre;
  * each ground bed: a point current source at the surface; beside the pipe it raises the soil by
    rho_bed I / (2 pi distance);
  * at every node the soil delivers J = G (drive - eta) through the coating, drive = soil potential - pipe
    potential, and the steel/soil interface carries J = K (exp(eta / b) - 1) (cathodic Tafel branch);
    G and K are the section values times 10 m;
  * a shorted casing drains a further G_S * drive at its node, without polarization;
  * the polarized (IR-free) potential is E_native - eta, in volts against copper/copper-sulphate.
The nonlinear node equations are solved by Newton's method on the pipe potential, with the interface
polarization found by bisection inside each step.
"""
import math

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

NODE_M = 10.0


def _sections(edges, values, x):
    idx = np.searchsorted(np.asarray(edges, float), x, side="right") - 1
    return np.asarray(values, float)[idx]


def _interface(drive, G, K, b):
    low, high = np.zeros_like(drive), drive.copy()
    for _ in range(64):
        mid = 0.5 * (low + high)
        surplus = G * (drive - mid) - K * np.expm1(np.minimum(mid / b, 60.0))
        low = np.where(surplus > 0, mid, low)
        high = np.where(surplus > 0, high, mid)
    return 0.5 * (low + high)


def polarized_profile(line, currents):
    """Polarized potentials (V CSE) at every node for the currents (dict bed name -> A)."""
    x = np.arange(0.0, line["length_m"] + 1e-9, NODE_M)
    n = len(x)
    r = line["steel_resistivity"] / (math.pi * (line["od_m"] - line["wall_m"]) * line["wall_m"])
    G = _sections(line["section_starts_m"], line["g_S_per_m"], x) * NODE_M
    K = _sections(line["section_starts_m"], line["k_A_per_m"], x) * NODE_M
    short = np.zeros(n)
    for point in line["casing_shorts"]:
        short[int(round(point["chainage_m"] / NODE_M))] += point["G_S"]
    native = np.interp(x, line["native_chainage_m"], line["native_V"])
    soil = np.zeros(n)
    drain = np.zeros(n)
    for bed in line["beds"]:
        amps = currents[bed["name"]]
        soil += bed["rho_ohm_m"] * amps / (2 * math.pi * np.hypot(x - bed["bed_chainage_m"], bed["offset_m"]))
        drain[int(round(bed["drain_chainage_m"] / NODE_M))] -= amps
    link = 1.0 / (r * NODE_M)
    centre = np.zeros(n)
    centre[:-1] += link
    centre[1:] += link
    ladder = diags([-link * np.ones(n - 1), centre, -link * np.ones(n - 1)], [-1, 0, 1], format="csc")
    b = line["tafel_b_V"]
    pipe = -0.3 - soil
    for _ in range(100):
        drive = np.maximum(soil - pipe, 0.0)
        eta = _interface(drive, G, K, b)
        pickup = G * (drive - eta) + short * (soil - pipe)
        mismatch = ladder @ pipe - (pickup + drain)
        slope = K / b * np.exp(np.minimum(eta / b, 60.0))
        tangent = ladder + diags(G * slope / (G + slope) * (soil - pipe > 0) + short, 0, format="csc")
        delta = spsolve(tangent, -mismatch)
        pipe = pipe + delta
        if np.max(np.abs(delta)) < 1e-11:
            break
    eta = _interface(np.maximum(soil - pipe, 0.0), G, K, b)
    return x, native - eta
