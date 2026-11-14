#!/usr/bin/env python3
"""Collision cross sections of the two analyte ions of one drift-tube acquisition, in He and in N2.

    python3 ccs_reduce.py <acquisition folder> <answer.json>

1. Arrival time of every ion in every frame: centroid of a Gaussian fitted to its arrival-time distribution, with a
   flat baseline.
2. Each gas run, each ion: t_A = t0 + L / v_d over the six drift voltages, with v_d = K E / sqrt(1 + M v_d^2 / 3kT)
   (the ion's energy gain in the field) and K = K0 (760 / P) (T / 273.15) from the frame's own pressure and
   temperature; t0 and K0 are fitted.
3. Mason-Schamp turns K0 into a cross section with the gas temperature of the run.
4. The two reference ions, whose cross sections in each pure gas are tabulated, check the run. A reading error of
   pressure or temperature scales every ion's result by one factor. Foreign gas in the drift tube does not: by Blanc's
   law each ion's mobility in a mixture is the harmonic mean of its pure-gas mobilities, which changes with the ion.
   The N2 run gets one factor from the two references. The He run gets the factor and the N2 fraction x that make
   both references right, and each analyte's He mobility is recovered from its mixed mobility with its N2 mobility
   from the N2 run.
"""
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

E = 1.602176634e-19
KB = 1.380649e-23
AMU = 1.66053906660e-27
TORR = 133.322368
N0 = 760 * TORR / (KB * 273.15)
GAS_MASS = {"He": 4.002602, "N2": 28.0134}
TABLE = {"TBA": {"He": 113.6, "N2": 160.4}, "THA": {"He": 165.8, "N2": 213.9}}


def ms_k0(omega, z, mass, gas, t):
    """Reduced mobility (m^2/Vs) of an ion of cross section omega (A^2) by Mason-Schamp."""
    mu = mass * GAS_MASS[gas] / (mass + GAS_MASS[gas]) * AMU
    return 3 * z * E / (16 * N0) * math.sqrt(2 * math.pi / (mu * KB * t)) / (omega * 1e-20)


def ms_omega(k0, z, mass, gas, t):
    mu = mass * GAS_MASS[gas] / (mass + GAS_MASS[gas]) * AMU
    return 3 * z * E / (16 * N0) * math.sqrt(2 * math.pi / (mu * KB * t)) / k0 * 1e20


def centroid(t_ms, counts):
    t = np.asarray(t_ms) * 1e-3
    y = np.asarray(counts, float)
    k = int(np.argmax(y))
    p0 = [y[k] - np.median(y), t[k], max((t[-1] - t[0]) / 12, 1e-6), np.median(y)]

    def res(p):
        model = p[0] * np.exp(-0.5 * ((t - p[1]) / p[2]) ** 2) + p[3]
        return (model - y) / np.sqrt(np.maximum(model, 1.0))
    r = least_squares(res, p0, x_scale=[p0[0], p0[2] / 10, p0[2], 1.0])
    return r.x[1], abs(r.x[2])


def drift_time(k0, volts, length_m, p_torr, t_k, gas, heating=True):
    k = k0 * (760.0 / p_torr) * (t_k / 273.15)
    e_field = volts / length_m
    m = GAS_MASS[gas] * AMU
    v = k * e_field
    for _ in range(40 if heating else 0):
        v = k * e_field / math.sqrt(1.0 + m * v * v / (3 * KB * t_k))
    return length_m / v


def fit_run(rows, length_m, gas, heating=True):
    """rows: list of (volts, p, t, t_arrival); returns K0 and t0 from the stepped-field fit."""
    V, P, T, TA = map(np.array, zip(*rows))
    slope = np.polyfit((P / T) / V, TA, 1)[0]                       # linear start
    k0_guess = length_m ** 2 * 273.15 / (760.0 * slope)

    def res(q):
        t0, k0 = q[0], q[1] * k0_guess
        return np.array([t0 + drift_time(k0, v, length_m, p, t, gas, heating) for v, p, t in zip(V, P, T)]) - TA
    r = least_squares(res, [float(np.polyfit((P / T) / V, TA, 1)[1]), 1.0], x_scale=[1e-4, 1e-3])
    return r.x[1] * k0_guess, r.x[0]


def width_excess(rows, t0, z, gate_s):
    """Median ratio of the fitted width to the width a single ion species can have: the gate and the diffusion
    along the tube, sigma^2 = gate^2 + 2 D t_d / v_d^2 with D = K k_B T / (z e) and K = v_d / E."""
    ratios = []
    for v, p, t, ta, w in rows:
        td = ta - t0
        L = LENGTH_HOLDER[0]
        vd = L / td
        k = vd / (v / L)
        d = k * KB * t / (z * E)
        ratios.append(w / math.sqrt(gate_s ** 2 + 2 * d * td / vd ** 2))
    return float(np.median(ratios))


LENGTH_HOLDER = [0.7824]


def reduce(folder, references="full", heating=True, width_check=True):
    """references: "full" (scale and N2 fraction), "scale" (one factor per run) or "none" (first principles);
    the switches exist for the evidence scripts."""
    folder = Path(folder)
    sheet = json.loads((folder / "acquisition.json").read_text())
    length_m = sheet["drift_length_cm"] / 100.0
    ions = {i["ion"]: i for i in sheet["ions"]}
    frames = {}
    with open(folder / "frames.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            frames[(r["gas"], int(r["frame"]))] = (float(r["drift_v"]), float(r["pressure_torr"]),
                                                   float(r["temperature_k"]))
    atd = defaultdict(lambda: ([], []))
    with open(folder / "atd.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            key = (r["gas"], int(r["frame"]), r["ion"])
            atd[key][0].append(float(r["t_ms"]))
            atd[key][1].append(float(r["counts"]))
    k0 = {}
    temps = {}
    broad = {}
    LENGTH_HOLDER[0] = length_m
    gate_s = sheet["trap_gate_us"] * 1e-6 / math.sqrt(12.0)
    for gas in ("N2", "He"):
        temps[gas] = float(np.mean([v[2] for (g, _), v in frames.items() if g == gas]))
        for ion in ions:
            rows = []
            for (g, j), (v, p, t) in frames.items():
                if g == gas:
                    rows.append((v, p, t) + centroid(*atd[(gas, j, ion)]))
            k0[(gas, ion)], t0 = fit_run([r[:4] for r in rows], length_m, gas, heating)
            broad[(gas, ion)] = width_excess(rows, t0, ions[ion]["z"], gate_s)

    def mass(ion):
        return ions[ion]["mz"] * ions[ion]["z"]

    tn, th = temps["N2"], temps["He"]
    refs = list(TABLE)

    def table_k0(i, gas, t):
        return ms_k0(TABLE[i][gas], ions[i]["z"], mass(i), gas, t)

    def solve_run(gas, t, other_pure):
        """Factor s and foreign-gas fraction x of one run from the two references (Blanc's law)."""
        other = "N2" if gas == "He" else "He"

        def mixed(i, x):
            return 1.0 / ((1.0 - x) / table_k0(i, gas, t) + x / table_k0(i, other, t))

        def ratio(x):
            return (k0[(gas, refs[0])] / mixed(refs[0], x)) / (k0[(gas, refs[1])] / mixed(refs[1], x)) - 1.0
        x = 0.0
        clean = [r for r in refs if broad[(gas, r)] < 1.12 or not width_check]
        if references == "full" and len(clean) == 1:
            # one reference carries an unresolved neighbour (its peak is wider than one species can be): the
            # other one, with the logged conditions taken as they stand, fixes the foreign-gas fraction
            r1 = clean[0]
            target = k0[(gas, r1)]
            f = lambda x: 1.0 / ((1.0 - x) / table_k0(r1, gas, t) + x / table_k0(r1, other, t)) - target
            x = brentq(f, 0.0, 0.3) if f(0.0) * f(0.3) < 0 else 0.0
            s = 1.0
        elif references == "full" and abs(ratio(0.0)) >= 4e-4 and ratio(0.0) * ratio(0.3) < 0:
            x = brentq(ratio, 0.0, 0.3)
        if references == "full" and len(clean) == 1:
            pass
        elif references == "none":
            s = 1.0
        elif references == "scale":
            s = float(np.mean([k0[(gas, r)] / table_k0(r, gas, t) for r in refs]))
        elif len(clean) != 1:
            s = k0[(gas, refs[0])] / mixed(refs[0], x)
        pure = {}
        for i in ("a", "b"):
            k_mix = k0[(gas, i)] / s
            k_other = other_pure[i] * math.sqrt(other_pure["t"] / t)
            pure[i] = (1.0 - x) / (1.0 / k_mix - x / k_other)
        pure["t"] = t
        return pure, s, x

    # each run's analytes in the pure gas need the other gas's pure mobility; two passes settle it
    pure_he = {"a": k0[("He", "a")], "b": k0[("He", "b")], "t": th}
    pure_n2 = {"a": k0[("N2", "a")], "b": k0[("N2", "b")], "t": tn}
    for _ in range(3):
        pure_n2, s_n2, x_n2 = solve_run("N2", tn, pure_he)
        pure_he, s_he, x_he = solve_run("He", th, pure_n2)
    out = {}
    for tag in ("a", "b"):
        z, m = ions[tag]["z"], mass(tag)
        out[f"ccs_he_{tag}_A2"] = round(ms_omega(pure_he[tag], z, m, "He", th), 3)
        out[f"ccs_n2_{tag}_A2"] = round(ms_omega(pure_n2[tag], z, m, "N2", tn), 3)
    out["n2_fraction_in_he_run"] = round(x_he, 5)
    out["he_fraction_in_n2_run"] = round(x_n2, 5)
    out["scale_n2_run"] = round(float(s_n2), 5)
    out["scale_he_run"] = round(float(s_he), 5)
    return out


def main():
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1]), indent=1) + "\n")


if __name__ == "__main__":
    main()
