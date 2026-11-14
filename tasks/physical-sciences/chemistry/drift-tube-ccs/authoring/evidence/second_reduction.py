#!/usr/bin/env python3
"""A second reduction of one drift-tube acquisition, written apart from solution/ccs_reduce.py.

    python second_reduction.py <acquisition folder> <answer.json>

- Arrival times: first moment of each distribution above its baseline (median of the outer bins), over the bins
  within three widths of the maximum.
- Stepped field: straight line of t_A against (P / T) / V for every ion and run (the method note: low fields), so
  K0 = L^2 273.15 / (760 slope).
- Peak width: a reference whose distribution is wider than one species can be (the gate and the diffusion along
  the tube, sigma^2 = gate^2 + t_d^2 2 k T / (z e V)) is left out, and its run keeps the logged conditions.
- Both runs at once: each run has a reading factor s (pressure, temperature, length) and a foreign-gas fraction
  (N2 in the He run, He in the N2 run). The four reference equations (two ions, two runs) fix the four unknowns
  with scipy.optimize.fsolve. For each analyte, Blanc's law in the two runs is linear in 1/K_He and 1/K_N2, so its
  two pure-gas mobilities come from one 2x2 solve.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

QE, KB, U = 1.602176634e-19, 1.380649e-23, 1.66053906660e-27
N0 = 2.686780e25
GAS = {"He": 4.002602, "N2": 28.0134}
REF = {"TBA": (113.6, 160.4), "THA": (165.8, 213.9)}


def k0_of(omega, z, m, gas, t):
    mu = m * GAS[gas] / (m + GAS[gas]) * U
    return 3 * z * QE / (16 * N0) * math.sqrt(2 * math.pi / (mu * KB * t)) / (omega * 1e-20)


def omega_of(k0, z, m, gas, t):
    return k0_of(1.0, z, m, gas, t) / k0


def main(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "acquisition.json").read_text())
    L = sheet["drift_length_cm"] / 100
    ions = {d["ion"]: (d["z"], d["mz"] * d["z"]) for d in sheet["ions"]}
    fr = np.genfromtxt(folder / "frames.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    at = np.genfromtxt(folder / "atd.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    k0, temp, wide = {}, {}, {}
    gate = sheet["trap_gate_us"] * 1e-6 / math.sqrt(12)
    for gas in ("He", "N2"):
        sel = fr["gas"] == gas
        temp[gas] = float(fr["temperature_k"][sel].mean())
        for ion in ions:
            xs, ys, sds, vs, ts = [], [], [], [], []
            for row in fr[sel]:
                m = (at["gas"] == gas) & (at["frame"] == row["frame"]) & (at["ion"] == ion)
                t, c = at["t_ms"][m] * 1e-3, at["counts"][m].astype(float)
                base = np.median(np.concatenate([c[:8], c[-8:]]))
                y = np.clip(c - base, 0, None)
                mean = np.sum(t * y) / np.sum(y)
                sd = math.sqrt(np.sum((t - mean) ** 2 * y) / np.sum(y))
                win = np.abs(t - t[np.argmax(c)]) <= 3 * sd
                xs.append(row["pressure_torr"] / row["temperature_k"] / row["drift_v"])
                ys.append(np.sum(t[win] * y[win]) / np.sum(y[win]))
                sds.append(math.sqrt(np.sum((t[win] - ys[-1]) ** 2 * y[win]) / np.sum(y[win])))
                vs.append(row["drift_v"])
                ts.append(row["temperature_k"])
            slope, icpt = np.polyfit(xs, ys, 1)
            k0[(gas, ion)] = L * L * 273.15 / (760 * slope)
            # resolving-power check: a single species has sigma^2 = gate^2 + t_d^2 * 2 k T / (z e V), whatever K is
            # (the moment over +-3 sd of a Gaussian keeps 0.986 of its width)
            z = ions[ion][0]
            pred = [math.sqrt(gate ** 2 + (ta - icpt) ** 2 * 2 * KB * tt / (z * QE * v)) for ta, v, tt in zip(ys, vs, ts)]
            wide[(gas, ion)] = float(np.median(np.array(sds) / (0.986 * np.array(pred))))

    def pure(ion, gas):
        z, m = ions[ion]
        om = REF[ion][0 if gas == "He" else 1]
        return k0_of(om, z, m, gas, temp[gas])

    usable = {gas: [r for r in REF if wide[(gas, r)] < 1.12] for gas in ("He", "N2")}

    def equations(q):
        s_n2, y_he, s_he, x_n2 = q
        out = []
        # a run that lost a reference to an unresolved neighbour keeps its logged conditions (factor 1)
        if len(usable["N2"]) < 2:
            out.append(s_n2 - 1.0)
        if len(usable["He"]) < 2:
            out.append(s_he - 1.0)
        for ion in REF:
            inv_n2run = (1 - y_he) / pure(ion, "N2") + y_he / (pure(ion, "He") * math.sqrt(temp["He"] / temp["N2"]))
            inv_herun = (1 - x_n2) / pure(ion, "He") + x_n2 / (pure(ion, "N2") * math.sqrt(temp["N2"] / temp["He"]))
            if ion in usable["N2"]:
                out.append(k0[("N2", ion)] * inv_n2run / s_n2 - 1.0)
            if ion in usable["He"]:
                out.append(k0[("He", ion)] * inv_herun / s_he - 1.0)
        return out
    s_n2, y_he, s_he, x_n2 = fsolve(equations, [1.0, 0.0, 1.0, 0.0], xtol=1e-13)
    answer = {}
    for tag in ("a", "b"):
        z, m = ions[tag]
        r = math.sqrt(temp["He"] / temp["N2"])
        # unknowns u = 1/K_He(at T_He), w = 1/K_N2(at T_N2); K at the other run's temperature scales with sqrt(T)
        A = np.array([[y_he / r, 1 - y_he], [1 - x_n2, x_n2 * r]])
        b = np.array([s_n2 / k0[("N2", tag)], s_he / k0[("He", tag)]])
        u, w = np.linalg.solve(A, b)
        answer[f"ccs_he_{tag}_A2"] = omega_of(1 / u, z, m, "He", temp["He"])
        answer[f"ccs_n2_{tag}_A2"] = omega_of(1 / w, z, m, "N2", temp["N2"])
    return answer


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(main(sys.argv[1])) + "\n")
