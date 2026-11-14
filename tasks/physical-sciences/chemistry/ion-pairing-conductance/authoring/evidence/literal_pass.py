"""The literal one-pass reduction: every documented field is read, every reading is taken at face value.

Temperature from the bath thermistor; R_inf by a straight line in 1/f; cell constant from the mean of every
standard measured that day; solvent conductivity from the pure-solvent fill (the blank for a dilution run);
molality from the weighings with the cell taken as closed; the guide's equation fitted for L0 and K_A.

    python authoring/evidence/literal_pass.py <series folder> <answer.json>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def lam(c, L0, KA, eps, eta_P, T):
    e = eps * T
    S, A = 8.204e5 / e ** 1.5 * L0 + 82.50 / (eta_P * e ** 0.5), 1.8246e6 / e ** 1.5
    out = []
    for ci in c:
        a = 1.0
        for _ in range(80):
            q = KA * ci * 10 ** (-2 * A * (a * ci) ** 0.5)
            a = 2.0 / (1.0 + (1.0 + 4.0 * q) ** 0.5)
        out.append(a * (L0 - S * (a * ci) ** 0.5))
    return np.array(out)


def main(folder, out):
    d = Path(folder)
    sheet = json.loads((d / "series.json").read_text())
    br = list(csv.DictReader(open(d / "bridge.csv")))
    T = float(np.mean([float(r["bath_C"]) for r in br]))
    g = {}
    for r in br:
        g.setdefault((r["fill"], int(r["step"])), []).append((1.0 / float(r["f_Hz"]), float(r["R_ohm"])))
    rinf = {k: np.polyfit([p[0] for p in v], [p[1] for p in v], 1)[1] for k, v in g.items()}
    std = list(csv.DictReader(open(d / "standards.csv")))
    tc = [float(r["t_C"]) for r in std]
    Ks = [np.interp(T, tc, [float(r[sheet["standards"][f.split(":")[1]]["column"]]) for r in std]) * 1e-6 * R
          for (f, _), R in rinf.items() if f.startswith("standard:")]
    K = float(np.mean(Ks))
    ms = {(r["entry"], r["quantity"]): float(r["grams"]) for r in csv.DictReader(open(d / "masses.csv"))}
    M = sheet["salt"]["molar_mass_g_per_mol"]; D = sheet["salt"]["solution_density_slope_kg_dm3_per_mol_kg"]
    sv = list(csv.DictReader(open(d / "solvent_properties.csv")))
    ts = [float(r["t_C"]) for r in sv]
    rho0, eps, eta = (np.interp(T, ts, [float(r[k]) for r in sv]) for k in ("density_g_cm3", "permittivity", "viscosity_mPa_s"))
    mol, solv, steps = 0.0, 0.0, {}
    if sheet["procedure"] == "concentration":
        w = ms[("stock", "salt_g")] / (ms[("stock", "salt_g")] + ms[("stock", "solvent_g")])
        solv = ms[("cell", "solvent_g")]
        i = 1
        while (f"add{i}", "stock_g") in ms:
            mol += ms[(f"add{i}", "stock_g")] * w / M; solv += ms[(f"add{i}", "stock_g")] * (1 - w)
            steps[i] = mol * 1000 / solv; i += 1
        ks = K / rinf[("series", 0)]
    else:
        mol, solv = ms[("cell", "salt_g")] / M, ms[("cell", "solvent_g")]
        steps[0] = mol * 1000 / solv
        i = 1
        while (f"add{i}", "solvent_g") in ms:
            solv += ms[(f"add{i}", "solvent_g")]; steps[i] = mol * 1000 / solv; i += 1
        ks = K / rinf[("blank", 0)]
    c = np.array([m * (rho0 + D * m) / (1 + m * M / 1000) for m in steps.values()])
    L = np.array([1000 * (K / rinf[("series", s)] - ks) / ci for s, ci in zip(steps, c)])
    fit = least_squares(lambda p: lam(c, p[0], p[1], eps, eta * 0.01, T + 273.15) - L, [L.max() * 1.05, 100.0],
                        bounds=([1, 0], [1e4, 1e6]), x_scale=[10, 100])
    Path(out).write_text(json.dumps({"lambda0_S_cm2_per_mol": fit.x[0], "ka_dm3_per_mol": fit.x[1]}) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
