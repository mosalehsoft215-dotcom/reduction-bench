"""The plain one-pass reading: documented fields used as written, nothing cross-checked.

Time zero at the shutter command; the cell thermistor as the liquid temperature; the certified alpha of every
reference; eta averaged over the references; alpha = P_abs / (eta P l) (absorbed power proportional to alpha).

    python authoring/proofs/plain_reading.py <day folder> <answer.json>
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

LP = 632.8e-9


def rows(p):
    return list(csv.DictReader(open(p, newline="")))


def main(folder, out):
    d = Path(folder)
    s = json.loads((d / "setup.json").read_text())
    pr = s["probe"]
    w0 = pr["waist_radius_um"] * 1e-6
    zc = math.pi * w0 ** 2 / LP
    z1 = (pr["sample_position_mm"] - pr["waist_position_mm"]) * 1e-3
    we = s["excitation"]["waist_radius_um_at_sample"] * 1e-6
    m, V = w0 ** 2 * (1 + (z1 / zc) ** 2) / we ** 2, z1 / zc
    f = lambda t, th, tc: (1 - th / 2 * np.arctan(2 * m * V / (((1 + 2 * m) ** 2 + V ** 2) * tc / (2 * t) + 1 + 2 * m + V ** 2))) ** 2
    props = rows(d / "properties.csv")
    def kd(sol, T):
        r = [x for x in props if x["solvent"] == sol]
        x = [float(a["t_C"]) for a in r]
        return np.interp(T, x, [float(a["k_W_per_m_K"]) for a in r]), np.interp(T, x, [float(a["dn_dT_per_K"]) for a in r])
    etas, sam = [], []
    for r in rows(d / "runs.csv"):
        tr = rows(d / "transients" / f"{r['run']}.csv")
        t = np.array([float(x["t_ms"]) for x in tr]) / 1000
        p = np.array([float(x["probe_V"]) for x in tr])
        (th, tc), _ = curve_fit(f, t[t > 0], p[t > 0] / p[t < 0].mean(), p0=[0.3, 0.008], bounds=([0, 1e-5], [5, 1]))
        liq = s["liquids"][r["liquid"]]
        k, dn = kd(liq["solvent"], float(r["cell_C"]))
        P_abs = -th * k * LP / dn
        P, l = float(r["pump_meter_mW"]) * 1e-3, float(r["path_cm"])
        if liq["role"] == "reference":
            etas.append(P_abs / (P * liq["certified_alpha_per_cm"] * l))
        else:
            sam.append((P_abs, P, l, tc))
    eta = float(np.mean(etas))
    alpha = float(np.mean([pa / (eta * P * l) for pa, P, l, _ in sam]))
    D = float(np.mean([we ** 2 / (4 * tc) for *_, tc in sam])) * 1e6
    Path(out).write_text(json.dumps({"alpha_per_cm": alpha, "diffusivity_mm2_per_s": D}) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
