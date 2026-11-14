"""Second, separately written reduction of a hole-drilling record (same call as the reference).

Instead of hold means and strain combinations it fits the raw gauge samples at once: every sample from 20 s into each
hold onwards is modelled per gauge as an offset, plus the relieved strain of that gauge written out from (P, Q, T)
through the rosette angles and the coefficient functions at the measured diameter, plus the apparent strain
kappa (T_specimen - T_start), with kappa taken from the S-T-C number and the specimen's expansion.  The problem is
linear in (P, Q, T, three offsets), solved by least squares.
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


def reduce_hole(folder):
    folder = Path(folder)
    meta = json.loads((folder / "sheet.json").read_text())
    rec = np.genfromtxt(folder / "bridge.csv", delimiter=",", names=True)
    depth, temp = rec["depth_mm"], rec["specimen_c"]
    strains = np.column_stack([rec["gauge1_ue"], rec["gauge2_ue"], rec["gauge3_ue"]])
    stc = float(re.search(r"S-T-C\s*(\d+)", meta["rosette"]).group(1)) * 1.8
    kappa = meta["thermal_expansion_ppm_per_k"] - stc
    strains = strains - kappa * (temp - temp[0])[:, None]
    # samples belonging to a hold, at least 20 s after its depth was reached
    since = np.zeros(depth.size)
    for i in range(1, depth.size):
        since[i] = since[i - 1] + 1 if depth[i] == depth[i - 1] else 0
    keep = (since >= 20) | (depth == 0)
    h = depth[keep]
    d = meta["hole_diameter_measured_mm"]
    scale = (d / 1.8) ** 2
    a = 0.155 * scale * (1 - np.exp(-h / 0.42))
    b = 0.32 * scale * (1 - np.exp(-h / 0.30))
    e = meta["youngs_modulus_gpa"] * 1e-3              # MPa -> microstrain factor is 1/e
    nu = meta["poisson_ratio"]
    # gauge strains: e1 = p - q, e2 = p - t, e3 = p + q with p = -P(1+nu)a/E, q = -Q b/E, t = -T b/E
    rows, rhs = [], []
    n = h.size
    for g, (cq, ct) in enumerate([(-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)]):
        block = np.zeros((n, 6))
        block[:, 0] = -(1 + nu) * a / e
        block[:, 1] = -cq * b / e
        block[:, 2] = -ct * b / e
        block[:, 3 + g] = 1.0
        rows.append(block)
        rhs.append(strains[keep, g])
    sol = np.linalg.lstsq(np.vstack(rows), np.concatenate(rhs), rcond=None)[0]
    P, Q, T = sol[:3]
    r = math.hypot(Q, T)
    return {"max_principal_mpa": float(P + r), "min_principal_mpa": float(P - r)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_hole(sys.argv[1])))
