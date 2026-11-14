"""Reduces one incremental hole-drilling record to the principal residual stresses.

    python3 hole_reduce.py <hole directory> <result.json>

1. Steps: the record is split into its depth holds; the last 30 s of every hold (after the cutting heat has gone)
   are used.
2. Temperature: the rosette is compensated for one expansion coefficient; on a specimen with another, every gauge
   reads (alpha_specimen - alpha_compensation) microstrain per kelvin of the specimen thermocouple, removed from
   every reading relative to the start.
3. Strains: p = (e1 + e3)/2, q = (e3 - e1)/2, t = (e3 + e1 - 2 e2)/2 at each depth, from the hold means.
4. Stresses: least squares of p = -P (1 + nu) a(h) / E, q = -Q b(h) / E, t = -T b(h) / E over the 20 depths, the
   coefficient functions taken at the measured hole diameter; sigma = P +- sqrt(Q^2 + T^2).
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


def coeffs(h, d):
    s = (d / 1.8) ** 2
    return 0.155 * s * (1 - np.exp(-h / 0.42)), 0.32 * s * (1 - np.exp(-h / 0.30))


def reduce_hole(folder, thermal=True, diameter="measured"):
    folder = Path(folder)
    sh = json.loads((folder / "sheet.json").read_text())
    with open(folder / "bridge.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    depth = np.array([float(r["depth_mm"]) for r in rows])
    air = np.array([float(r["specimen_c"]) for r in rows])
    gauges = np.array([[float(r[f"gauge{j}_ue"]) for j in (1, 2, 3)] for r in rows])
    # holds: runs of constant depth; keep their last 30 samples
    edges = np.flatnonzero(np.diff(depth) != 0) + 1
    runs = [r for r in np.split(np.arange(depth.size), edges) if len(r) >= 30]
    holds = [r[-30:] for r in runs]
    # a self-temperature-compensated gauge on a specimen whose expansion differs from its compensation reads an
    # apparent strain (alpha_specimen - alpha_compensation) per kelvin of specimen temperature change
    stc = int(re.search(r"S-T-C\s*(\d+)", sh["rosette"]).group(1))   # maker's number: matched expansion in ppm/F
    mismatch = sh["thermal_expansion_ppm_per_k"] - 1.8 * stc
    kappa = np.full(3, mismatch if thermal else 0.0)
    corrected = gauges - np.outer(air - air[0], kappa)
    level = np.array([corrected[r].mean(axis=0) for r in holds])
    h = np.array([depth[r[0]] for r in holds])
    zero = level[h == 0].mean(axis=0) if np.any(h == 0) else np.zeros(3)
    rel = level[h > 0] - zero
    h = h[h > 0]
    e1, e2, e3 = rel.T
    p, q, t = (e1 + e3) / 2, (e3 - e1) / 2, (e3 + e1 - 2 * e2) / 2
    d = sh["hole_diameter_measured_mm"] if diameter == "measured" else sh["cutter_diameter_mm"]
    a, b = coeffs(h, d)
    e = sh["youngs_modulus_gpa"] * 1e3 * 1e-6          # MPa per microstrain
    nu = sh["poisson_ratio"]
    P = -np.dot(a, p) / np.dot(a, a) * e / (1 + nu)
    Q = -np.dot(b, q) / np.dot(b, b) * e
    T = -np.dot(b, t) / np.dot(b, b) * e
    r = math.hypot(Q, T)
    return {"max_principal_mpa": float(P + r), "min_principal_mpa": float(P - r)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_hole(sys.argv[1])))
