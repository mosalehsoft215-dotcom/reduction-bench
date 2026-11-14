"""Second, separately written MHE reduction of one lot (same call as the reference, different numerics).

    python second_headspace.py <lot folder> <answer.json>

- Baseline of each trace: the straight line through the medians of its first and last 5 s, taken off before fitting.
- Peak width: from the standard vial's traces alone (second central moment of the baseline-free signal within 6 s of
  the maximum, corrected for that truncation by a fitted Gaussian), shared by every peak of the session.
- Lot traces: the co-eluting product's offset comes from the last extraction (where the product is largest), fitted
  with the styrene place, both heights and the offset free; it is kept when it lowers that trace's misfit by more than
  a third and the product's area is positive, and is then held for every extraction, whose styrene place and the two
  areas are refitted (least squares, Levenberg-Marquardt).
- Leak: area x start pressure / end pressure.
- MHE: A_i = a Q^(i-1) fitted to the areas themselves by least squares with relative weights (not the log-linear
  line of the reference); area sum a / (1 - Q).
- Amount and partition coefficient from the same relations as the notes.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def load(folder):
    sheet = json.loads((folder / "lot_sheet.json").read_text())
    rows = list(csv.reader(open(folder / "chromatograms.csv", newline="")))[1:]
    tr = {}
    for v, e, s, y in rows:
        tr.setdefault((v, int(e)), []).append((float(s), float(y)))
    pr = {(v, int(e)): (float(a), float(b)) for v, e, a, b in list(csv.reader(open(folder / "pressure_checks.csv")))[1:]}
    return sheet, {k: np.array(x) for k, x in tr.items()}, pr


def flatten(d):
    s, y = d[:, 0], d[:, 1]
    lo, hi = s < s[0] + 5.0, s > s[-1] - 5.0
    x0, x1 = s[lo].mean(), s[hi].mean()
    y0, y1 = np.median(y[lo]), np.median(y[hi])
    return s, y - (y0 + (y1 - y0) * (s - x0) / (x1 - x0))


def g(s, mu, sig, a):
    return a * np.exp(-0.5 * ((s - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))


def width(traces):
    sig = []
    for d in traces:
        s, y = flatten(d)
        mu0 = s[np.argmax(y)]
        sol = least_squares(lambda p: g(s, p[0], p[1], p[2]) - y, [mu0, 1.2, y.max() * 3.0], method="lm")
        sig.append(abs(sol.x[1]))
    return float(np.median(sig))


def fit_trace(s, y, sig, offset=None):
    mu0 = s[np.argmax(y)]
    if offset is None:
        p0 = [mu0, y.max() * 3.0]
        model = lambda p: g(s, p[0], sig, p[1])  # noqa: E731
    else:
        p0 = [mu0, y.max() * 3.0, y.max() * 0.3]
        model = lambda p: g(s, p[0], sig, p[1]) + g(s, p[0] + offset, sig, p[2])  # noqa: E731
    sol = least_squares(lambda p: model(p) - y, p0, method="lm")
    return sol.x, float(np.sum(sol.fun ** 2))


def styrene_areas(traces, sig, product):
    if product:
        s, y = flatten(traces[-1])
        one, ss1 = fit_trace(s, y, sig)
        sol = least_squares(lambda p: g(s, p[0], sig, p[1]) + g(s, p[0] + p[3], sig, p[2]) - y,
                            [one[0], one[1] * 0.85, one[1] * 0.15, 3.0], method="lm")
        ss2 = float(np.sum(sol.fun ** 2))
        if ss2 < ss1 * 2.0 / 3.0 and sol.x[2] > 0:
            off = float(sol.x[3])
            return np.array([fit_trace(*flatten(d), sig, off)[0][1] for d in traces])
    return np.array([fit_trace(*flatten(d), sig)[0][1] for d in traces])


def mhe_sum(areas):
    i = np.arange(areas.size)
    q0 = np.exp(np.polyfit(i, np.log(areas), 1)[0])
    sol = least_squares(lambda p: (p[0] * p[1] ** i - areas) / areas, [areas[0], q0], method="lm")
    a, q = sol.x
    return a / (1 - q), q


def run(folder):
    folder = Path(folder)
    sheet, tr, pr = load(folder)
    n, baro = int(sheet["extractions"]), sheet["barometric_pressure_kpa"]
    sig = width([tr[("standard", i)] for i in range(1, n + 1)])
    res = {}
    for vial in ("standard", "lot"):
        areas = styrene_areas([tr[(vial, i)] for i in range(1, n + 1)], sig, vial == "lot")
        p = np.array([pr[(vial, i)] for i in range(1, n + 1)])
        total, q = mhe_sum(areas * p[:, 0] / p[:, 1])
        res[vial] = (total, q, float(p[:, 0].mean()))
    vp = sheet["sample_mass_g"] / sheet["polymer_density_g_per_ml_at_oven"]
    vv = sheet["vial_volume_ml"]
    vg = vv - vp
    (ts, _, ps), (tl, ql, pl) = res["standard"], res["lot"]
    ug = sheet["standard_styrene_ug"] * tl / ts * vg / vv * (1 - baro / pl) / (1 - baro / ps)
    k = vg / vp * ((1 - baro / pl) / (1 - ql) - 1)
    return {"styrene_mg_per_kg": float(ug / sheet["sample_mass_g"]), "partition_coefficient": float(k)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(run(sys.argv[1])))
