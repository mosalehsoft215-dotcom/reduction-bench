"""Residual styrene in one polystyrene lot by multiple headspace extraction (MHE).

    python3 styrene_mhe.py <lot folder> <answer.json>

1. Peaks: every chromatogram of a vial is fitted together with its sister extractions: a straight baseline and the
   styrene peak (Gaussian, one width for the vial, its own place per extraction) and, where it lowers the misfit
   clearly, a second Gaussian at one fixed offset from styrene for the whole vial (a co-eluting product); areas and
   baselines are linear, width, offset and places nonlinear.
2. Pressure: gas that left the vial during the hold took styrene with it, so each area is put back to the
   concentration before the hold: A x start / end (the logged pressures).
3. MHE: ln A against the extraction number is a straight line (ratio Q per extraction); the area sum is A1 / (1 - Q).
4. Amount: area sums are proportional to amount / (V_gas (1 - p_a / p_start)); the empty standard vial (V_gas = the
   vial) gives the scale, the lot's own gas volume V_vial - m / rho its share.
5. Partition: Q = 1 - beta / (beta + K) (1 - p_a / p_start), beta = V_gas / V_polymer.
"""
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402


def read(folder):
    sheet = json.loads((folder / "lot_sheet.json").read_text())
    traces = {}
    with open(folder / "chromatograms.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            traces.setdefault((r["vial"], int(r["extraction"])), []).append((float(r["seconds"]), float(r["fid_pa"])))
    press = {}
    with open(folder / "pressure_checks.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            press[(r["vial"], int(r["extraction"]))] = (float(r["start_kpa"]), float(r["end_kpa"]))
    return sheet, {k: np.array(v) for k, v in traces.items()}, press


def gauss(s, mu, sig):
    return np.exp(-0.5 * ((s - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))


def vial_areas(traces, n, two):
    """Styrene areas of the n extractions of one vial (and the residual sum of squares of the joint fit)."""
    data = [traces[i + 1] for i in range(n)]
    tops = [d[np.argmax(d[:, 1]), 0] for d in data]

    def solve(p):
        sig, off, mus = p[0], p[1], p[2:]
        res, areas = [], []
        for d, mu in zip(data, mus):
            s, y = d[:, 0], d[:, 1]
            cols = [np.ones_like(s), s - 310.0, gauss(s, mu, sig)]
            if two:
                cols.append(gauss(s, mu + off, sig))
            B = np.vstack(cols).T
            c = np.linalg.lstsq(B, y, rcond=None)[0]
            res.append(B @ c - y)
            areas.append(c[2])
        return np.concatenate(res), np.array(areas)

    best = None
    for off0 in ((2.0, 3.0, 4.5) if two else (3.0,)):
        p0 = np.r_[1.2, off0, tops]
        lo = np.r_[0.5, 1.0, np.array(tops) - 3]
        hi = np.r_[3.0, 8.0, np.array(tops) + 3]
        sol = least_squares(lambda p: solve(p)[0], p0, bounds=(lo, hi), x_scale="jac")
        if best is None or sol.cost < best.cost:
            best = sol
    return solve(best.x)[1], 2 * best.cost


def area_sum(areas):
    i = np.arange(areas.size)
    slope, icept = np.polyfit(i, np.log(areas), 1)
    q = float(np.exp(slope))
    return float(np.exp(icept)) / (1 - q), q


def reduce_lot(folder, leak=True, deconvolve=True, gas_volume=True):
    folder = Path(folder)
    sheet, traces, press = read(folder)
    n = int(sheet["extractions"])
    baro = sheet["barometric_pressure_kpa"]
    out = {}
    for vial in ("standard", "lot"):
        tr = {i: traces[(vial, i)] for i in range(1, n + 1)}
        one, ss1 = vial_areas(tr, n, False)
        areas = one
        if deconvolve and vial == "lot":
            pair, ss2 = vial_areas(tr, n, True)
            noise = np.var(np.concatenate([np.diff(tr[i][:60, 1]) for i in tr])) / 2
            if ss1 - ss2 > 25.0 * n * noise:             # a co-eluting product: keep the two-peak fit
                areas = pair
        start = np.array([press[(vial, i)][0] for i in range(1, n + 1)])
        end = np.array([press[(vial, i)][1] for i in range(1, n + 1)])
        if leak:
            areas = areas * start / end
        total, q = area_sum(areas)
        out[vial] = (total, q, float(np.mean(start)))
    v_s = sheet["sample_mass_g"] / sheet["polymer_density_g_per_ml_at_oven"]
    v_vial = sheet["vial_volume_ml"]
    v_g = v_vial - v_s if gas_volume else v_vial
    (t_std, _, p_std), (t_lot, q_lot, p_lot) = out["standard"], out["lot"]
    amount = sheet["standard_styrene_ug"] * (t_lot / t_std) * (v_g / v_vial) * (1 - baro / p_lot) / (1 - baro / p_std)
    beta = (v_vial - v_s) / v_s
    k = beta * ((1 - baro / p_lot) / (1 - q_lot) - 1)
    return {"styrene_mg_per_kg": float(amount / sheet["sample_mass_g"]), "partition_coefficient": float(k)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_lot(sys.argv[1])) + "\n")
