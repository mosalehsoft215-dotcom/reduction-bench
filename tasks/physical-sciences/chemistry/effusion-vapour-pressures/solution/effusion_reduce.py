"""Knudsen effusion series -> sublimation enthalpy and log10 vapour pressure at 298.15 K.

    python3 effusion_reduce.py <series folder> <out.json>

1. Thermocouple contact losses in the holds (a sudden drop, later a sudden rise back) are bridged first.  Each run's
   loss is the Hertz-Knudsen rate integrated over the cell's own temperature record, from the start of
   pumping to venting: dm = k exp(a - b (1/T - 1/T0)) / sqrt(T) dt with k = W A sqrt(M / (2 pi R)), W the Clausing
   factor of the orifice.  For a trial b the level a follows in closed form; b is found by a bounded 1-D search on the
   weighted misfit of ln(loss).
2. A setpoint that was run twice is checked: when the earlier run lost clearly more than the later one predicts
   (more than four standard errors of the pair, plus half a per cent), the earlier run is left out and the fit is
   repeated.
3. dH = b R; log10 p(298.15 K) from the fitted line.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

RG = 8.314462618
SIG_MG = 0.0035          # weighing noise of a loss (two weighings)


def repair(minute, temp, run):
    """During the hold the cell temperature changes slowly; a sudden drop of more than 1.5 K followed later in the hold
    by a sudden rise back is the thermocouple losing and regaining contact with the cell base: the readings in between
    are replaced by the line joining the last reading before the drop and the first after the recovery."""
    temp = temp.copy()
    inside = np.flatnonzero((minute > run["hold_start_min"] + 2.0) & (minute <= run["hold_end_min"]))
    if inside.size < 3:
        return temp
    step = np.diff(temp)
    i = inside[0]
    while i < inside[-1]:
        if step[i] < -1.5:
            ups = [j for j in range(i + 1, inside[-1]) if step[j] > 1.5]
            if ups:
                j = ups[0]
                temp[i + 1:j + 1] = np.interp(minute[i + 1:j + 1], [minute[i], minute[j + 1]], [temp[i], temp[j + 1]])
                i = j
        i += 1
    return temp


def load(folder, fix=True):
    folder = Path(folder)
    sheet = json.loads((folder / "series.json").read_text())
    with open(folder / "runs.csv", newline="") as f:
        runs = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]
    for r in runs:
        data = np.loadtxt(folder / "traces" / f"run-{int(r['run']):02d}.csv", delimiter=",", skiprows=1)
        r["t"], r["T"] = data[:, 0] * 60.0, (repair(data[:, 0], data[:, 1], r) if fix else data[:, 1])
        r["loss"] = r["mass_before_mg"] - r["mass_after_mg"]
    return sheet, runs


def reduce_series(folder, trace=True, repeats=True, fix=True):
    sheet, runs = load(folder, fix)
    d = sheet["orifice_diameter_mm"]
    W = 1.0 / (1.0 + 3.0 * sheet["foil_thickness_um"] * 1e-3 / (8.0 * d / 2.0))
    k = W * math.pi * (d * 1e-3) ** 2 / 4.0 * math.sqrt(sheet["molar_mass_g_mol"] * 1e-3 / (2 * math.pi * RG)) * 1e6
    T0 = 1.0 / np.mean([1.0 / r["setpoint_k"] for r in runs])
    for r in runs:
        if trace:
            dt = np.gradient(r["t"])
            r["u"], r["w"] = 1.0 / r["T"] - 1.0 / T0, dt / np.sqrt(r["T"])
        else:
            hold = (r["hold_end_min"] - r["hold_start_min"]) * 60.0
            r["u"], r["w"] = np.array([1.0 / r["setpoint_k"] - 1.0 / T0]), np.array([hold / math.sqrt(r["setpoint_k"])])
        r["sig"] = SIG_MG / r["loss"]

    def resid(b, use):
        z = np.array([math.log(r["loss"] / (k * np.sum(r["w"] * np.exp(-b * r["u"])))) for r in use])
        wt = np.array([1.0 / r["sig"] ** 2 for r in use])
        a = np.sum(wt * z) / np.sum(wt)
        return a, z - a, wt

    def fit(use):
        b = minimize_scalar(lambda b: float(np.sum(resid(b, use)[2] * resid(b, use)[1] ** 2)),
                            bounds=(2000.0, 30000.0), method="bounded", options={"xatol": 1e-3}).x
        return b, resid(b, use)[0]

    use = list(runs)
    b, a = fit(use)
    if repeats:
        for _ in range(3):
            drop = set()
            by_set = {}
            for r in use:
                by_set.setdefault(round(r["setpoint_k"], 2), []).append(r)
            for group in by_set.values():
                if len(group) < 2:
                    continue
                group.sort(key=lambda r: r["run"])
                z = {id(r): math.log(r["loss"] / (k * np.sum(r["w"] * np.exp(-b * r["u"])))) for r in group}
                last = group[-1]
                for r in group[:-1]:
                    if z[id(r)] - z[id(last)] > 4.0 * math.hypot(r["sig"], last["sig"]) + 0.005:
                        drop.add(int(r["run"]))
            if not drop:
                break
            use = [r for r in use if int(r["run"]) not in drop]
            b, a = fit(use)
    ln298 = a - b * (1.0 / 298.15 - 1.0 / T0)
    return {"sublimation_enthalpy_kj_mol": b * RG / 1000.0, "log10_pressure_298_pa": ln298 / math.log(10.0)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_series(sys.argv[1])) + "\n")
