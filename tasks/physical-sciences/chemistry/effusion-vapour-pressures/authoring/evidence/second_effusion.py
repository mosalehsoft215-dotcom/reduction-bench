"""Second, separately written reduction of a Knudsen effusion series (effective times and a straight line).

While the furnace ramps and holds, the cell only warms, so a hold reading more than 1 K under the highest reading
logged before it is a lost thermocouple contact; such readings are replaced by interpolation in time between the good
readings around them.  For a trial slope b every run gets an effective time at its setpoint Ts,
t_eff = integral of exp(-b (1/T - 1/Ts)) sqrt(Ts / T) dt over its temperature record, and a pressure
p = loss / (W A sqrt(M / (2 pi R Ts)) t_eff); a weighted straight line ln p against 1/Ts gives the next b, until b
stops moving.  Where a setpoint was run twice and the earlier pressure lies above the later one by more than three
standard errors of the pair plus one per cent, the earlier run is set aside.
"""
import json
import math
from pathlib import Path

import numpy as np

R = 8.314462618


def reduce_series(folder):
    folder = Path(folder)
    s = json.loads((folder / "series.json").read_text())
    tab = np.genfromtxt(folder / "runs.csv", delimiter=",", names=True)
    recs = [np.genfromtxt(folder / "traces" / f"run-{int(n):02d}.csv", delimiter=",", skip_header=1)
            for n in tab["run"]]
    for rc, h0, h1 in zip(recs, tab["hold_start_min"], tab["hold_end_min"]):
        # heating (pump-down, ramp, hold) is monotone for a lagged ramp-and-hold programme: a hold reading that sits
        # more than 1 K under the running maximum of the record so far is a lost contact
        heat = np.flatnonzero(rc[:, 0] <= h1)
        if heat.size < 8:
            continue
        top = np.maximum.accumulate(rc[heat, 1])
        bad = (rc[heat, 1] < top - 1.0) & (rc[heat, 0] > h0)
        if bad.any() and (~bad).sum() >= 2:
            good = heat[~bad]
            rc[heat[bad], 1] = np.interp(rc[heat[bad], 0], rc[good, 0], rc[good, 1])
    r_o = s["orifice_diameter_mm"] * 5e-4
    clausing = 1.0 / (1.0 + 3.0 * s["foil_thickness_um"] * 1e-6 / (8.0 * r_o))
    ts = tab["setpoint_k"]
    loss = (tab["mass_before_mg"] - tab["mass_after_mg"]) * 1e-6                  # kg
    rel = 3.5e-9 / loss
    conv = clausing * math.pi * r_o ** 2 * np.sqrt(s["molar_mass_g_mol"] * 1e-3 / (2 * math.pi * R * ts))

    def pressures(b):
        teff = np.array([np.trapezoid(np.exp(-b * (1 / rc[:, 1] - 1 / t)) * np.sqrt(t / rc[:, 1]), rc[:, 0] * 60.0)
                         for rc, t in zip(recs, ts)])
        return np.log(loss / (conv * teff))

    keep = np.ones(len(ts), bool)
    b = 12000.0
    for _ in range(4):
        for _ in range(60):
            y = pressures(b)
            x, w = 1.0 / ts[keep], 1.0 / rel[keep] ** 2
            xm, ym = np.sum(w * x) / np.sum(w), np.sum(w * y[keep]) / np.sum(w)
            slope = np.sum(w * (x - xm) * (y[keep] - ym)) / np.sum(w * (x - xm) ** 2)
            if abs(-slope - b) < 1e-4:
                break
            b = -slope
        y = pressures(b)
        new = keep.copy()
        for t in np.unique(ts):
            idx = np.flatnonzero(np.isclose(ts, t))
            if idx.size > 1:
                last = idx[-1]
                for i in idx[:-1]:
                    if y[i] - y[last] > 3.0 * math.hypot(rel[i], rel[last]) + 0.01:
                        new[i] = False
        if np.array_equal(new, keep):
            break
        keep = new
    lnp298 = ym - b * (1.0 / 298.15 - xm)
    return {"sublimation_enthalpy_kj_mol": b * R / 1000.0, "log10_pressure_298_pa": lnp298 / math.log(10.0)}


if __name__ == "__main__":
    import sys
    Path(sys.argv[2]).write_text(json.dumps(reduce_series(sys.argv[1])) + "\n")
