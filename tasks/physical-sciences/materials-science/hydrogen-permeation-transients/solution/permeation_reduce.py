"""Lattice hydrogen diffusivity at 25 C and sub-surface hydrogen concentration from one permeation test.

    python3 /app/solution/permeation_reduce.py <specimen directory> <result file>

1. Background: a straight line through the exit current before charging starts, carried over the whole test.
2. Transients: the decay (charging off) and the second rise (charging on again) are fitted together with Fick's
   solutions for a membrane whose entry is switched between C0 and 0 and whose exit is held at 0:
   rise(t) = I_inf [1 + 2 sum (-1)^n exp(-n^2 pi^2 D t / L^2)], decay(t) = I_inf - rise(t); least squares in D and
   I_inf. The first rise is left out: in a cold-worked specimen it is slowed while irreversible traps fill, and only
   after that do the transients follow lattice diffusion.
3. Temperature: the mean of the cell thermocouple; D is referred to 25 C with the lattice activation energy,
   D25 = D exp(Ea/R (1/T - 1/298.15)).
4. Thickness: the mean of the micrometer readings of the tested membrane.
5. Sub-surface concentration: C0 = I_inf L / (F A D) at the test temperature.
PERMEATION_ROUTE (comma-separated) switches single steps off; authoring/evidence/bench.py replays them.
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

F, R_GAS = 96485.33212, 8.314462618
ROUTE = set(filter(None, os.environ.get("PERMEATION_ROUTE", "").split(",")))
N_TERMS = np.arange(1, 120)


def rise_shape(t, d, L):
    """Rise of the exit current over its steady value, from a membrane free of hydrogen."""
    t = np.maximum(np.asarray(t, float), 0.0)
    return 1 + 2 * np.sum(((-1.0) ** N_TERMS)[None, :] * np.exp(-np.outer(t, N_TERMS ** 2) * math.pi ** 2 * d / L ** 2), 1)


def reduce(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "specimen.json").read_text())
    data = np.loadtxt(folder / "transient.csv", delimiter=",", skiprows=1)
    t, cur, temp = data[:, 0], data[:, 1], data[:, 2]
    with open(folder / "events.csv", newline="") as fh:
        ev = [(float(r["time_s"]), r["event"]) for r in csv.DictReader(fh)]
    ons = [a for a, e in ev if e == "charging on"]
    off = [a for a, e in ev if e == "charging off"][0]
    pre = t < ons[0] - 5
    base = np.polyval(np.polyfit(t[pre], cur[pre], 1), t) if "no_background" not in ROUTE else np.full_like(t, np.median(cur[pre]))
    net = cur - base
    if "as_cut_thickness" in ROUTE:                # authoring shortcut: the sheet's as-cut thickness
        L = sheet["thickness_as_cut_mm"] * 1e-3
    else:                                          # the tested membrane: mean of the micrometer readings
        L = float(np.loadtxt(folder / "micrometer.csv", delimiter=",", skiprows=1)[:, 1].mean()) * 1e-3
    rise1 = (t > ons[0]) & (t < off)
    decay = (t > off) & (t < ons[1])
    rise2 = t > ons[1]

    def model(p, sel_rise1, sel_decay, sel_rise2):
        d, i_inf = p[0] * 1e-12, p[1]
        parts = []
        if sel_rise1:
            parts.append(i_inf * rise_shape(t[rise1] - ons[0], d, L) - net[rise1])
        if sel_decay:
            parts.append(i_inf * (1 - rise_shape(t[decay] - off, d, L)) - net[decay])
        if sel_rise2:
            parts.append(i_inf * rise_shape(t[rise2] - ons[1], d, L) - net[rise2])
        return np.concatenate(parts)

    use = (False, True, True)
    if "first_transient" in ROUTE:
        use = (True, False, False)
    if "all_transients" in ROUTE:
        use = (True, True, True)
    i_guess = float(np.median(net[rise2][-50:]))
    fit = least_squares(lambda p: model(p, *use), [500.0, i_guess], x_scale=[50.0, max(i_guess, 1e-3) * 0.1],
                        bounds=([1.0, 0.0], [1e5, np.inf]))
    d_test, i_inf = fit.x[0] * 1e-12, fit.x[1]
    T = (sheet["cell_setpoint_c"] if "setpoint_temperature" in ROUTE else float(np.mean(temp))) + 273.15
    ea = sheet["lattice_activation_energy_kj_per_mol"] * 1e3
    d25 = d_test * math.exp(ea / R_GAS * (1 / T - 1 / 298.15))
    c0 = i_inf * 1e-6 * L / (F * sheet["exposed_area_cm2"] * 1e-4 * d_test)
    return {"diffusivity_25c_um2_per_s": d25 * 1e12, "subsurface_hydrogen_mol_per_m3": c0}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + "\n", encoding="utf-8")
