"""Second reduction, written apart from the reference: every second of every run fitted with the disc's radial modes.

    python authoring/evidence/modal_solver.py <batch folder> <answer.json>

Differences from solution/absorption_calorimetry.py:
  * the heating stretch is used, not set aside: the rim response is computed from radial heat diffusion (Bessel
    modes of an insulated disc of the logged diameter, weights of a 2 mm Gaussian spot, or the uniform film);
  * the thermistor's scattered-light offset is a fitted term proportional to the monitor while the shutter is open;
  * the unheated disc relaxes towards the recorded wall from a fitted start; gamma, start, scatter and the rise per
    joule are fitted together with Levenberg-Marquardt on log gamma;
  * power = monitor times the median of the meter / monitor ratios (monitor median within 8 s of each reading);
  * windows: weighted least squares of absorptance against thickness (weights 1 / thickness), the mirror directly.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.special import j0, jn_zeros

KAPPA = 8.4e-7
BESSEL = np.r_[0.0, jn_zeros(1, 36)]


def mode_weights(radius, spot):
    if not spot:
        return np.r_[1.0, np.zeros(len(BESSEL) - 1)]
    w = 2.0e-3
    reach = 7 * w
    total = quad(lambda r: math.exp(-2 * r * r / w / w) * r, 0, reach)[0]
    vals = [quad(lambda r: math.exp(-2 * r * r / w / w) * j0(b * r / radius) * r, 0, reach, limit=300)[0] / total
            for b in BESSEL]
    return np.array(vals) / j0(BESSEL) ** 2


def response(t, power, gamma, radius, weights):
    """Rim temperature per unit (rise per joule) for the absorbed-power record `power`."""
    rates = KAPPA * BESSEL ** 2 / radius ** 2 + gamma
    keep = np.exp(-rates)
    feed = weights * (1 - keep) / rates
    rim = j0(BESSEL)
    state = np.zeros(len(BESSEL))
    out = np.empty(len(t))
    for i in range(len(t)):
        state = state * keep + power[i] * feed
        out[i] = state @ rim
    return out


def drift_follow(wall, start, gamma):
    out = np.empty(len(wall))
    out[0] = start
    k = math.exp(-gamma)
    for i in range(1, len(wall)):
        out[i] = wall[i] + (out[i - 1] - wall[i]) * k
    return out


def one_run(folder, run, leads, radius):
    with open(Path(folder) / "traces" / f"{run['run']}.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["t_s"]) for r in rows])
    y = np.array([float(r["T_rim_K"]) for r in rows])
    wall = np.array([float(r["T_chamber_K"]) for r in rows])
    heater = "heater_intervals_s" in run
    spans = run["heater_intervals_s"] if heater else run["shutter_intervals_s"]
    lit = np.zeros(len(t), bool)
    for a, b in spans:
        lit |= (t >= a) & (t < b)
    if heater:
        amps = run["supply_current_A"]
        power = np.where(lit, amps * (run["supply_voltage_V"] - amps * leads), 0.0)
        scatter = np.zeros(len(t))
    else:
        mon = np.array([float(r["monitor_mV"]) for r in rows])
        k = np.median([w / np.median(mon[np.abs(t - ts) <= 8]) for ts, w in run["power_meter_W"]])
        power = np.where(lit, mon * k * 0.9962, 0.0)
        scatter = np.where(lit, mon / mon.mean(), 0.0)
    weights = mode_weights(radius, not heater)
    cache = {}

    def model(x):
        per_j, log_g, start, step = x
        g = math.exp(log_g)
        if log_g not in cache:
            cache.clear()
            cache[log_g] = response(t, power, g, radius, weights)
        return drift_follow(wall, y[0] + start * 1e-3, g) + per_j * 1e-3 * cache[log_g] + step * 1e-4 * scatter

    guess = max(y.max() - np.median(y[: int(spans[0][0]) - 5]), 1e-9) / max(power.sum(), 1e-30) * 1e3
    fit = least_squares(lambda x: (model(x) - y) * 1e4, [guess, math.log(2.5e-4), 0.0, 0.0], method="lm",
                        x_scale=[0.05 * guess, 0.05, 0.01, 0.1], xtol=1e-12, ftol=1e-12, max_nfev=3000)
    return fit.x[0] * 1e-3


def main(folder, target):
    rec = json.loads((Path(folder) / "runs.json").read_text())
    spec = {s["id"]: s for s in rec["specimens"]}
    per = {r["specimen"]: one_run(folder, r, rec["heater_leads_ohm"], spec[r["specimen"]]["diameter_mm"] / 2000.0)
           for r in rec["runs"]}
    dummy = next(s for s in spec.values() if s["kind"] == "dummy")
    extra = 1.0 / per[dummy["id"]] - dummy["mass_g"] * 1e-3 * 740.0
    absorb = {sid: per[sid] * (spec[sid]["mass_g"] * 1e-3 * 740.0 + extra) for sid in per if spec[sid]["kind"] != "dummy"}
    wins = [(spec[s]["thickness_mm"] / 10.0, a) for s, a in absorb.items() if spec[s]["kind"] == "window"]
    L = np.array([w[0] for w in wins])
    A = np.array([w[1] for w in wins])
    sw = 1.0 / np.sqrt(L)
    (c0, c1), *_ = np.linalg.lstsq(np.c_[np.ones_like(L), L] * sw[:, None], A * sw, rcond=None)
    mirror = next(a for s, a in absorb.items() if spec[s]["kind"] == "mirror")
    Path(target).write_text(json.dumps({"alpha_bulk_ppm_per_cm": c1 * 1e6, "coating_absorptance_ppm_per_surface": c0 / 2 * 1e6,
                                        "mirror_absorptance_ppm": mirror * 1e6}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
