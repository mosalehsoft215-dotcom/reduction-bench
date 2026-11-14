"""Second, separately written reduction of one mercury-intrusion fill (same call as the reference).

    python second_porosimetry.py <fill folder> <result.json>

- Penetrometer: the blank run interpolated linearly in pressure (not in log pressure), its mercury compression moved
  from the blank's mercury volume to the fill's (weighings over the mercury density).
- Equilibration: one time constant for the fill by least squares over every step whose rate exceeds 2e-4 mL/s,
  minimising the misfit of v10 - v5 = rate tau (exp(5 / tau) - 1) (bounded scalar search), instead of the
  reference's rate-weighted median of per-step roots; each step gains rate x tau.
- Skeleton: the Theil-Sen slope of volume against pressure above 300 MPa, times the pressure, taken off.
- Low-pressure head as in the notes; Washburn diameters.
- The voids end where the three-point slope of volume against ln D is smallest between 3 and 40 um (on the measured
  steps, no resampling); mesopore volume and median by linear interpolation in ln D.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar


def run(folder):
    folder = Path(folder)
    sh = json.loads((folder / "fill_sheet.json").read_text())
    raw = np.genfromtxt(folder / "intrusion.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    bl = np.genfromtxt(folder / "blank_run.csv", delimiter=",", names=True)
    p = raw["pressure_mpa"].astype(float)
    a5, a10, r = raw["volume_ml_at_5s"].astype(float), raw["volume_ml_at_10s"].astype(float), raw["rate_ml_per_s_at_10s"]
    mercury = (sh["filled_mass_g"] - sh["assembly_mass_g"]) / sh["mercury_density_g_per_ml"]
    blank = np.interp(p, bl["pressure_mpa"], bl["apparent_volume_ml"]) + 4.05e-5 * p * (mercury - sh["blank_run_mercury_ml"])
    vol = a10 - blank
    sel = r > 2e-4
    if sel.sum() >= 3:
        cost = lambda t: float(np.sum(((a10 - a5)[sel] - r[sel] * t * (np.exp(5.0 / t) - 1.0)) ** 2))  # noqa: E731
        tau = minimize_scalar(cost, bounds=(0.1, 300.0), method="bounded", options={"xatol": 1e-4}).x
        vol = vol + r * tau
    hi = np.flatnonzero(p > 300.0)
    slopes = [(vol[j] - vol[i]) / (p[j] - p[i]) for k, i in enumerate(hi) for j in hi[k + 1:]]
    vol = vol - float(np.median(slopes)) * p
    low = raw["port"] == "low"
    h = np.clip(sh["head_at_full_stem_cm"] - vol[low] / sh["stem_bore_cm2"], 0.0, None)
    pe = p.copy()
    pe[low] = p[low] + sh["mercury_density_g_per_ml"] * 9.81e-5 * h
    dia = -4 * sh["mercury_surface_tension_n_per_m"] * np.cos(np.radians(sh["contact_angle_deg"])) / pe * 1e3   # nm
    order = np.argsort(-dia)
    ln, cum = np.log(dia[order]), np.maximum.accumulate(vol[order])
    win = np.flatnonzero((dia[order] < 40000.0) & (dia[order] > 3000.0))
    slope = [abs((cum[i + 1] - cum[i - 1]) / (ln[i + 1] - ln[i - 1])) for i in win if 0 < i < ln.size - 1]
    edge = win[int(np.argmin(slope))]
    at = lambda d: float(np.interp(-np.log(d), -ln, cum))  # noqa: E731
    start, stop = cum[edge], at(4.0)
    median = float(np.exp(np.interp(start + 0.5 * (stop - start), cum, ln)))
    return {"mesopore_volume_ml_per_g": (at(4.0) - at(50.0)) / sh["sample_mass_g"], "median_pore_diameter_nm": median}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(run(sys.argv[1])))
