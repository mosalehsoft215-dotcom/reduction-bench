"""Mesopore volume and median pore diameter of one mercury-intrusion fill.

    python3 porosimetry_reduce.py <fill folder> <answer.json>

1. Mercury in the fill from the weighings; the blank run gives the penetrometer's own apparent intrusion, of which
   the mercury compression (4.05e-5 of the mercury volume per MPa) is rescaled from the blank's mercury to the fill's.
2. Each step's volume is carried to equilibrium: within a hold the volume approaches its end exponentially, so the two
   readings and the rate give the time constant ((v10 - v5) / rate = tau (exp(5 / tau) - 1)) and the volume still to
   come (rate x tau); tau is the fill's rate-weighted median over the steps where the rate stands clear of the noise.
3. Low-pressure port: the sample feels the logged pressure plus rho g h of the mercury above it, h falling as the stem
   empties.  Washburn D = -4 gamma cos(theta) / P.
4. Skeleton compression: above 300 MPa (pores narrower than 5 nm) no pore volume is left, so the volume there rises only as the
   skeleton compresses, in proportion to the pressure; that slope times the pressure is taken off every step.
5. The inter-extrudate voids end at the minimum of the log-differential intrusion between 3 and 40 um; mesopore volume
   = intrusion between 50 nm and 4 nm per gram; median = the diameter at half of the intraparticle intrusion down to
   4 nm (log interpolation).
"""
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

KAPPA = 4.05e-5
G = 9.81


def table(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def tau_of(v5, v10, rate):
    r = (v10 - v5) / rate
    if not 5.02 < r < 1e3:                 # outside what an exponential approach within the hold can give
        return float("nan")
    return brentq(lambda t: t * (np.exp(5.0 / t) - 1.0) - r, 0.05, 1e4)


def reduce_fill(folder, equil=True, soft=True, scaled_blank=True, head=True):
    folder = Path(folder)
    sh = json.loads((folder / "fill_sheet.json").read_text())
    rows = table(folder / "intrusion.csv")
    port = np.array([r["port"] for r in rows])
    p = np.array([float(r["pressure_mpa"]) for r in rows])
    v5 = np.array([float(r["volume_ml_at_5s"]) for r in rows])
    v10 = np.array([float(r["volume_ml_at_10s"]) for r in rows])
    rate = np.array([float(r["rate_ml_per_s_at_10s"]) for r in rows])
    blank = np.array([[float(r["pressure_mpa"]), float(r["apparent_volume_ml"])] for r in table(folder / "blank_run.csv")])
    hg = (sh["filled_mass_g"] - sh["assembly_mass_g"]) / sh["mercury_density_g_per_ml"]
    b = np.interp(np.log(p), np.log(blank[:, 0]), blank[:, 1])
    if scaled_blank:
        b = b + KAPPA * p * (hg - sh["blank_run_mercury_ml"])
    v = v10 - b
    if equil:
        use = rate > 3e-4
        taus = np.array([tau_of(a, c, r) for a, c, r in zip(v5[use], v10[use], rate[use])])
        ok = np.isfinite(taus)
        taus, wts = taus[ok], rate[use][ok]
        tau = 0.0
        if taus.size:
            order = np.argsort(taus)
            cw = np.cumsum(wts[order])
            tau = float(taus[order][np.searchsorted(cw, cw[-1] / 2)])
        v = v + rate * tau
    if soft:
        top = p > 300.0
        slope = np.polyfit(p[top], v[top], 1)[0]
        v = v - slope * p
    ps = p.copy()
    if head:
        low = port == "low"
        h = np.clip(sh["head_at_full_stem_cm"] - v[low] / sh["stem_bore_cm2"], 0.0, None)
        ps[low] = p[low] + sh["mercury_density_g_per_ml"] * 1e3 * G * h / 100 * 1e-6
    d = -4 * sh["mercury_surface_tension_n_per_m"] * np.cos(np.deg2rad(sh["contact_angle_deg"])) / (ps * 1e6) * 1e9
    order = np.argsort(-d)
    d, v = d[order], np.maximum.accumulate(v[order])
    logd = np.log(d)
    grid = np.linspace(np.log(40000.0), np.log(3000.0), 200)
    cum = np.interp(-grid, -logd, v)
    dens = np.gradient(cum, -grid)
    cut = float(np.exp(grid[int(np.argmin(dens))]))
    at = lambda dd: float(np.interp(-np.log(dd), -logd, v))  # noqa: E731
    base, end = at(cut), at(4.0)
    meso = (at(4.0) - at(50.0)) / sh["sample_mass_g"]
    half = base + 0.5 * (end - base)
    median = float(np.exp(np.interp(half, v, logd))) if np.all(np.diff(v) >= 0) else float("nan")
    return {"mesopore_volume_ml_per_g": float(meso), "median_pore_diameter_nm": median}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_fill(sys.argv[1])) + "\n")
