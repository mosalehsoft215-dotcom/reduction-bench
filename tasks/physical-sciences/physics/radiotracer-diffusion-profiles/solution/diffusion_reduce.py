"""Radiotracer serial sectioning -> activation energy and log10 D(1200 K).

    python3 diffusion_reduce.py <series folder> <out.json>

1. Depths from the section masses (density, area); count rates net of background, divided by the section mass and
   corrected for decay to the end of the anneal.
2. The first two sections (the surface layer, where tracer that never entered the alloy can sit) are left out, and
   each profile is fitted, in counts with Poisson weights, by S0 [<exp(-x^2 / w)> + g <exp(-(x / lam)^1.2)>] over
   the sections (section averages, the first term exactly through erf), g >= 0 and the tail length between 0.8 and 15 times the first width estimate;
   w = 4 integral D dt.
3. Arrhenius from the specimens' temperature records: w / 4 = D(1200) * integral exp(-Q/R (1/T - 1/1200)) dt; for a
   trial Q, ln D(1200) follows in closed form and Q comes from a bounded 1-D search.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import erf

RG = 8.314462618


def table(path):
    with open(path, newline="") as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def fit_profile(sec, sheet, tail=True, surface=True):
    lam = math.log(2) / sheet["half_life_d"]
    m = np.array([r["mass_removed_mg"] for r in sec])
    n = np.array([r["counts"] for r in sec])
    live = np.array([r["live_s"] for r in sec])
    decay = np.exp(-lam * np.array([r["counted_day"] for r in sec]))
    thick = m * 1e-3 / sheet["density_g_cm3"] / sheet["specimen_area_cm2"] * 1e-2          # m
    xb = np.cumsum(thick)
    xa = xb - thick
    bkg = sheet["background_cps"]
    fine = np.linspace(0, 1, 41)

    def shape(p):
        lw, lg_, llam = p[1], p[2], p[3]
        w = math.exp(lw)
        s = math.sqrt(w)
        gauss = math.sqrt(math.pi) / 2 * s * (erf(xb / s) - erf(xa / s)) / thick
        if not tail:
            return gauss
        lam_t = math.exp(llam)
        xs = xa[:, None] + thick[:, None] * fine[None, :]
        t = np.exp(-(xs / lam_t) ** 1.2)
        tavg = np.trapezoid(t, fine, axis=1)
        return gauss + math.exp(lg_) * tavg

    keep = np.arange(len(sec)) >= (2 if surface else 0)

    def resid(p):
        mu = (math.exp(p[0]) * m * decay * shape(p) + bkg) * live
        return ((n - mu) / np.sqrt(np.maximum(mu, 1.0)))[keep]

    s_surf = (n[2] / live[2] - bkg) / (m[2] * decay[2])
    half = np.argmax((n / live - bkg) / (m * decay) < 0.37 * s_surf)
    w0 = max(xb[max(half, 1)] ** 2, 1e-14)
    p0 = np.array([math.log(s_surf), math.log(w0), math.log(1e-3), math.log(2.6 * math.sqrt(w0))])
    if not tail:
        sol = least_squares(lambda p: resid(np.r_[p, p0[2:]]), p0[:2], method="lm")
        return math.exp(sol.x[1])
    r0 = math.sqrt(w0)
    lo = [-np.inf, math.log(w0 / 25.0), math.log(1e-7), math.log(0.8 * r0)]
    hi = [np.inf, math.log(w0 * 25.0), math.log(0.2), math.log(15.0 * r0)]
    best = None
    for g0 in (1e-5, 1e-3, 1e-2):
        start = p0.copy()
        start[2] = math.log(g0)
        sol = least_squares(resid, start, bounds=(lo, hi), x_scale="jac", max_nfev=4000)
        if best is None or sol.cost < best.cost:
            best = sol
    return math.exp(best.x[1])


def reduce_series(folder, records=True, tail=True, surface=True):
    folder = Path(folder)
    sheet = json.loads((folder / "series.json").read_text())
    rows = table(folder / "anneals.csv")
    w, recs = [], []
    for r in rows:
        j = int(r["anneal"])
        w.append(fit_profile(table(folder / "sections" / f"anneal-{j}.csv"), sheet, tail, surface))
        if records:
            d = np.loadtxt(folder / "records" / f"anneal-{j}.csv", delimiter=",", skiprows=1)
            recs.append((d[:, 0] * 60.0, d[:, 1]))
        else:
            recs.append((np.array([0.0, (r["removed_min"] - r["inserted_min"]) * 60.0]),
                         np.array([r["setpoint_k"], r["setpoint_k"]])))
    y = np.log(np.array(w) / 4.0)

    def level(q):
        lnI = np.array([math.log(np.trapezoid(np.exp(-q / RG * (1.0 / T - 1.0 / 1200.0)), t)) for t, T in recs])
        a = np.mean(y - lnI)
        return a, float(np.sum((y - lnI - a) ** 2))

    q = minimize_scalar(lambda q: level(q)[1], bounds=(5e4, 6e5), method="bounded", options={"xatol": 1.0}).x
    return {"activation_energy_kj_mol": q / 1000.0, "log10_diffusivity_1200k": level(q)[0] / math.log(10.0)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_series(sys.argv[1])) + "\n")
