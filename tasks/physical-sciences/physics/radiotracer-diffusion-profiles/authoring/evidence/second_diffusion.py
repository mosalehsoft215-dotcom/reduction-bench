"""Second, separately written reduction of a radiotracer diffusion series (tail peeling and effective times).

Profile: the two surface sections are set aside (tracer can sit there without having entered the alloy); specific
activities at the depths where the current Gaussian equals its section average; the deep tail is fitted first as a straight line of ln S against
x^(6/5) over the deepest third of the sections (after removing the Gaussian of the previous pass), its extrapolation
is subtracted from the rest, and a weighted straight line ln S against x^2 over the shallower sections whose remainder
is still above 8 % of the third section's value gives the width w; the two fits alternate eight times.  When the deep third holds no significant tail (its summed activity is within eight standard errors of zero or
within 1.5 times the Gaussian of a first pass) nothing is subtracted.
Arrhenius: for a trial Q every anneal gets an effective time at its setpoint,
t_eff = integral exp(-Q/R (1/T - 1/Ts)) dt over its record, and a weighted straight line ln(w / (4 t_eff)) against
1/Ts gives the next Q until it settles.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import erf

R = 8.314462618


def width(folder, j, s):
    d = np.genfromtxt(folder / "sections" / f"anneal-{j}.csv", delimiter=",", names=True)
    lam = math.log(2) / s["half_life_d"]
    dz = d["mass_removed_mg"] / (1e3 * s["density_g_cm3"] * s["specimen_area_cm2"]) * 1e4          # um
    hi = np.cumsum(dz)
    lo_ = hi - dz
    x = hi - dz / 2
    x2 = x ** 2                      # depth^2 at which the Gaussian equals its section average (updated below)
    net = d["counts"] / d["live_s"] - s["background_cps"]
    sa = net / d["mass_removed_mg"] * np.exp(lam * d["counted_day"])
    err = np.sqrt(d["counts"]) / d["live_s"] / d["mass_removed_mg"] * np.exp(lam * d["counted_day"])

    def gauss(sel, y):
        w8 = (y[sel] / err[sel]) ** 2
        A = np.vstack([np.ones(sel.sum()), x2[sel]]).T * np.sqrt(w8)[:, None]
        c = np.linalg.lstsq(A, np.log(y[sel]) * np.sqrt(w8), rcond=None)[0]
        return c

    def eff(c):
        r = math.sqrt(-1.0 / c[1])
        avg = math.sqrt(math.pi) / 2 * r * (erf(hi / r) - erf(lo_ / r)) / dz
        return -(r ** 2) * np.log(avg)

    inner = np.arange(len(x)) >= 2                 # the two surface sections are not used
    top = sa[2]
    first = gauss((sa > 0.25 * top) & inner, sa)
    deep = np.arange(len(x)) >= 2 * len(x) // 3
    pred = np.exp(first[0] + first[1] * x2[deep])
    c = first
    signif = np.sum(sa[deep]) / math.sqrt(np.sum(err[deep] ** 2))
    if signif > 8.0 and np.sum(sa[deep]) > 1.5 * np.sum(pred):
        for _ in range(8):                      # alternate: tail on the deep third, Gaussian on the remainder
            tail_data = sa - np.exp(c[0] + c[1] * x2)
            tl = deep & (tail_data > 0)
            w8 = (tail_data[tl] / err[tl]) ** 2
            A = np.vstack([np.ones(tl.sum()), x[tl] ** 1.2]).T * np.sqrt(w8)[:, None]
            t = np.linalg.lstsq(A, np.log(tail_data[tl]) * np.sqrt(w8), rcond=None)[0]
            rem = sa - np.exp(t[0] + t[1] * x ** 1.2)
            sel = (rem > 0.08 * top) & ~deep & inner
            c = gauss(sel, np.where(rem > 0, rem, 1.0))
            x2 = eff(c)
    else:
        for _ in range(4):
            c = gauss((sa > 0.02 * top) & ~deep & inner, sa)
            x2 = eff(c)
    return -1e-12 / c[1]                                                             # m^2


def reduce_series(folder):
    folder = Path(folder)
    s = json.loads((folder / "series.json").read_text())
    an = np.genfromtxt(folder / "anneals.csv", delimiter=",", names=True)
    ts = an["setpoint_k"]
    w = np.array([width(folder, int(j), s) for j in an["anneal"]])
    recs = [np.genfromtxt(folder / "records" / f"anneal-{int(j)}.csv", delimiter=",", skip_header=1)
            for j in an["anneal"]]
    q = 2.0e5
    for _ in range(80):
        teff = np.array([np.trapezoid(np.exp(-q / R * (1 / r[:, 1] - 1 / t)), r[:, 0] * 60.0) for r, t in zip(recs, ts)])
        y = np.log(w / (4 * teff))
        x = 1.0 / ts
        slope, icpt = np.polyfit(x, y, 1)
        if abs(-slope * R - q) < 1e-3:
            break
        q = -slope * R
    lnd = icpt + slope / 1200.0
    return {"activation_energy_kj_mol": q / 1000.0, "log10_diffusivity_1200k": lnd / math.log(10.0)}


if __name__ == "__main__":
    import sys
    Path(sys.argv[2]).write_text(json.dumps(reduce_series(sys.argv[1])) + "\n")
