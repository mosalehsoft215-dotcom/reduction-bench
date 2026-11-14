"""Reduces one 22.2 GHz calibration track to the aperture efficiency and the source flux density.

    python3 track_reduce.py <track directory> <result.json>

1. Loads: G = (hot - cold counts) / (T_ambient - 80.0 K) and the diode's temperature from its deflection on both loads.
2. Tipping: the sky counts over G against airmass A = 1/sin(el), fitted as c - (0.95 T_ambient - 2.73) e^(-tau A).
   At every scan the off-source system temperature on the diode scale, put into the same relation with the tipping's
   c, gives the opacity at the time of that scan.
3. Scans: each scan is put on the kelvin scale by its own diode deflection (on against the mean of the off
   blocks before and after it), a linear baseline from its outer parts is
   removed, and a Gaussian is fitted near the peak; the peak is raised for the pointing offset that the scans along the
   other axis show, then corrected for e^(-tau A) and the gain curve at its elevation.
4. Planet: Rayleigh-Jeans flux of the disc with Omega = pi R_eq R_pol / d^2 and the disc dilution
   (1 - e^-x) / x, x = ln2 theta^2 / FWHM^2 with theta the equal-area diameter; eta = 2 k T* / (A_geo S K).
5. Source: S = 2 k T* (W_az W_el / FWHM^2) / (eta A_geo), the total flux of a source whose fitted widths W exceed
   the beam.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

C, KB, AU = 299792458.0, 1.380649e-23, 1.495978707e11
DISH, T_COLD, T_BG = 32.0, 80.0, 2.73
EL0, GA = 48.0, -7.5e-5
PLANETS = {"Jupiter": (71492e3, 66854e3, 136.0), "Venus": (6100e3, 6100e3, 480.0), "Mars": (3396.2e3, 3376.2e3, 205.0)}
L2 = math.log(2)


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def gain_curve(el):
    return 1 + GA * (el - EL0) ** 2


def load_scale(track, amb):
    lo = rows(track / "loads.csv")
    m = {(r["load"], int(r["diode"])): [] for r in lo}
    for r in lo:
        m[(r["load"], int(r["diode"]))].append(float(r["counts"]))
    m = {k: np.mean(v) for k, v in m.items()}
    g = (m[("hot", 0)] - m[("cold", 0)]) / (amb - T_COLD)
    tcal = 0.5 * ((m[("hot", 1)] - m[("hot", 0)]) + (m[("cold", 1)] - m[("cold", 0)])) / g
    return g, tcal


def tipping_tau(track, g, amb):
    tp = rows(track / "tipping.csv")
    els = sorted({float(r["elevation_deg"]) for r in tp})
    y = np.array([np.mean([float(r["counts"]) for r in tp if float(r["elevation_deg"]) == e]) / g for e in els])
    air = 1 / np.sin(np.radians(els))
    amp = 0.95 * amb - T_BG
    (c0, tau), _ = curve_fit(lambda A, c, t: c - amp * np.exp(-t * A), air, y, p0=(y.max(), 0.1))
    return tau, c0


def scan_peaks(track, tcal, fwhm_as, c0, amp):
    sc = rows(track / "scans.csv")
    out = []
    for k in sorted({int(r["scan"]) for r in sc}):
        rr = [r for r in sc if int(r["scan"]) == k]
        c = np.array([float(r["counts"]) for r in rr])
        d = np.array([int(r["diode"]) for r in rr])
        n_on = np.flatnonzero(d == 1)
        before, after = c[n_on[0] - 20:n_on[0]], c[n_on[-1] + 1:n_on[-1] + 21]
        k_per_count = tcal / (c[n_on].mean() - 0.5 * (before.mean() + after.mean()))
        # the off-source system temperature at the scan start gives the opacity at the time of the scan
        def opacity(block):
            el_b = np.mean([float(rr[i]["elevation_deg"]) for i in block])
            t_sys = c[block].mean() * k_per_count
            return -math.log((c0 - t_sys) / amp) * math.sin(math.radians(el_b))
        head = np.arange(n_on[-1] + 1, n_on[-1] + 21)
        tail = np.arange(len(rr) - 20, len(rr))
        times = np.array([float(r["time_s"]) for r in rr])
        t_head, t_tail = times[head].mean(), times[tail].mean()
        tau_head, tau_tail = opacity(head), opacity(tail)
        body = slice(n_on[-1] + 21, None)
        x = np.array([float(r["offset_arcsec"]) for r in rr])[body]
        el = np.array([float(r["elevation_deg"]) for r in rr])[body]
        y = c[body]
        outer = np.abs(x) > 2.5 * fwhm_as
        base = np.polyval(np.polyfit(x[outer], y[outer], 1), x)
        temp = (y - base) * k_per_count
        i0 = int(np.argmax(temp))
        near = np.abs(x - x[i0]) < (0.9 if rr[0]["target"] == "planet" else 2.2) * fwhm_as
        gauss = lambda u, a, u0, w: a * np.exp(-4 * L2 * (u - u0) ** 2 / w ** 2)
        (a, u0, w), _ = curve_fit(gauss, x[near], temp[near], p0=(temp[i0], x[i0], fwhm_as))
        out.append({"target": rr[0]["target"], "axis": rr[0]["axis"], "peak": a, "centre": u0, "width": abs(w),
                    "el": float(np.interp(u0, x, el)),
                    "tau": float(np.interp(np.interp(u0, x, times[body]), [t_head, t_tail], [tau_head, tau_tail]))})
    return out


def corrected(peaks, target):
    p = [q for q in peaks if q["target"] == target]
    dx = np.mean([q["centre"] for q in p if q["axis"] == "az"])
    dy = np.mean([q["centre"] for q in p if q["axis"] == "el"])
    vals = []
    for q in p:
        other = dy if q["axis"] == "az" else dx
        w = np.mean([s["width"] for s in p if s["axis"] != q["axis"]])
        up = q["peak"] * math.exp(4 * L2 * other ** 2 / w ** 2)
        vals.append(up * math.exp(q["tau"] / math.sin(math.radians(q["el"]))) / gain_curve(q["el"]))
    return float(np.mean(vals))


def reduce_track(track):
    track = Path(track)
    head = json.loads((track / "track.json").read_text())
    nu = head["frequency_ghz"] * 1e9
    fwhm = 1.16 * C / nu / DISH
    fwhm_as = math.degrees(fwhm) * 3600
    amb = head["ambient_k"]
    g, tcal = load_scale(track, amb)
    tau, c0 = tipping_tau(track, g, amb)
    peaks = scan_peaks(track, tcal, fwhm_as, c0, 0.95 * amb - T_BG)
    rp, rq, tb = PLANETS[head["planet"]]
    d = head["planet_distance_au"] * AU
    s_planet = 2 * KB * nu ** 2 / C ** 2 * tb * math.pi * rp * rq / d ** 2
    x = L2 * (2 * math.sqrt(rp * rq) / d) ** 2 / fwhm ** 2
    dilution = (1 - math.exp(-x)) / x
    a_geo = math.pi * (DISH / 2) ** 2
    eta = 2 * KB * corrected(peaks, "planet") / (a_geo * s_planet * dilution)
    # total flux of a Gaussian source: peak times (observed width / beam width) along both axes
    src = [q for q in peaks if q["target"] == "source"]
    spread = np.prod([np.mean([q["width"] for q in src if q["axis"] == ax]) / fwhm_as for ax in ("az", "el")])
    if spread < 1.03 ** 2:
        spread = 1.0  # widths at the beam within their scatter: an unresolved source
    flux = 2 * KB * corrected(peaks, "source") * spread / (eta * a_geo) / 1e-26
    return {"aperture_efficiency": eta, "source_flux_jy": flux}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_track(sys.argv[1])))
