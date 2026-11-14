"""Mean sheet resistance and its non-uniformity from one wafer's four-point-probe site log.

Run as: python3 sheet_map.py <wafer folder> <result.json>

1. Temperature.  A wafer from the coater rests on pins above the chuck and settles toward the chuck temperature with
   one time constant, so every reading carries the factor 1 + A exp(-t / tau) (t from the start of the map).  A and
   tau come from all readings at once: inside a touchdown the readings follow the wafer's cooling, and site 1's
   closing touchdown must return site 1's value (the centre is far from every edge, so the infinite-sheet dual
   relation compares the two touchdowns).  A wafer that was already at the chuck temperature gives A = 0.
2. Readings.  For each touchdown and routing the voltage is fitted against the signed current.  A constant absorbs
   the thermoelectric offset (it does not follow the current direction).  When the source stepped through several
   currents, the contact spots warm the film and the resistance rises as I^2, so V = V0 + (R I + c I^3) f(t) with
   f = 1 + A exp(-t / tau), and R is the film's resistance at the chuck temperature (zero current).
3. Geometry.  The film fills a round wafer with a primary flat, on an insulating substrate.  A Moebius map sends the
   two ends of the flat to 0 and infinity, turning the wafer into a wedge, and a power map opens the wedge to a
   half-plane, where a needle's Neumann potential is the needle plus its mirror image.  Each touchdown is solved with
   that exact potential for the sheet resistance and the head's middle-gap excess e (needles at gaps 1, 1 + e, 1 mm
   about the site, along +x): R_A / R_B fixes e, R_A then gives the sheet resistance.
4. Wafer.  A site's value is the mean of its touchdowns; the mean of the sites and 100 times their sample standard
   deviation over the mean.
"""
import cmath
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares


class Wafer:
    """Round wafer of `radius` with a flat of `flat` mm on the -y side (flat = 0: a plain disk)."""

    def __init__(self, radius, flat):
        self.a = radius
        self.flat = flat
        if flat <= 0:
            return
        c = flat / 2
        self.h = math.sqrt(radius ** 2 - c ** 2)
        self.A, self.B = complex(-c, -self.h), complex(c, -self.h)
        arc = cmath.phase(self._mobius(complex(0.0, radius))) % (2 * math.pi)
        lo, hi = sorted((arc, math.pi))
        centre = cmath.phase(self._mobius(0j)) % (2 * math.pi)
        self.start, self.alpha = (lo, hi - lo) if lo < centre < hi else (hi, 2 * math.pi - (hi - lo))

    def _mobius(self, z):
        return (z - self.A) / (z - self.B)

    def _map(self, z):
        w = self._mobius(np.asarray(z, dtype=complex))
        turn = (np.angle(w) - self.start) % (2 * np.pi)
        return np.abs(w) ** (np.pi / self.alpha) * np.exp(1j * turn * np.pi / self.alpha)

    def potential(self, point, z_in, z_out):
        """Unit current in at z_in, out at z_out, unit sheet resistance."""
        if self.flat <= 0:                         # disk: image of the same sign at a^2 / conj(z)
            img = lambda s: self.a ** 2 / np.conj(s)
            term = lambda s: math.log(abs(point - s) * abs(point - img(s)))
        else:
            p = self._map(point)
            term = lambda s: (lambda q: math.log(abs(p - q) * abs(p - np.conj(q))))(self._map(s))
        return (term(z_out) - term(z_in)) / (2 * math.pi)


def read_log(folder):
    """Readings grouped by (touchdown, routing) as rows (current A, volts, t s); touchdown -> site; site -> position."""
    groups, site_of, where = defaultdict(list), {}, {}
    with open(Path(folder) / "sites.csv", newline="") as fh:
        for n, row in enumerate(csv.DictReader(fh)):
            site = int(row["site"])
            touch = int(row.get("touchdown") or site)
            site_of[touch] = site
            where[site] = complex(float(row["x_mm"]), float(row["y_mm"]))
            t = float(row["t_s"]) if row.get("t_s") else 0.0
            groups[(touch, row["routing"])].append((float(row["current_mA"]) * 1e-3, float(row["volts"]), t))
    return {k: np.array(v) for k, v in groups.items()}, site_of, where


def fit_reading(arr, lift, tau, warming=True, lowest=False):
    """Zero-current resistance at the chuck temperature of one touchdown and routing, and the fit residuals."""
    i, v, t = arr.T
    if lowest:                                     # wrong route: only the smallest current
        keep = np.abs(i) <= np.abs(i).min() * 1.0001
        i, v, t = i[keep], v[keep], t[keep]
    f = 1 + lift * np.exp(-t / tau)
    levels = len(np.unique(np.round(np.abs(i), 9)))
    m = np.column_stack([np.ones_like(i), i * f] + ([i ** 3 * f] if warming and levels >= 2 else []))
    coef = np.linalg.lstsq(m, v, rcond=None)[0]
    return coef[1], v - m @ coef


def cooling(groups, site_of, warming=True):
    """Lift A and time constant tau of the wafer's settling, from every reading and the repeated sites."""
    repeats = defaultdict(list)
    for touch, site in site_of.items():
        repeats[site].append(touch)
    repeats = [sorted(v) for v in repeats.values() if len(v) > 1]

    def residual(p):
        lift, tau = p[0], math.exp(p[1])
        r, out = {}, []
        for key, arr in groups.items():
            r[key], res = fit_reading(arr, lift, tau, warming)
            out.append(res / 1e-6)
        for touches in repeats:
            first = dual(r[(touches[0], "A")], r[(touches[0], "B")])
            out.append(np.array([dual(r[(k, "A")], r[(k, "B")]) / first - 1 for k in touches[1:]]) / 3e-5)
        return np.concatenate(out)
    best = None
    for lift0 in (0.0, 0.05):
        for tau0 in (30.0, 100.0, 300.0):
            fit = least_squares(residual, [lift0, math.log(tau0)], x_scale=[0.01, 0.3],
                                bounds=([-0.3, math.log(3.0)], [0.6, math.log(1e4)]))
            if best is None or fit.cost < best.cost:
                best = fit
    return best.x[0], math.exp(best.x[1])


def factors(wafer, centre, excess, spacing):
    half, inner = spacing * (1.5 + excess / 2), spacing * (0.5 + excess / 2)
    z = centre + np.array([-half, -inner, inner, half])
    f = wafer.potential
    return f(z[1], z[0], z[3]) - f(z[2], z[0], z[3]), f(z[1], z[0], z[2]) - f(z[3], z[0], z[2])


def site_value(wafer, centre, r_a, r_b, spacing):
    gap = lambda e: (lambda a, b: a / b)(*factors(wafer, centre, e, spacing)) - r_a / r_b
    excess = brentq(gap, -0.4, 0.6, xtol=1e-13)
    return r_a / factors(wafer, centre, excess, spacing)[0]


def dual(r_a, r_b):
    """Rymaszewski: exp(2 pi R_A / Rs) - exp(2 pi R_B / Rs) = 1 for any needle spacings on an infinite sheet."""
    f = lambda rs: math.exp(2 * math.pi * r_a / rs) - math.exp(2 * math.pi * r_b / rs) - 1
    guess = math.pi / math.log(2) * r_a
    return brentq(f, guess * 0.5, guess * 2.0, xtol=1e-12 * guess)


def measure(folder, *, warming=True, lowest=False, flat=True, edge=True, cool="fit"):
    """cool: fit (the reference); none, start, end (wrong routes: no settling, or a straight-line drift through
    site 1's two touchdowns referred to the first or to the closing touchdown)."""
    folder = Path(folder)
    groups, site_of, where = read_log(folder)
    timed = any(arr[:, 2].max() != arr[:, 2].min() for arr in groups.values())
    lift, tau = cooling(groups, site_of, warming) if cool == "fit" and timed else (0.0, 1.0)
    r = {key: fit_reading(arr, lift, tau, warming, lowest)[0] for key, arr in groups.items()}
    meta = json.loads((folder / "wafer.json").read_text())
    touches = sorted(site_of)
    if edge:
        wafer = Wafer(meta["diameter_mm"] / 2, meta.get("primary_flat_mm", 0.0) if flat else 0.0)
        value = {k: site_value(wafer, where[site_of[k]], r[(k, "A")], r[(k, "B")], meta["nominal_spacing_mm"])
                 for k in touches}
    else:
        value = {k: dual(r[(k, "A")], r[(k, "B")]) for k in touches}
    if cool in ("start", "end"):
        when = {k: np.mean([groups[(k, q)][:, 2].mean() for q in ("A", "B")]) for k in touches}
        pair = [k for k in touches if site_of[k] == site_of[touches[0]]]
        slope = (value[pair[-1]] / value[pair[0]] - 1) / (when[pair[-1]] - when[pair[0]])
        ref = when[pair[0]] if cool == "start" else when[pair[-1]]
        value = {k: v / (1 + slope * (when[k] - ref)) for k, v in value.items()}
    sites = sorted(set(site_of.values()))
    values = np.array([np.mean([value[k] for k in touches if site_of[k] == s]) for s in sites])
    return {"sheet_resistance_ohm_sq": float(values.mean()),
            "nonuniformity_percent": float(100 * values.std(ddof=1) / values.mean())}


if __name__ == "__main__":
    result = measure(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result, allow_nan=False) + "\n", encoding="utf-8")
