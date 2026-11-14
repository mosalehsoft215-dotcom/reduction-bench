"""Independent reduction of one wafer log, written apart from solution/sheet_map.py.

Settling: a warm wafer cools toward the chuck as 1 + A exp(-t / tau).  Each touchdown's readings are fitted with a
local drift term D I (t - tmean); the relative drift rates of all touchdowns give the decay, and site 1's first and
closing touchdowns (geometry solved) anchor the lift.  Every reading is then divided by its factor.
Readings: for each touchdown, routing and current level, half the difference of the two current directions over the
current (the thermoelectric offset cancels); the warm-up of the film then follows a straight line in I^2, whose value
at I = 0 is the resistance at the chuck temperature.  Geometry: the method of fundamental solutions.  The wafer's
boundary (arc plus primary flat) is sampled densely; the potential of the two current needles is the free-sheet
logarithm plus a field of sources placed just outside the wafer, whose strengths make the normal current vanish at
every boundary sample (least squares).  The sheet resistance and the head's middle-gap excess come from the two
routings together by a two-dimensional Newton solve (scipy fsolve).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve, least_squares


def load(folder):
    raw = np.genfromtxt(Path(folder) / "sites.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    names = raw.dtype.names
    touch = raw["touchdown"] if "touchdown" in names else raw["site"]
    t = raw["t_s"].astype(float) if "t_s" in names else np.arange(raw.size, dtype=float)
    return raw, touch, t


def drift_rates(raw, touch, t):
    """Relative resistance drift (1/s) and mean time of every touchdown."""
    out = []
    for k in np.unique(touch):
        rates, weights = [], []
        for routing in ("A", "B"):
            pick = (touch == k) & (raw["routing"] == routing)
            i, v, tt = raw["current_mA"][pick] * 1e-3, raw["volts"][pick].astype(float), t[pick]
            if np.ptp(tt) == 0 or i.size < 4:
                continue
            cols = [np.ones_like(i), i, i * (tt - tt.mean())]
            if len(set(np.round(np.abs(i), 9))) >= 2:
                cols.append(i ** 3)
            m = np.column_stack(cols)
            coef = np.linalg.lstsq(m, v, rcond=None)[0]
            rates.append(coef[2] / coef[1])
        if rates:
            out.append((float(t[touch == k].mean()), float(np.mean(rates))))
    return np.array(out)


def readings(folder, lift=0.0, tau=1.0):
    raw, touch, t = load(folder)
    out, where, site_of = {}, {}, {}
    for k in np.unique(touch):
        pick_t = touch == k
        site = int(raw["site"][pick_t][0])
        site_of[int(k)] = site
        where[site] = complex(float(raw["x_mm"][pick_t][0]), float(raw["y_mm"][pick_t][0]))
        for routing in ("A", "B"):
            pick = pick_t & (raw["routing"] == routing)
            i = raw["current_mA"][pick] * 1e-3
            v = raw["volts"][pick].astype(float) / (1 + lift * np.exp(-t[pick] / tau))
            levels = sorted(set(np.round(np.abs(i), 9)))
            per = []
            for lev in levels:
                up = v[np.isclose(i, lev)].mean()
                down = v[np.isclose(i, -lev)].mean()
                per.append((lev, (up - down) / (2 * lev)))
            per = np.array(per)
            if len(per) == 1:
                out[(int(k), routing)] = per[0, 1]
            else:
                out[(int(k), routing)] = np.polyfit(per[:, 0] ** 2, per[:, 1], 1)[1]
    return out, where, site_of


class Boundary:
    def __init__(self, radius, flat, points=900, lift=1.2):
        c = flat / 2
        h = math.sqrt(radius ** 2 - c ** 2)
        t0 = math.atan2(-h, c)                       # the arc runs from the right end of the flat round to the left
        arc_t = np.linspace(t0, t0 + 2 * math.pi - 2 * math.asin(c / radius), int(points * 0.8))
        arc = radius * np.exp(1j * arc_t)
        line = np.linspace(-c, c, int(points * 0.2)) - 1j * h
        self.z = np.concatenate([arc, line])
        self.n = np.concatenate([np.exp(1j * arc_t), np.full(line.size, -1j)])
        self.sources = self.z + lift * self.n        # just outside the wafer
        d = self.z[:, None] - self.sources[None, :]
        # the source-to-boundary normal-gradient matrix depends only on the wafer: invert it once
        self.solver = np.linalg.pinv(np.real(np.conj(d) * self.n[:, None]) / np.abs(d) ** 2, rcond=1e-12)

    def correction(self, z_in, z_out):
        """Strengths of the outside sources that cancel the free field's normal current at the boundary."""
        def normal_grad(z, s):                       # d/dn of ln|z - s|
            d = z - s
            return np.real(np.conj(d) * self.n[:, None] if np.ndim(s) else np.conj(d) * self.n) / np.abs(d) ** 2
        free = normal_grad(self.z, z_out) - normal_grad(self.z, z_in)
        return -self.solver @ free

    def potential(self, point, z_in, z_out, strengths):
        return (math.log(abs(point - z_out)) - math.log(abs(point - z_in))
                + float(strengths @ np.log(np.abs(point - self.sources)))) / (2 * math.pi)


def factors(bound, centre, excess, spacing):
    z = centre + spacing * np.array([-1.5, -0.5, 0.5, 1.5]) + spacing * excess * np.array([-0.5, -0.5, 0.5, 0.5])
    qa = bound.correction(z[0], z[3])
    qb = bound.correction(z[0], z[2])
    fa = bound.potential(z[1], z[0], z[3], qa) - bound.potential(z[2], z[0], z[3], qa)
    fb = bound.potential(z[1], z[0], z[2], qb) - bound.potential(z[3], z[0], z[2], qb)
    return fa, fb


def solve_site(bound, centre, r_a, r_b, spacing):
    start = math.pi / math.log(2) * r_a

    def equations(p):
        fa, fb = factors(bound, centre, p[1], spacing)
        return [p[0] * fa / r_a - 1, p[0] * fb / r_b - 1]
    return fsolve(equations, [start, 0.0], xtol=1e-12)[0]


def reduce(folder):
    folder = Path(folder)
    meta = json.loads((folder / "wafer.json").read_text())
    bound = Boundary(meta["diameter_mm"] / 2, meta.get("primary_flat_mm", 0.0))
    spacing = meta["nominal_spacing_mm"]
    raw, touch, t = load(folder)
    rates = drift_rates(raw, touch, t)
    lift, tau = 0.0, 1.0
    if rates.size and np.abs(rates[:, 1]).max() > 2e-5:
        r0, where, site_of = readings(folder)
        first = [k for k in sorted(site_of) if site_of[k] == 1]
        t1, t2 = (float(t[touch == k].mean()) for k in (first[0], first[-1]))
        target = 1.0                                   # corrected site-1 ratio still to be explained
        for _ in range(4):                             # the anchor is re-read through the current correction
            rc, _, _ = readings(folder, lift, tau)
            val = [solve_site(bound, where[1], rc[(k, "A")], rc[(k, "B")], spacing) for k in (first[0], first[-1])]
            ratio = val[0] / val[1] * (1 + lift * math.exp(-t1 / tau)) / (1 + lift * math.exp(-t2 / tau))

            def residual(p):
                a, tau_ = p[0], math.exp(p[1])
                e = a * np.exp(-rates[:, 0] / tau_)
                model = -(e / tau_) / (1 + e)
                anchor = (1 + a * math.exp(-t1 / tau_)) / (1 + a * math.exp(-t2 / tau_)) - ratio
                return np.concatenate([(model - rates[:, 1]) / 2e-5, [anchor / 3e-5]])
            fits = [least_squares(residual, [ratio - 1, math.log(tau0)], bounds=([-0.3, 1.0], [0.6, 9.0]))
                    for tau0 in (40.0, 120.0, 400.0)]
            best = min(fits, key=lambda f: f.cost)
            lift, tau = best.x[0], math.exp(best.x[1])
    r, where, site_of = readings(folder, lift, tau)
    per_touch = {k: solve_site(bound, where[site_of[k]], r[(k, "A")], r[(k, "B")], spacing) for k in sorted(site_of)}
    values = np.array([np.mean([v for k, v in per_touch.items() if site_of[k] == s]) for s in sorted(where)])
    mean = values.mean()
    spread = math.sqrt(((values - mean) ** 2).sum() / (values.size - 1))
    return {"sheet_resistance_ohm_sq": float(mean), "nonuniformity_percent": float(100 * spread / mean)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + chr(10), encoding="utf-8")
