#!/usr/bin/env python3
"""Governor droop and deadband from one machine's primary-frequency-response trend.

usage: python3 governor_response.py <machine folder> <result file>

The governor answers x = f - f0 + g u (f the real frequency, u what the test set sends, g the gain of the
input it is wired to) with nothing inside +-db and -K (|x| - db) sign(x) outside, K = P_rated / (droop f0).

Signal choices
  frequency  the governor's own shaft-speed record is the real frequency (the shaft runs in step with the
             grid). A DCS copy of the speed can carry its loop's zero error, so the bus transducer is used
             instead; with that transducer out, the copy is used and its zero is a fitted unknown.
  power      mean of the two generator transducers, rescaled by the revenue meter. At the generator
             terminals the scale is the ratio of energies. On the HV side the meter misses the auxiliary
             load and transformer losses, a near-constant deficit, so the scale is the slope of the register
             on the transducers' running energy with that deficit fitted beside it.
  samples    a step (or rest) during which the exhaust comes within 1 degC of its control limit is left out
             with the 30 s after it, and so is the first minute.

Model: P = lag(load signal) - K lag(band(x; db)) + c - M df/dt, a one-pole lag; with the plant
controller's set point logged in place of the reference, two more columns carry its integral trim
(lagged running sum of power, lagged time). The lag and db start from a grid with the linear part solved
exactly, then everything is refined by least squares. A response limiter (a bound on the demand) and the
input gain g are each added only when they cut the misfit by more than 25 noise variances.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lfilter, lfiltic

ACCEPT = 25.0


class Trend:
    def __init__(self, folder):
        folder = Path(folder)
        self.sheet = json.loads((folder / "machine.json").read_text())
        self.setup = self.sheet.get("setup", {})
        with open(folder / "trend_1s.csv", newline="") as fh:
            rows = list(csv.DictReader(fh))
        self.col = {k: np.array([float(r[k]) if r[k].strip() else np.nan for r in rows]) for k in rows[0]}
        self.n = len(rows)

    def frequency(self):
        """(frequency in Hz, True when the copy's zero has still to be found)."""
        hz_per_rpm = self.sheet.get("poles", 2) / 120.0
        shaft = self.col["shaft_rpm"] * hz_per_rpm
        if self.setup.get("shaft_rpm_from") == "governor":
            return shaft, False
        bus = self.col["bus_hz"]
        if self.setup.get("bus_frequency") == "recorded" and not np.isnan(bus).any():
            return bus, False
        return shaft, True

    def power(self):
        mean = 0.5 * (self.col["gen_mw_1"] + self.col["gen_mw_2"])
        reg = self.col["meter_mwh"]
        if self.setup.get("revenue_meter") == "generator terminals":
            return mean * (reg[-1] - reg[0]) * 3600.0 / mean[1:].sum()
        running = np.r_[0.0, np.cumsum(mean[1:])] / 3600.0
        hours = np.arange(self.n) / 3600.0
        design = np.column_stack([running, hours, np.ones(self.n)])
        return mean * np.linalg.lstsq(design, reg, rcond=None)[0][0]

    def free(self):
        near = self.col["exhaust_degc"] >= self.col["exhaust_limit_degc"] - 1.0
        cuts = np.flatnonzero(np.diff(self.col["injected_mhz"])) + 1
        ok = np.ones(self.n, bool)
        for seg in np.split(np.arange(self.n), cuts):
            if near[seg].any():
                ok[seg[0]:seg[-1] + 31] = False
        ok[:60] = False
        return ok


def one_pole(u, tau):
    a = 1.0 / tau
    num, den = [a], [1.0, a - 1.0]
    out, _ = lfilter(num, den, u, zi=lfiltic(num, den, [u[0]], [u[0]]))
    return out


def band(x, db):
    return np.sign(x) * np.maximum(np.abs(x) - db, 0.0)


class Fit:
    """Parameter vector: [tau, db, K, (zero), (gain), (limit), c, M, (trim terms)]."""

    def __init__(self, f, u, p, load, free, f0, trimmed):
        self.f, self.u, self.p, self.load, self.free, self.f0 = f, u, p, load, free, f0
        self.rate = np.r_[0.0, np.diff(f)]
        self.trim = []
        if trimmed:
            self.trim = [np.r_[0.0, np.cumsum(p[:-1])] / 1000.0, np.arange(p.size) / 1000.0]

    def columns(self, tau, db, zero=0.0, gain=1.0, limit=None, k=1.0):
        d = band(self.f - self.f0 - zero + gain * self.u, db)
        if limit is not None:
            d = np.clip(d, -limit / k, limit / k)
        a = np.column_stack([-one_pole(d, tau), np.ones_like(d), -self.rate] + [one_pole(t, tau) for t in self.trim])
        return a, self.p - one_pole(self.load, tau)

    def exact(self, tau, db, zero=0.0):
        a, y = self.columns(tau, db, zero)
        a, y = a[self.free], y[self.free]
        coef, *_ = np.linalg.lstsq(a, y, rcond=None)
        r = a @ coef - y
        return float(r @ r), coef

    def unpack(self, v, has):
        tau, db, k = v[:3]
        i, extra = 3, {}
        for name in ("zero", "gain", "limit"):
            if has[name]:
                extra[name] = v[i]
                i += 1
        return tau, db, k, extra, v[i:]

    def misfit(self, v, has):
        tau, db, k, extra, rest = self.unpack(v, has)
        a, y = self.columns(tau, db, extra.get("zero", 0.0), extra.get("gain", 1.0), extra.get("limit"), k)
        return (a @ np.r_[k, rest] - y)[self.free]

    def refine(self, v, has):
        lo, hi = [1.0, 0.0, 1.0], [30.0, 0.1, 1e4]
        for name, l, h in (("zero", -0.05, 0.05), ("gain", 0.5, 1.5), ("limit", 0.1, 1e3)):
            if has[name]:
                lo.append(l)
                hi.append(h)
        lo += [-np.inf] * (len(v) - len(lo))
        hi += [np.inf] * (len(v) - len(hi))
        v = np.minimum(np.maximum(v, np.array(lo) + 1e-9), np.array(hi) - 1e-9)
        sol = least_squares(self.misfit, v, args=(has,), bounds=(lo, hi), x_scale="jac", diff_step=1e-6)
        return sol.x, float(sol.fun @ sol.fun)


def settle(fit, zero_unknown, look_for_limit=True, look_for_gain=True):
    start = min((fit.exact(t, d)[0], t, d) for t in np.arange(2.0, 14.01, 0.5) for d in np.arange(0.0, 0.0601, 0.001))
    _, tau, db = start
    zero = 0.0
    if zero_unknown:                              # steps fix db + zero; the far edge parts them
        _, db, zero = min((fit.exact(tau, d, z)[0], d, z) for z in np.arange(-0.03, 0.03001, 0.0005)
                          for d in np.arange(0.0, 0.0601, 0.001))
    coef = fit.exact(tau, db, zero)[1]
    has = {"zero": zero_unknown, "gain": False, "limit": False}
    head = [tau, db, coef[0]] + ([zero] if zero_unknown else [])
    best_v, best_s = fit.refine(np.array(head + list(coef[1:])), has)
    var = best_s / fit.free.sum()
    best_has = dict(has)

    def split(v, h):
        n_head = 3 + h["zero"]
        n_mid = h["gain"] + h["limit"]
        return list(v[:n_head]), list(v[n_head:n_head + n_mid]), list(v[n_head + n_mid:])

    if look_for_limit:
        head, _, tail = split(best_v, best_has)
        a, _ = fit.columns(head[0], head[1], head[3] if best_has["zero"] else 0.0)
        top = float(np.abs(a[fit.free, 0]).max()) * head[2]
        trial = dict(best_has, limit=True)
        for share in (0.9, 0.7, 0.5):
            v, s = fit.refine(np.array(head + [share * top] + tail), trial)
            if s < best_s - ACCEPT * var:
                best_v, best_s, best_has = v, s, trial
    if look_for_gain:
        head, mid, tail = split(best_v, best_has)
        trial = dict(best_has, gain=True)
        v, s = fit.refine(np.array(head + [1.0] + mid + tail), trial)
        if s < best_s - ACCEPT * var:
            best_v, best_s, best_has = v, s, trial
    return best_v[1], best_v[2]


def solve(folder, **options):
    tr = Trend(folder)
    f, zero_unknown = tr.frequency()
    trimmed = tr.setup.get("load_set_mw_is") == "plant controller set point" and tr.setup.get("plant_load_controller") == "on"
    fit = Fit(f, tr.col["injected_mhz"] / 1000.0, tr.power(), tr.col["load_set_mw"], tr.free(),
              tr.sheet["f_nominal_hz"], trimmed)
    db, k = settle(fit, zero_unknown, **options)
    droop = 100.0 * tr.sheet["rated_mw"] / (k * tr.sheet["f_nominal_hz"])
    return {"droop_percent": round(float(droop), 3), "deadband_millihertz": round(1000.0 * float(db), 2)}


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__.split("\n")[2])
    Path(sys.argv[2]).write_text(json.dumps(solve(sys.argv[1])) + "\n")
