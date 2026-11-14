"""Second, separately written reduction of a DVS run: the two pools as first-order filters, Cauchy loss.

Each water pool is a first-order recursive filter of its target sequence over the one-minute readings
(x_i = a x_(i-1) + (1 - a) y_i, a = exp(-1 min / tau), y its share f or 1 - f of the GAB content of the current RH),
starting from free initial contents, with an optional solvent mass S exp(-t / tau_s) on top.  All parameters (m0, Wm, C, K, f, tau1, tau2, the
initial contents, S, tau_s) are fitted to every logged mass with a Cauchy loss of scale 0.3 ug, so that apparent-mass
excursions (a charged pan) barely count; several starts, the solvent term kept when it at least halves the
robust misfit.  Stretches whose 21-minute running median residual exceeds 0.15 ug are then dropped (widened by 21
minutes) and the fit repeated by plain least squares on the rest, three times at most.  100 Wm and 100 W(0.75) come from the fitted GAB parameters.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import lfilter, medfilt


def gab(a, wm, c, k):
    return wm * c * k * a / ((1 - k * a) * (1 - k * a + c * k * a))


def reduce_run(folder):
    folder = Path(folder)
    log = np.genfromtxt(folder / "mass_log.csv", delimiter=",", names=True)
    t, a_w, mass = log["minute"], log["rh_target_pct"] / 100.0, log["mass_mg"]
    dt = np.diff(np.r_[t[0] - 1.0, t])

    def simulate(q, solvent):
        # each pool is a first-order filter of its target sequence (one-minute readings): x_i = a x_(i-1) + (1 - a) y_i
        m0, wm, c, k, f, l1, l2, p1, p2 = q[:9]
        target = gab(a_w, wm, c, k) * m0
        out = np.full(t.size, m0)
        for share, tau, x0 in ((f, math.exp(l1), p1), (1 - f, math.exp(l1) + math.exp(l2), p2)):
            a = math.exp(-1.0 / tau)
            out = out + lfilter([1 - a], [1, -a], share * target, zi=[a * x0])[0]
        if solvent:
            out = out + q[9] * np.exp(-t / math.exp(q[10]))
        return out

    ref = float(mass[np.argmax(a_w > 0) - 1])

    keep = np.ones(t.size, bool)

    def fit(solvent, starts=None, loss="cauchy"):
        lo = [0.5 * ref, 1e-3, 0.5, 0.3, 0.05, math.log(0.3), math.log(0.5), 0.0, 0.0]
        hi = [1.5 * ref, 0.3, 300.0, 0.99, 0.99, math.log(60.0), math.log(3000.0), 0.3 * ref, 0.3 * ref]
        if solvent:
            lo, hi = lo + [0.0, math.log(150.0)], hi + [0.03 * ref, math.log(20000.0)]
        best = None
        trial = starts or [[ref, 0.04, 10.0, 0.8, f0, math.log(3.0), l2, 0.02 * ref, 0.01 * ref]
                           + ([0.005 * ref, math.log(600.0)] if solvent else [])
                           for l2 in (math.log(8.0), math.log(150.0), math.log(400.0)) for f0 in (0.6, 0.85)]
        for q0 in trial:
            if True:
                sol = least_squares(lambda q: (simulate(q, solvent) - mass)[keep], np.clip(q0, lo, hi), bounds=(lo, hi),
                                    x_scale="jac", loss=loss, f_scale=3e-4, max_nfev=400)
                if best is None or sol.cost < best.cost:
                    best = sol
        return best

    plain = fit(False)
    with_s = fit(True)
    solvent = with_s.cost < 0.5 * plain.cost
    best = with_s if solvent else plain
    # stretches whose 21-minute running median residual exceeds 0.15 ug (over five times its noise) are dropped,
    # widened by 21 minutes, and the fit repeated from where it stood, three times at most
    for _ in range(3):
        r = simulate(best.x, solvent) - mass
        bad = np.abs(medfilt(r, 21)) > 1.5e-4
        bad = np.convolve(bad.astype(float), np.ones(21), "same") > 0
        if not bad.any() or bad.mean() > 0.5 or np.array_equal(~bad, keep):
            break
        keep[:] = ~bad
        best = fit(solvent, [best.x], loss="linear")
    _, wm, c, k = best.x[:4]
    return {"monolayer_water_pct": float(100 * wm), "water_at_75rh_pct": float(100 * gab(0.75, wm, c, k))}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_run(sys.argv[1])) + "\n")
