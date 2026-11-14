"""Second, separately written reduction of a dual-label vial (same call as the reference).

Differences from solution/lsc_reduce.py:
- Standards on the day from the decay constant ln 2 / half-life with the elapsed time in years of 365.25 d.
- Quench curves: mis-dispensed vials are found by least trimmed squares: of the fifteen four-vial subsets of a
  series, the one whose quadratic (on the main window, A for 3H and B for 14C) has the smallest sum of squared
  relative residuals is taken, and each of the other two vials is kept only when it lies within 3 % of that quadratic;
  the kept vials' quadratics give the curves of both windows.
- The sample's model is made linear: for a fixed colour fraction f, glass tSIE q and glow time constant tau, the
  window rates and the SIS relation rewritten as
      (SIS_k - s_H) n_H + (SIS_k - s_C) n_C + (SIS_k - 4) g_k = n_bkg (SIS_bkg - SIS_k)
  are linear in H, C and L0, and are solved by bounded linear least squares (L0 >= 0), each row over its scatter
  (Poisson counts; 0.3 % of SIS times the count's total rate).  f (0 to 1 by 0.02, glass at the mean tSIE), q (in
  1.0 steps from the mean tSIE less 90 to less 40, polyethylene, f 0 or 1) and tau (40 values from 5 to 600 min) are
  searched on grids, the polyethylene choices carrying 2 more in the sum of squares; the best grid point is refined
  around it (tau on 11 values across the neighbouring grid steps, f in steps of 0.004 within 0.02, q in steps of
  0.1 within 1).
- A sample on counter B: the window efficiencies scaled by the monthly normalization's net-rate ratios B/A of the same
  unquenched sources, the background taken from B's blank there.
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

YEAR = 365.25
HALF = {"H3": 12.32, "C14": 5700.0}


def reduce_vial(folder):
    folder = Path(folder)
    meta = json.loads((folder / "vial.json").read_text())
    tab = np.genfromtxt(folder / "counts.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    lab = np.array([str(x) for x in tab["label"]])
    minutes = tab["count_min"].astype(float)
    ra, rb = tab["counts_a"] / minutes, tab["counts_b"] / minutes
    tsie, sis = tab["tsie"].astype(float), tab["sis"].astype(float)
    b = lab == "background"
    bga, bgb, bgs = float(ra[b][0]), float(rb[b][0]), float(sis[b][0])
    nb = bga + bgb
    day = np.datetime64(meta["counted"])
    eff, sline = {}, {}
    for kind in ("chemical", "colour"):
        for nuc in ("H3", "C14"):
            sel = np.char.startswith(lab, f"std-{nuc}-{kind}-")
            std = meta["standards"][nuc]
            years = (day - np.datetime64(std["reference_date"])).astype(float) / YEAR
            act = std["activity_bq"] * np.exp(-np.log(2) / HALF[nuc] * years)
            q = tsie[sel]
            ea, eb = (ra[sel] - bga) / (60 * act), (rb[sel] - bgb) / (60 * act)
            main = ea if nuc == "H3" else eb
            best = None
            for sub in itertools.combinations(range(q.size), q.size - 2):
                sub = list(sub)
                coef = np.polyfit(q[sub], main[sub], 2)
                sse = float(np.sum((main[sub] / np.polyval(coef, q[sub]) - 1) ** 2))
                if best is None or sse < best[0]:
                    best = (sse, sub, coef)
            keep = sorted(best[1] + [i for i in range(q.size) if i not in best[1]
                                     and abs(main[i] / np.polyval(best[2], q[i]) - 1) <= 0.03])
            eff[(kind, nuc, "a")] = np.polyfit(q[keep], ea[keep], 2)
            eff[(kind, nuc, "b")] = np.polyfit(q[keep], eb[keep], 2)
            net = (ra[sel] - bga) + (rb[sel] - bgb)
            alone = (sis[sel] * (net + nb) - nb * bgs) / net
            sline[(kind, nuc)] = np.polyfit(q, alone, 1)
    cyc = np.char.startswith(lab, "sample-cycle-")
    t0, dur = tab["start_min"][cyc].astype(float), minutes[cyc]
    scale = {"a_h": 1.0, "a_c": 1.0, "b_c": 1.0}
    ctr = str(tab["counter"][cyc][0])
    if ctr != "A":
        nm = json.loads((folder / "counters.json").read_text())
        per = nm["minutes"]
        r_ = lambda blk, c, w: blk[c][w] / per  # noqa: E731
        scale = {"a_h": (r_(nm["tritium_source"], ctr, "A") - r_(nm["blank"], ctr, "A"))
                 / (r_(nm["tritium_source"], "A", "A") - r_(nm["blank"], "A", "A")),
                 "a_c": (r_(nm["carbon14_source"], ctr, "A") - r_(nm["blank"], ctr, "A"))
                 / (r_(nm["carbon14_source"], "A", "A") - r_(nm["blank"], "A", "A")),
                 "b_c": (r_(nm["carbon14_source"], ctr, "B") - r_(nm["blank"], ctr, "B"))
                 / (r_(nm["carbon14_source"], "A", "B") - r_(nm["blank"], "A", "B"))}
        bga, bgb = r_(nm["blank"], ctr, "A"), r_(nm["blank"], ctr, "B")
        nb = bga + bgb
    ya, yb, s = ra[cyc] - bga, rb[cyc] - bgb, sis[cyc]
    sa, sb = np.sqrt(tab["counts_a"][cyc]) / dur, np.sqrt(tab["counts_b"][cyc]) / dur
    total = ra[cyc] + rb[cyc]
    q_glass = float(np.mean(tsie[cyc]))

    def mix(f, key, q):
        lo, hi = ("chemical",) + key, ("colour",) + key
        table = eff if len(key) == 2 else sline
        return (1 - f) * np.polyval(table[lo], q) + f * np.polyval(table[hi], q)

    def solve(f, q, tau):
        e_ha = mix(f, ("H3", "a"), q) * scale["a_h"]
        e_ca, e_cb = mix(f, ("C14", "a"), q) * scale["a_c"], mix(f, ("C14", "b"), q) * scale["b_c"]
        s_h, s_c = mix(f, ("H3",), q), mix(f, ("C14",), q)
        phi = tau / dur * (np.exp(-t0 / tau) - np.exp(-(t0 + dur) / tau))
        z = np.zeros_like(phi)
        rows_a = np.column_stack([60 * e_ha + z, 60 * e_ca + z, phi]) / sa[:, None]
        rows_b = np.column_stack([z, 60 * e_cb + z, z]) / sb[:, None]
        sig = 0.002 * s * total
        rows_s = np.column_stack([(s - s_h) * 60 * e_ha, (s - s_c) * 60 * (e_ca + e_cb), (s - 4.0) * phi]) / sig[:, None]
        mat = np.vstack([rows_a, rows_b, rows_s])
        rhs = np.r_[ya / sa, yb / sb, nb * (bgs - s) / sig]
        res = lsq_linear(mat, rhs, bounds=([0.0, 0.0, 0.0], [np.inf, np.inf, np.inf]))
        return float(np.sum((mat @ res.x - rhs) ** 2)), res.x

    taus = np.geomspace(5.0, 600.0, 40)
    trials = []
    for f in np.arange(0.0, 1.0001, 0.02):
        for tau in taus:
            trials.append((solve(f, q_glass, tau)[0], f, q_glass, tau, False))
    for f in (0.0, 1.0):
        for q in np.arange(q_glass - 90.0, q_glass - 39.99, 1.0):
            for tau in taus:
                trials.append((solve(f, q, tau)[0] + 2.0, f, q, tau, True))
    _, f, q, tau, plastic = min(trials, key=lambda r: r[0])
    step_tau = taus[1] / taus[0]
    fine = []
    for tf in np.geomspace(tau / step_tau, tau * step_tau, 11):
        if plastic:
            for qq in np.arange(q - 1.0, q + 1.001, 0.1):
                if q_glass - 90.0 <= qq <= q_glass - 40.0:
                    fine.append((solve(f, qq, tf)[0], f, qq, tf))
        else:
            for ff in np.arange(max(f - 0.02, 0.0), min(f + 0.02, 1.0) + 1e-9, 0.004):
                fine.append((solve(ff, q, tf)[0], ff, q, tf))
    _, f, q, tau = min(fine, key=lambda r: r[0])
    h, c, _ = solve(f, q, tau)[1]
    return {"tritium_bq": float(h), "carbon14_bq": float(c)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_vial(sys.argv[1])))
