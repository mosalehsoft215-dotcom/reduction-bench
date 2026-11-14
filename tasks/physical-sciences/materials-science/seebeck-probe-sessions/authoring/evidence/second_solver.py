"""A second reduction written apart from solution/seebeck_session.py.

    python authoring/evidence/second_solver.py <session folder> <answer.json>

Differences from the reference: cubic-spline inversion of the lot table; the bar on the rig from the log through a
searchsorted walk; V fitted on dT over t >= 330 s; the steady rise from scipy curve_fit; beta = 0.0755 (my own pairwise
estimate from D-03 and D-08); alpha per item set point, carried to the sample with a monotone PCHIP in T; the sample's
kappa from its heat flow by root finding; S(T) by a cubic in (T - 300) / 100.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import brentq, curve_fit

BETA = 0.0755


def table(path):
    with open(path) as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    return {h: np.array([float(r[i]) for r in rows[1:]]) for i, h in enumerate(head) if h != "mount"}


def main(folder, out):
    f = Path(folder)
    rig = json.loads((f / "rig.json").read_text())
    sc = table(f / "scan.csv")
    lot = table(f / rig["wire_lot"]["table"])
    to_T = CubicSpline(lot["E_AB_uV"], lot["T_K"])
    to_E = CubicSpline(lot["T_K"], lot["E_AB_uV"])
    with open(f / "log.csv") as fh:
        ev = sorted((float(r["t_s"]), r["event"], r["item"]) for r in csv.DictReader(fh) if r["event"] in ("MOUNT", "UNMOUNT"))
    stamps = np.array([e[0] for e in ev])
    def on_rig(t):
        k = np.searchsorted(stamps, t, side="right") - 1
        return ev[k][2] if k >= 0 and ev[k][1] == "MOUNT" else None

    t = sc["t_s"]
    ecj = to_E(sc["cj_K"])
    b1, b2 = to_T(sc["E1_uV"] + ecj), to_T(sc["E2_uV"] + ecj)
    leg = rig["voltage_channel"]["legs"]
    sign = 1.0 if rig["voltage_channel"]["polarity"] == "probe2_minus_probe1" else -1.0
    hot = sc["heater_mW"] > 0
    edges = np.flatnonzero(np.diff(np.r_[0, hot.astype(int), 0]))
    ramps = []
    for a, b in zip(edges[::2], edges[1::2]):
        tt = t[a:b] - t[a]
        dt = b2[a:b] - b1[a:b]
        v = sign * sc["V_uV"][a:b]
        use = tt >= 330.0
        slope = np.polyfit(dt[use], v[use], 1)[0]
        (y0, rise, tau), _ = curve_fit(lambda x, y0, r, ta: y0 + r * (1 - np.exp(-x / ta)), tt, dt, p0=[0.0, dt[-1], 150.0], maxfev=20000)
        Tm = float(np.mean((b1[a:b] + b2[a:b])[use]) / 2)
        sl = float(np.interp(Tm, lot["T_K"], lot[f"S_{leg}_uV_per_K"]))
        ramps.append(dict(bar=on_rig(t[a]), T=Tm, raw=sl - slope, rise=rise, P=sc["heater_mW"][a] * 1e-3))
    bars = rig["bars"]
    d = rig["probes"]["spacing_mm"] * 1e-3
    sec = lambda n: bars[n]["section_mm"][0] * bars[n]["section_mm"][1] * 1e-6
    alphas, etas = [], []
    for r in ramps:
        b = bars.get(r["bar"], {})
        fn = b.get("certificate") or b.get("table")
        if not fn:
            continue
        cert = table(f / fn)
        c = np.interp(r["T"], cert["T_K"], cert["S_uV_per_K"]) / r["raw"]
        k = np.interp(r["T"], cert["T_K"], cert["kappa_W_per_mK"])
        alphas.append((r["T"], 1 / c - 1 - BETA / k))
        etas.append(k * sec(r["bar"]) * r["rise"] / (c * r["P"] * d))
    eta = float(np.median(etas))
    grid = {}
    for T, a in alphas:
        grid.setdefault(round(T / 10) * 10, []).append((T, a))
    xs = np.array([np.mean([p[0] for p in v]) for _, v in sorted(grid.items())])
    ys = np.array([np.mean([p[1] for p in v]) for _, v in sorted(grid.items())])
    alpha_of = (lambda T: float(ys[0])) if len(xs) == 1 else PchipInterpolator(xs, ys, extrapolate=True)
    sample = next(n for n, b in bars.items() if b["role"] == "sample")
    pts = []
    for r in ramps:
        if r["bar"] != sample:
            continue
        al = float(alpha_of(r["T"]))
        g = lambda k: k - eta * r["P"] * d / (sec(sample) * r["rise"]) / (1 + al + BETA / k)
        kap = brentq(g, 0.05, 50.0)
        pts.append((r["T"], r["raw"] / (1 + al + BETA / kap)))
    pts.sort()
    x = (np.array([p[0] for p in pts]) - 300) / 100
    coef = np.polyfit(x, [p[1] for p in pts], 3)
    Path(out).write_text(json.dumps({"S_350K_uV_per_K": float(np.polyval(coef, 0.5)), "S_550K_uV_per_K": float(np.polyval(coef, 2.5))}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
