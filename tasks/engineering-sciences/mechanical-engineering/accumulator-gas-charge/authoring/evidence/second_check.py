"""Independent reduction of one accumulator unit (does not import the reference).

Differences from the reference: the gas temperature is carried through the whole record instead of a free amplitude
per pause.  Inside a pause (fixed volume) the gas follows the shell exactly,
T(t) = Tamb + D' e^(-t/ts) + (T(t0) - Tamb - D' e^(-t0/ts)) e^(-(t-t0)/tau) with D' = D ts / (ts - tau);
through each release the gas energy n cv dT = -p dV + (n cv / tau)(Ts - T) dt is stepped at 0.1 s with the volume
taken from the pulse totals second by second.  Unknowns
(n, gas volume at the first pause, D, ts, tau, gas temperature at the first sample) by least squares on every pause
sample; the meter factor is the booked one.

    python second_check.py <unit folder> <output file>
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

R, CV, A, B, T20 = 8.314462618, 20.8, 0.1370, 3.87e-5, 293.15


def bv(T):
    return B - A / (R * T)


def blocks(folder):
    sheet = json.loads((folder / "unit_sheet.json").read_text())
    with open(folder / "oil_side_log.csv") as f:
        rows = [(float(r["time_s"]), float(r["oil_bar_abs"]), int(r["gear_pulses"])) for r in csv.DictReader(f)]
    t, p, n = (np.array(c) for c in zip(*rows))
    k = sheet["meter"]["k_factor_pulses_per_L"]
    end = int(np.nonzero(p > 5.0)[0][-1]) + 1
    moving = np.r_[n[0] != 0, np.diff(n) != 0][:end]
    edges = np.r_[0, np.nonzero(np.diff(moving.astype(int)))[0] + 1, end]
    out = [(bool(moving[a]), t[a:b], p[a:b] * 1e5, n[a:b] / k / 1000.0) for a, b in zip(edges[:-1], edges[1:])]
    return sheet, out, n[-1] / k / 1000.0


def model(x, segs, Tamb):
    """predicted pressure on every pause sample; the gas state is carried through the releases"""
    n, v0, D, ts, tau, T0 = x
    Dp = D * ts / (ts - tau)
    shell = lambda tt: Tamb + Dp * np.exp(-tt / ts)
    T, tc, vc, out = T0, 0.0, v0, []
    for moving, tt, _, drained in segs:
        if not moving:
            v = v0 + drained[0]
            Tp = shell(tt) + (T - shell(tc)) * np.exp(-(tt - tc) / tau)
            out.append(n * R * Tp / (v - n * bv(Tp)))
            T, tc, vc = Tp[-1], tt[-1], v
            continue
        for ti, di in zip(tt, drained):
            dv, dt = (v0 + di - vc) / 10.0, (ti - tc) / 10.0
            for _ in range(10):
                pg = n * R * T / (vc - n * bv(T))
                T += (-pg * dv + n * CV * (shell(tc) - T) * dt / tau) / (n * CV)
                vc, tc = vc + dv, tc + dt
    return out


def residual(x, segs, Tamb):
    quiet = [s for s in segs if not s[0]]
    return np.concatenate([(m[2:] - s[2][2:]) / 1e5 for m, s in zip(model(x, segs, Tamb), quiet)])


def main(folder, dest):
    sheet, segs, total = blocks(Path(folder))
    Tamb = sheet["ambient_C"] + 273.15
    p_last = segs[0][2][-1]
    best = None
    for D0 in (-15.0, 0.0, 15.0):
        g = p_last * 0.01 / (R * (Tamb + D0 * 0.6))
        x0 = [g, 0.01, D0, 4000.0, 300.0, Tamb + D0 + 4.0]
        fit = least_squares(residual, args=(segs, Tamb),
                            x0=x0, x_scale=[g, 0.01, 5.0, 1000.0, 100.0, 5.0],
                            bounds=([0, 1e-4, -40, 600, 20, 200], [np.inf, 1.0, 40, 40000, 2000, 400]))
        if best is None or fit.cost < best.cost:
            best = fit
    n, v0 = best.x[:2]
    v = v0 + total
    Path(dest).write_text(json.dumps({"n2_precharge_bar": n * R * T20 / (v - n * bv(T20)) / 1e5,
                                      "gas_volume_litres": v * 1000.0}))


if __name__ == "__main__":
    main(*sys.argv[1:3])
