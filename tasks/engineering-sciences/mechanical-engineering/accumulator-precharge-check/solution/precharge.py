"""Nitrogen precharge at 20 C and gas volume of a piston accumulator from one discharge check.

    python3 precharge.py <bottle folder> <output file>

1. Meter factor from the prover check: pulses per litre = mean prover count / 2.000 L; drained oil volume from the
   pulse count.
2. Rests: stretches with a constant pulse count and the oil side above 2 bar.  In a rest the gas volume is fixed and
   its temperature relaxes to the wall, so each rest is fitted with p = p_inf + A exp(-t/tau), with tau
   taken from the longest rest (scanned), then p_inf and A per rest by linear least squares.
3. The rest end-points p_inf are isothermal states at the ambient temperature: p = n R T / (V - n B(T)) with the van
   der Waals second virial coefficient B(T) = b - a/(R T) and V = V_start + drained volume.  n and V_start come from a
   least-squares fit.
4. Gas volume on the stop = V_start + total drained volume; precharge at 20 C = n R T20 / (V - n B(T20)).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

R, A_VDW, B_VDW, T20 = 8.314462618, 0.1370, 3.87e-5, 293.15


def virial_b(T):
    return B_VDW - A_VDW / (R * T)


def rests(t, p, pulses):
    out, start = [], 0
    for i in range(1, len(t) + 1):
        if i == len(t) or pulses[i] != pulses[start]:
            seg = slice(start, i)
            if i - start > 60 and p[start] > 2.0:
                out.append((t[seg][10:] - t[seg][10], p[seg][10:], pulses[start]))
            start = i
    return out


def wall_time_constant(segment):
    tt, pp, _ = segment
    best = None
    for tau in np.geomspace(40, 2000, 160):
        X = np.c_[np.ones_like(tt), np.exp(-tt / tau)]
        coef, *_ = np.linalg.lstsq(X, pp, rcond=None)
        res = float(np.sum((X @ coef - pp) ** 2))
        if best is None or res < best[0]:
            best = (res, tau)
    return best[1]


def asymptotes(segments):
    tau = wall_time_constant(max(segments, key=lambda s: len(s[0])))
    ends = []
    for tt, pp, _ in segments:
        X = np.c_[np.ones_like(tt), np.exp(-tt / tau)]
        ends.append(np.linalg.lstsq(X, pp, rcond=None)[0][0])
    return ends


def reduce(folder):
    sheet = json.loads((folder / "bottle.json").read_text())
    with open(folder / "discharge.csv") as f:
        rows = [(float(r["t_s"]), float(r["oil_bar"]), int(r["meter_pulses"])) for r in csv.DictReader(f)]
    t, p, pulses = (np.array(c) for c in zip(*rows))
    live = np.nonzero(p > 5.0)[0][-1] + 1                     # the piston reaches its stop after this row
    t, p, pulses_all, pulses = t[:live], p[:live], pulses, pulses[:live]
    k = np.mean(sheet["prover_check"]["pulses"]) / sheet["prover_check"]["vessel_L"]
    segs = rests(t, p, pulses)
    p_inf = np.array(asymptotes(segs)) * 1e5
    drained = np.array([s[2] for s in segs]) / k / 1000.0
    Tw = sheet["ambient_C"] + 273.15

    def resid(x):
        n, v_start = x
        return (n * R * Tw / (v_start + drained - n * virial_b(Tw)) - p_inf) / 1e5

    guess_n = p_inf[0] * 0.01 / (R * Tw)
    fit = least_squares(resid, x0=[guess_n, 0.01], x_scale=[guess_n, 0.01])
    n, v_start = fit.x
    v_stop = v_start + pulses_all[-1] / k / 1000.0
    p20 = n * R * T20 / (v_stop - n * virial_b(T20))
    return p20 / 1e5, v_stop * 1000.0


def main():
    p20, vol = reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps({"precharge_bar_20c": round(p20, 4), "gas_volume_l": round(vol, 4)}) + "\n")


if __name__ == "__main__":
    main()
