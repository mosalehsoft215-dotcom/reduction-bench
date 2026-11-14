"""Nitrogen precharge at 20 C and gas volume of a piston accumulator from one stepped discharge.

    python3 charge_audit.py <check folder> <output file>

Model: shell temperature Ts(t) = Tamb + dT exp(-t/taus); in every rest
p(t) = n R Ts(t) / (V - n B(Ts)) + A exp(-(t - t_rest)/tau); globals (n, Vstart, dT, taus, tau) by least squares with
the per-rest amplitudes A solved linearly."""
import csv, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

R, A_VDW, B_VDW, T20 = 8.314462618, 0.1370, 3.87e-5, 293.15


def b(T):
    return B_VDW - A_VDW / (R * T)


def reduce(folder):
    sheet = json.loads((folder / "accu_sheet.json").read_text())
    rows = [(float(r["t_sec"]), float(r["fluid_bara"]), int(r["meter_count"])) for r in csv.DictReader(open(folder / "drain_log.csv"))]
    t, p, pu = (np.array(c) for c in zip(*rows))
    live = np.nonzero(p > 5.0)[0][-1] + 1
    total = pu[-1]
    t, p, pu = t[:live], p[:live], pu[:live]
    k = np.mean(sheet["prover_check"]["pulses"]) / sheet["prover_check"]["vessel_L"]
    segs, start = [], 0
    for i in range(1, len(t) + 1):
        if i == len(t) or pu[i] != pu[start]:
            if i - start > 60:
                segs.append((t[start + 10:i], p[start + 10:i] * 1e5, pu[start] / k / 1000.0))
            start = i
    Tamb = sheet["ambient_C"] + 273.15

    def resid(x):
        n, vs, dT, taus, tau = x
        out = []
        for tt, pp, d in segs:
            Ts = Tamb + dT * np.exp(-tt / taus)
            base = n * R * Ts / (vs + d - n * b(Ts))
            e = np.exp(-(tt - tt[0]) / tau)
            amp = np.dot(e, pp - base) / np.dot(e, e)
            out.append((base + amp * e - pp) / 1e5)
        return np.concatenate(out)

    p0 = segs[0][1][-1]
    g = p0 * 0.01 / (R * Tamb)
    fit = least_squares(resid, x0=[g, 0.01, 0.0, 3000.0, 300.0], x_scale=[g, 0.01, 5.0, 1000.0, 100.0],
                        bounds=([0, 1e-4, -30, 300, 20], [np.inf, 1.0, 30, 30000, 3000]))
    n, vs, dT, taus, tau = fit.x
    v = vs + total / k / 1000.0
    return n * R * T20 / (v - n * b(T20)) / 1e5, v * 1000.0


def main():
    p20, vol = reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps({"nitrogen_charge_bar": round(p20, 4), "gas_space_l": round(vol, 4)}) + "\n")


if __name__ == "__main__":
    main()
