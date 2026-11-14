"""Independently written reduction of a tape piece (does not import the reference).

Differences from solution/tape_ic.py: the fitted shunt is chosen by least squares of the clamp against every register
resistance; the thermoelectric settling is fitted with scipy curve_fit; Ic and n come from log-log interpolation of
the 0.1 and 1 uV/cm crossings instead of a straight-line fit; the bubble point is solved with brentq on one equation in
T (x eliminated); the wet test uses a 0.5 cm margin and the thermometers within 30 s of the crossing.

    python second_reduction.py <tape folder> <answer file> <rig_tables.json>
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, curve_fit

GAIN, OFFSET = 1.0042, 0.10          # clamp, read off the worked pieces (clamp / (shunt_mV / booked R))


def hhmm(s):
    a, b = s.split(":")
    return 60 * int(a) + int(b)


def main(folder, answer, tables_path):
    tab = json.loads(Path(tables_path).read_text())
    folder = Path(folder)
    sheet = json.loads((folder / "tape.json").read_text())
    with open(folder / "ramps.csv") as f:
        rows = list(csv.DictReader(f))
    col = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}

    # current: which register resistance makes the clamp agree best
    clamp_true = (col["clamp_A"] - OFFSET) / GAIN
    keep = clamp_true > 0.3 * clamp_true.max()
    cands = {sheet["shunt"]["serial"]: sheet["shunt"]["resistance_ohm"]}
    cands.update({r["serial"]: r["resistance_ohm"] for r in tab["shunt_register"]})
    score = {s: float(np.mean((col["shunt_mV"][keep] * 1e-3 / R - clamp_true[keep]) ** 2)) for s, R in cands.items()}
    booked = sheet["shunt"]["serial"]
    pick = min(score, key=score.get)
    R = cands[pick] if score[pick] < 0.25 * score[booked] else cands[booked]
    I = col["shunt_mV"] * 1e-3 / R

    # settling of the outer-pair offset, from the zero-current dwells
    log = (folder / "bench_log.txt").read_text().splitlines()
    t_low = hhmm(next(l for l in log if "lowered" in l)[:5])
    t_go = hhmm(next(l for l in log if "ramps started" in l)[:5])
    s_since = col["t_s"] + 60 * (t_go - t_low)
    z = np.abs(I) < 0.3
    model = lambda s, c, a, k: c + a * np.exp(-k * s)
    try:
        p, _ = curve_fit(model, s_since[z], col["v_outer_uV"][z], p0=(0.0, 1.0, 1 / 400), maxfev=20000)
        v = col["v_outer_uV"] - model(s_since, *p)
    except RuntimeError:
        v = col["v_outer_uV"] - np.median(col["v_outer_uV"][z])

    # bath
    with open(folder / "dewar_log.csv") as f:
        dl = list(csv.DictReader(f))
    with open(folder / "room_log.csv") as f:
        rl = list(csv.DictReader(f))
    lt = np.array([hhmm(r["time"]) for r in dl], float)
    lv = np.array([float(r["level_cm"]) for r in dl])
    y = np.median([float(r["vent_o2_pct"]) for r in dl]) / 100
    elev = tab["sites"][sheet["site"]]["elevation_m"]
    P = np.median([float(r["pressure_hPa"]) for r in rl]) * (1 - 2.25577e-5 * elev) ** 5.25588
    vp = tab["vapour_pressure_hPa"]
    lnN2 = lambda T: np.interp(T, vp["T_K"], np.log(vp["n2"]))
    lnO2 = lambda T: np.interp(T, vp["T_K"], np.log(vp["o2"]))
    T_boil = brentq(lambda T: (1 - y * P / math.exp(lnO2(T))) * math.exp(lnN2(T)) + y * P - P, 66.0, 90.0)

    lf = tab["rb4_lift"]
    out_ic, out_n = [], []
    for k in sorted(set(col["ramp"])):
        m = col["ramp"] == k
        Ik, vk, tk = I[m], v[m], col["t_s"][m]
        top = int(np.argmax(Ik))
        Ik, vk, tk = Ik[:top + 1], vk[:top + 1], tk[:top + 1]
        Imax = Ik.max()
        lo = (Ik > 0.15 * Imax / 1.13) & (Ik < 0.45 * Imax / 1.13)
        b = np.polyfit(Ik[lo], vk[lo], 1)
        E = (vk - np.polyval(b, Ik)) / sheet["taps"]["outer_cm"]
        # current at a given field: first upward crossing, interpolated in log-log
        def cross(level):
            idx = np.where((E[:-1] < level) & (E[1:] >= level) & (Ik[:-1] > 0.6 * Imax))[0]
            j = idx[-1]
            f = (math.log(level) - math.log(max(E[j], 1e-9))) / (math.log(E[j + 1]) - math.log(max(E[j], 1e-9)))
            return math.exp(math.log(Ik[j]) + f * (math.log(Ik[j + 1]) - math.log(Ik[j]))), tk[j]
        # average plateau readings first so noise does not move the crossings
        order = np.argsort(Ik)
        Is, Es = Ik[order], E[order]
        win = max(5, len(Is) // 200)
        Ism = np.convolve(Is, np.ones(win) / win, mode="same")
        Esm = np.convolve(Es, np.ones(win) / win, mode="same")
        good = (Esm > 0.05) & (Esm < 2.0)
        I1 = math.exp(np.interp(0.0, np.log(Esm[good]), np.log(Ism[good])))
        I01 = math.exp(np.interp(math.log(0.1), np.log(Esm[good]), np.log(Ism[good])))
        n = math.log(10.0) / math.log(I1 / I01)
        tc = float(np.interp(I1, Is, col["t_s"][m][:top + 1][order]))
        lev = float(np.interp(t_go + tc / 60, lt, lv))
        if lev > sheet["holder"]["thermometer_height_cm"] + 0.5:
            near = m & (np.abs(col["t_s"] - tc) <= 30)
            T = float(np.mean((col["pt100_K"][near] + col["diode_K"][near]) / 2))
        else:
            T = T_boil
        out_ic.append(I1 / np.interp(T, lf["T_K"], lf["ic_ratio"]))
        out_n.append(n / np.interp(T, lf["T_K"], lf["n_ratio"]))
    Path(answer).write_text(json.dumps({"ic77_a": float(np.mean(out_ic)), "n77": float(np.mean(out_n))}))


if __name__ == "__main__":
    main(*sys.argv[1:4])
