"""A second blow planner, written apart from solution/blow_plan.py.

Where it differs:
  * trial values are medians of each window, not means;
  * SH1 outlet temperature: median of the two test thermocouples and the station one (no screening rule);
    feedwater flow likewise: median of the two nozzle transmitters and the clamp-on meter where one was fitted;
  * drum inventory: water and steam mass of the drum segment, from the true level at the first and the last
    minute of each trial (the DP transmitter undone sample by sample), not a level slope;
  * feedwater density correction with the window's median feed state;
  * counter-flow banks by composite Simpson quadrature of dQ / (T_gas - T_steam), scalar IF97 backward
    equations, and an outer fixed-point loop on the gas temperature between SH2 and SH1;
  * spray that passes the shut valve: from the attemperator balance of each trial, then the median;
  * design by a bracketed root search on the drum pressure with K on 31 points per section, aiming at K = 1.425.
python second_planner.py <records dir> <answer.json>
"""
import csv
import json
import math
import re
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from iapws import IAPWS97
from iapws.iapws97 import _Backward2_T_Ph, _Region1, _Region2, _TSat_P
from scipy.optimize import brentq

AIM = 1.425
N = 31


@lru_cache(maxsize=None)
def _water_steam(P):
    Ts = _TSat_P(P / 10.0)
    return Ts - 273.15, _Region2(Ts, P / 10.0), _Region1(Ts, P / 10.0)


def water_steam(P):
    return _water_steam(round(P, 7))


def T_of(P, h):
    Ts, g, _ = water_steam(P)
    return Ts if h <= g["h"] else _Backward2_T_Ph(P / 10.0, h) - 273.15


def v_of(P, h):
    Ts, g, f = water_steam(P)
    if h > g["h"]:
        return _Region2(_Backward2_T_Ph(P / 10.0, h), P / 10.0)["v"]
    return f["v"] + (h - f["h"]) / (g["h"] - f["h"]) * (g["v"] - f["v"])


def h_of(P, T):
    return IAPWS97(P=P / 10.0, T=T + 273.15).h


def dens(P, T):
    return 1.0 / IAPWS97(P=P / 10.0, T=T + 273.15).v


def segment_mass(drum, P, level):
    r, L = drum["inner_diameter_m"] / 2.0, drum["length_m"]
    y = min(max(level, -r), r)
    depth = r + y
    area = r * r * math.acos((r - depth) / r) - (r - depth) * math.sqrt(max(2 * r * depth - depth * depth, 0.0))
    _, g, f = water_steam(P)
    return area * L / f["v"] + (math.pi * r * r - area) * L / g["v"]


def read(folder):
    sheet = json.loads((folder / "circuit.json").read_text())
    with open(folder / "trial_blows.csv") as fh:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    return sheet, rows


def windows(sheet, rows):
    ins = sheet["instruments"]
    fixed = re.search(r"that of feedwater at ([\d.]+) bar abs and ([\d.]+) C", ins["F_fw_A_kg_s"])
    lv = re.search(r"taps (\d+) mm above and (\d+) mm below.*?about ([\d.]+) C.*?at ([\d.]+) bar abs", ins["L_drum_mm"])
    up, down, legT, rng = (float(x) for x in lv.groups())
    up, down = up / 1000.0, down / 1000.0
    head = re.search(r"([\d.]+) m below its tapping.*?about ([\d.]+) C", ins["P_sh1_out_bar"])
    out = []
    for n in sorted({r["trial"] for r in rows}):
        part = [r for r in rows if r["trial"] == n]
        med = {k: float(np.median([p[k] for p in part])) for k in part[0]}
        P = med["P_drum_bar"]
        _, gc, fc = water_steam(rng)
        _, g, f = water_steam(P)
        leg = dens(P, legT)

        def level(shown_mm):
            dp = (leg - 1 / gc["v"]) * (up + down) - (1 / fc["v"] - 1 / gc["v"]) * (shown_mm / 1000.0 + down)
            return ((leg - 1 / g["v"]) * (up + down) - dp) / (1 / f["v"] - 1 / g["v"]) - down
        first = [p for p in part if p["t_s"] <= part[0]["t_s"] + 60]
        last = [p for p in part if p["t_s"] >= part[-1]["t_s"] - 60]
        m0 = segment_mass(sheet["drum"], P, np.median([level(p["L_drum_mm"]) for p in first]))
        m1 = segment_mass(sheet["drum"], P, np.median([level(p["L_drum_mm"]) for p in last]))
        dt = np.median([p["t_s"] for p in last]) - np.median([p["t_s"] for p in first])
        scale = 1.0
        if fixed:
            scale = math.sqrt(dens(med["P_fw_bar"], med["T_fw_C"]) / dens(float(fixed.group(1)), float(fixed.group(2))))
        meters = [scale * med["F_fw_A_kg_s"], scale * med["F_fw_B_kg_s"]]
        if "F_fw_clamp_kg_s" in med:
            meters.append(med["F_fw_clamp_kg_s"])
        feed = float(np.median(meters))
        P1 = med["P_sh1_out_bar"]
        if head:
            P1 -= dens(P1, float(head.group(2))) * 9.80665 * float(head.group(1)) / 1e5
        T1 = float(np.median([med["T_sh1_out_A_C"], med["T_sh1_out_B_C"], med["T_sh1_out_C"]]))
        h1, hm = h_of(P1, T1), h_of(P1, med["T_sh2_in_C"])
        hw = h_of(med["P_fw_bar"], med["T_fw_C"])
        r = max((h1 - hm) / (hm - hw), 0.0)
        m = (feed - med["F_cbd_kg_s"] + (m0 - m1) / dt) / (1 + r)
        out.append(dict(P=P, P1=P1, h0=g["h"], h1=h1, hm=hm, m=m, leak=m * r, T_fw=med["T_fw_C"],
                        dp_fw=med["P_fw_bar"] - P, P2=med["P_sh2_out_bar"], h2=h_of(med["P_sh2_out_bar"],
                                                                                    med["T_sh2_out_C"]),
                        P3=med["P_ms_end_bar"], h3=h_of(med["P_ms_end_bar"], med["T_ms_end_C"]),
                        P4=med["P_tp_end_bar"], h4=h_of(med["P_tp_end_bar"], med["T_tp_end_C"])))
    return out


def needed_UA(Pa, Pb, ha, hb, m, gas_hot, cap):
    """Composite Simpson (32 intervals per part) of m dh / (T_gas - T_steam), split at saturation."""
    gas_cold = gas_hot - m * (hb - ha) / cap
    hg = water_steam(Pa)[1]["h"]
    parts = [(ha, hg, True), (hg, hb, False)] if ha < hg < hb else [(ha, hb, hb <= hg)]
    total = 0.0
    for lo, hi, wet in parts:
        n = 32
        step = (hi - lo) / n
        acc = 0.0
        for i in range(n + 1):
            h = lo + i * step
            P = Pa + (Pb - Pa) * (h - ha) / (hb - ha)
            T = _TSat_P(P / 10.0) - 273.15 if wet else _Backward2_T_Ph(P / 10.0, max(h, hg)) - 273.15
            gap = gas_cold + m * (h - ha) / cap - T
            if gap <= 0:
                return 1e12
            acc += (1 if i in (0, n) else 4 if i % 2 else 2) * m / gap
        total += acc * step / 3.0
    return total


def kcoef(Pa, Pb, ha, hb, m):
    return (Pa - Pb) * 1e5 / (m * m * v_of(0.5 * (Pa + Pb), 0.5 * (ha + hb)))


def fit(sheet, rows):
    gs = sheet["gas_side"]
    cap, Tg = gs["exhaust_flow_kg_s"] * gs["exhaust_cp_kJ_kgK"], gs["gas_temperature_entering_SH2_C"]
    per = []
    for w in windows(sheet, rows):
        m, m2 = w["m"], w["m"] + w["leak"]
        mid = Tg - m2 * (w["h2"] - w["hm"]) / cap
        per.append(dict(UA1=needed_UA(w["P"], w["P1"], w["h0"], w["h1"], m, mid, cap),
                        UA2=needed_UA(w["P1"], w["P2"], w["hm"], w["h2"], m2, Tg, cap),
                        k1=kcoef(w["P"], w["P1"], w["h0"], w["h1"], m), k2=kcoef(w["P1"], w["P2"], w["hm"], w["h2"], m2),
                        k3=kcoef(w["P2"], w["P3"], w["h2"], w["h3"], m2), k4=kcoef(w["P3"], w["P4"], w["h3"], w["h4"], m2),
                        C=m2 / math.sqrt(w["P4"] * 1e5 / v_of(w["P4"], w["h4"])), leak=w["leak"], T_fw=w["T_fw"],
                        dp_fw=w["dp_fw"]))
    model = {k: float(np.median([p[k] for p in per])) for k in per[0]}
    model.update(cap=cap, Tg=Tg)
    return model


def pdrop(Pa, ha, hb, m, k):
    Pb = Pa
    for _ in range(300):
        nxt = Pa - k * m * m * v_of(0.5 * (Pa + Pb), 0.5 * (ha + hb)) / 1e5
        if nxt < 0.5:
            return None
        if abs(nxt - Pb) < 1e-10:
            return nxt
        Pb = nxt
    return None


def heated(Pa, ha, m, UA, k, gas_hot, cap):
    Pb = Pa
    for _ in range(40):
        top = ha + 10.0
        while needed_UA(Pa, Pb, ha, top, m, gas_hot, cap) < UA:
            top += 0.6 * (top - ha)
        hb = brentq(lambda x: needed_UA(Pa, Pb, ha, x, m, gas_hot, cap) - UA, ha + 1e-6, top, xtol=1e-6)
        nxt = pdrop(Pa, ha, hb, m, k)
        if nxt is None:
            return None
        if abs(nxt - Pb) < 1e-6:
            return nxt, hb
        Pb = nxt
    return Pb, hb


def march(md, Pd, m):
    h0 = water_steam(Pd)[1]["h"]
    m2 = m + md["leak"]
    hw = h_of(Pd + md["dp_fw"], md["T_fw"])
    mid = md["Tg"] - 45.0
    for _ in range(100):
        one = heated(Pd, h0, m, md["UA1"], md["k1"], mid, md["cap"])
        if one is None:
            return None
        hm = (m * one[1] + md["leak"] * hw) / m2
        two = heated(one[0], hm, m2, md["UA2"], md["k2"], md["Tg"], md["cap"])
        if two is None:
            return None
        nxt = md["Tg"] - m2 * (two[1] - hm) / md["cap"]
        if abs(nxt - mid) < 1e-4:
            break
        mid = nxt
    P3 = pdrop(two[0], two[1], two[1], m2, md["k3"])
    P4 = P3 and pdrop(P3, two[1], two[1], m2, md["k4"])
    if not P4:
        return None
    return dict(gap=m2 - md["C"] * math.sqrt(P4 * 1e5 / v_of(P4, two[1])), m=m,
                spans=[("SH1", Pd, one[0], h0, one[1], m), ("SH2", one[0], two[0], hm, two[1], m2),
                       ("MS", two[0], P3, two[1], two[1], m2)])


def k_smallest(sheet, md, Pd):
    def gap(m):
        st = march(md, Pd, m)
        return 1e3 if st is None else st["gap"]
    lo = 4.0
    while gap(lo * 1.2) < 0:
        lo *= 1.2
    st = march(md, Pd, brentq(gap, lo, lo * 1.2, xtol=1e-5))
    steps = {s["section"]: s for s in sheet["steam_path"]}
    worst = math.inf
    for name, Pa, Pb, ha, hb, m in st["spans"]:
        sec = steps[name]
        a, b = sec["mcr_inlet"], sec["mcr_outlet"]
        hA = water_steam(a["P_bar"])[1]["h"] if a["T_C"] is None else h_of(a["P_bar"], a["T_C"])
        hB = h_of(b["P_bar"], b["T_C"])
        for i in range(N):
            f = i / (N - 1)
            vm = IAPWS97(P=(a["P_bar"] + f * (b["P_bar"] - a["P_bar"])) / 10.0, h=hA + f * (hB - hA)).v
            worst = min(worst, (m / sec["mcr_flow_kg_s"]) ** 2 * v_of(Pa + f * (Pb - Pa), ha + f * (hb - ha)) / vm)
    return worst, st["m"]


def main():
    records, target = Path(sys.argv[1]), Path(sys.argv[2])
    answer = {}
    for folder in sorted(p for p in records.iterdir() if (p / "circuit.json").is_file()):
        sheet, rows = read(folder)
        md = fit(sheet, rows)
        lo = max(r["P_drum_bar"] for r in rows)
        hi = lo + 8.0
        while k_smallest(sheet, md, hi)[0] < AIM:
            lo, hi = hi, hi + 8.0
        hi = brentq(lambda p: k_smallest(sheet, md, p)[0] - AIM, lo, hi, xtol=5e-3)
        k, flow = k_smallest(sheet, md, hi)
        answer[folder.name + "_drum_pressure_bar"] = round(hi, 3)
        answer[folder.name + "_drum_steam_flow_kg_s"] = round(flow, 3)
        print(folder.name, round(hi, 3), round(flow, 3), round(k, 4), flush=True)
    target.write_text(json.dumps(answer, indent=1) + "\n")


if __name__ == "__main__":
    main()
