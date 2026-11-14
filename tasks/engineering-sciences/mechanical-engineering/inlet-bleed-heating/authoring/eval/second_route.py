"""A second, separately written reduction of the inlet bleed heating tests (authoring evidence only).

    python3 second_route.py <unit folder> <answer.json>

Different choices from the reference at each step: no-bleed minutes are found from the valve position signal
(shut for over five minutes) with medians, the block valve's shut minutes taking precedence; a frozen element
is one whose value repeats in most samples; for each trial line UA the energy balance is solved for the
delivered bleed at every point by bisection, and the UA kept is the one under which delivered / venturi is the
same ratio at every point; that ratio scales the venturi flow at B.
"""
import csv
import json
import math
import sys
from datetime import datetime

import numpy as np

R, K = 287.05, 1.4


def h(t):
    return 1002.5 * t + 2.75e-4 * (t - 200.0) ** 3 / 3.0


def c_p(t):
    return 1002.5 + 2.75e-4 * (t - 200.0) ** 2


def run(folder):
    unit = json.load(open(f"{folder}/unit.json"))
    rows = list(csv.DictReader(open(f"{folder}/gt.csv", newline="")))
    site = list(csv.DictReader(open(f"{folder}/site.csv", newline="")))
    log = list(csv.DictReader(open(f"{folder}/events.csv", newline="")))
    origin = datetime.fromisoformat(rows[0]["time"])
    clock = lambda s: (datetime.fromisoformat(s) - origin).total_seconds()
    t = np.array([clock(r["time"]) for r in rows])
    get = lambda c: np.array([float(r[c]) for r in rows])
    pb = 100.0 * float(np.median([float(r["BARO_HPA"]) for r in site]))

    cti = np.stack([get(f"CTI_{k:02d}_C") for k in range(1, 13)])
    for k in range(12):          # a sample over 5 K above the last good one and above the next ones is a spike
        x = cti[k]
        good = x[0]
        for i in range(1, len(x)):
            ahead = np.median(x[i + 1: i + 6]) if i + 1 < len(x) else good
            if x[i] - good > 5.0 and x[i] - ahead > 5.0:
                x[i] = good
            else:
                good = x[i]
    good = []
    for k in range(12):
        _, counts = np.unique(cti[k], return_counts=True)
        if counts.max() < 0.5 * len(t):
            good.append(k)
    mix = cti[good].mean(axis=0) + 273.15
    amb = np.stack([get(f"AMB_{k}_C") for k in range(1, 5)]).mean(axis=0) + 273.15
    diff = mix - amb

    ev = lambda name: [clock(e["time"]) for e in log if e["event"] == name]
    end_of = {e["item"]: clock(e["time"]) for e in log if e["event"] == "point_end"}
    last10 = lambda p: (t > end_of[p] - 600.0) & (t <= end_of[p])

    # no bleed: block valve shut, else the valve shut for more than five minutes (its settled part)
    calm = np.zeros_like(t, dtype=bool)
    for a, b in zip(ev("block_valve_shut"), ev("block_valve_open")):
        calm |= (t >= a + 150.0) & (t <= b)
    if not calm.any():
        shut = get("FV_POS_PCT") <= 0.0
        edges = np.flatnonzero(np.diff(np.r_[0, shut.astype(int), 0]))
        for a, b in zip(edges[0::2], edges[1::2]):
            if t[b - 1] - t[a] > 200.0:
                calm[a:b] |= t[a:b] >= t[a] + 0.5 * (t[b - 1] - t[a])
    zero = float(np.median(diff[calm])) if calm.any() else 0.0
    rise_a = float(np.mean(diff[last10("A")])) - zero

    fe = unit["flow_element"]
    beta = fe["throat_mm"] / fe["pipe_bore_mm"]
    area = math.pi / 4.0 * (fe["throat_mm"] / 1000.0) ** 2
    dp = 50.0 * (get("FE_DP_1_MBAR") + get("FE_DP_2_MBAR"))
    p1 = 1e5 * get("FE_P_BARG") + pb
    t1 = get("FE_T_C") + 273.15
    r = np.clip((p1 - dp) / p1, 0.5, 1 - 1e-12)
    eps = np.sqrt(K * r ** (2 / K) / (K - 1) * (1 - beta ** 4) / (1 - beta ** 4 * r ** (2 / K))
                  * (1 - r ** ((K - 1) / K)) / (1 - r))
    q_fe = fe["discharge_coefficient"] * eps * area * np.sqrt(2 * dp * p1 / (R * t1) / (1 - beta ** 4))
    t_in = mix - zero
    bm = unit["bellmouth"]
    q_c = bm["flow_coefficient"] * bm["area_m2"] * np.sqrt(
        2 * 100.0 * np.clip(get("BM_DP_MBAR"), 0, None) * (pb - 100.0 * get("FILTER_DP_MBAR")) / (R * t_in))

    pts = [p for p in ("A", "C", "D", "B") if p in end_of]
    avg = {p: [float(np.mean(x[last10(p)])) for x in (q_fe, q_c, t_in, amb, t1)] for p in pts}

    def delivered(p, ua):
        qf, qc, ti, ta, tf = avg[p]
        lo, hi = 1e-3, 0.5 * qc
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            tm = ta + (tf - ta) * math.exp(-ua / (mid * c_p(tf)))
            if mid * (h(tm) - h(ti)) > (qc - mid) * (h(ti) - h(ta)):
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    best = None
    for ua in np.arange(0.0, 400.0, 0.5):
        ratio = [delivered(p, ua) / avg[p][0] for p in pts]
        spread = float(np.std(ratio))
        if best is None or spread < best[0]:
            best = (spread, ua, float(np.mean(ratio)))
    return {"bleed_b_kg_s": round(best[2] * avg["B"][0], 4), "rise_a_k": round(rise_a, 4)}


if __name__ == "__main__":
    json.dump(run(sys.argv[1]), open(sys.argv[2], "w"))
