#!/usr/bin/env python3
"""Reduce one performance test of a sulfur-burning double-absorption acid plant.

usage: acid_test.py PLANT_DIR OUT_JSON

Writes {"ipa_conversion_pct": ..., "so2_kg_per_t": ...} for the logged test window.
"""
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

O2_AIR = 0.20946
M_S, M_SO2, M_H2SO4, M_AIR = 32.06, 64.066, 98.079, 28.9647
R = 8.314462


# ---------------------------------------------------------------- reading
def read_trend(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    t0 = datetime.fromisoformat(rows[1][0])
    t = np.array([(datetime.fromisoformat(r[0]) - t0).total_seconds() for r in rows[1:]])
    cols = {}
    for j, name in enumerate(head[1:], start=1):
        cols[name] = np.array([float(r[j]) if r[j] != "" else np.nan for r in rows[1:]])
    return t0, t, cols


def read_log(path, t0):
    out = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = [p.strip() for p in line.split("|")]
        s = (datetime.fromisoformat(parts[0]) - t0).total_seconds()
        kv = {}
        words = []
        for p in parts[1:]:
            if "=" in p:
                k, v = p.split("=", 1)
                kv[k] = v
            else:
                words.append(p)
        out.append({"t": s, "words": words, "kv": kv})
    return out


def read_cylinders(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["cylinder"]: float(r["certified"]) for r in csv.DictReader(fh)}


# ---------------------------------------------------------------- analyzers
class Correction:
    """Zero and gain of an analyzer as functions of time, from its logged adjustments and checks."""

    def __init__(self):
        self.zero_pts, self.gain_pts = [], []

    def value(self, pts, t, default):
        if not pts:
            return np.full_like(t, default)
        pts = sorted(pts)
        ts = np.array([p[0] for p in pts])
        vs = np.array([p[1] for p in pts])
        return np.interp(t, ts, vs)

    def true(self, reading, t):
        z = self.value(self.zero_pts, t, 0.0)
        g = self.value(self.gain_pts, t, 1.0)
        return (reading - z) / g


def corrections(log, cyl):
    cors = {}
    for e in log:
        if not e["words"] or e["words"][0] != "CAL":
            continue
        tag = e["words"][1]
        kv = e["kv"]
        c = cors.setdefault(tag, Correction())
        kind = kv["kind"]
        gas = kv["gas"]
        if kind.endswith("adjust"):
            if kind.startswith("zero"):
                c.zero_pts.append((e["t"], 0.0))
            else:
                c.gain_pts.append((e["t"], 1.0))
            continue
        found = float(kv["as_found"])
        if kind == "zero-check":
            c.zero_pts.append((e["t"], found))
        elif kind == "span-check":
            ref = 100 * O2_AIR if gas == "ambient-air" else cyl[gas]
            z = c.value(c.zero_pts, np.array([e["t"]]), 0.0)[0]
            c.gain_pts.append((e["t"], (found - z) / ref))
    return cors


def invalid_mask(t, log, tags):
    """Samples of these analyzers that show calibration gas or a refilling line."""
    bad = np.zeros(t.size, bool)
    begin = None
    for e in log:
        w = e["words"]
        if not w:
            continue
        if w[0] == "CAL" and w[1] in tags:
            bad |= (t >= e["t"] - 60) & (t < e["t"] + 240)
        if w[0] == "CAL-PERIOD" and any(x in w[1].split(",") for x in tags):
            if w[2] == "begin":
                begin = e["t"]
            else:
                bad |= (t >= begin - 60) & (t < e["t"] + 180)
    return bad


def fill(x, bad):
    y = x.copy()
    good = ~bad & np.isfinite(x)
    idx = np.arange(x.size)
    y[~good] = np.interp(idx[~good], idx[good], x[good])
    return y


def scrubbed(inst):
    return "acid scrubber" in inst.get("sampling", "")


# ---------------------------------------------------------------- reduction
def reduce(plant_dir):
    d = Path(plant_dir)
    plant = json.loads((d / "plant.json").read_text(encoding="utf-8"))
    t0, t, c = read_trend(d / "trend.csv")
    log = read_log(d / "log.txt", t0)
    cyl = read_cylinders(d / "cylinders.csv")
    inst = {i["tag"]: i for i in plant["instruments"]}
    cors = corrections(log, cyl)
    win = [e["t"] for e in log if e["words"][:1] == ["TEST"]]
    w0, w1 = min(win), max(win)
    inw = (t >= w0) & (t < w1)

    def analyzer(tag):
        x = cors[tag].true(c[tag], t) if tag in cors else c[tag].copy()
        return fill(x, invalid_mask(t, log, {tag}))

    w_s = plant["sulfur"]["S_mass_share"]
    n_s = c["FT-101"] * w_s / M_S  # kmol/h
    ven = plant["air_venturi"]
    rho1 = c["PT-201"] * M_AIR / (R * (c["TT-201"] + 273.15))
    F_ven = ven["C_kg_h"] * np.sqrt(np.maximum(c["PDT-201"], 0.0) * rho1) / M_AIR

    # trip periods: sulfur feed stopped
    stops = [e["t"] for e in log if "sulfur-feed-stopped" in e["words"]]
    starts = [e["t"] for e in log if "sulfur-feed-restarted" in e["words"]]
    trips = list(zip(stops, starts))
    steady = np.ones(t.size, bool)
    for a, b in trips:
        steady &= ~((t >= a - 120) & (t < b + 1200))
    for e in log:
        if "setpoint-change" in e["words"]:
            steady &= ~((t >= e["t"]) & (t < e["t"] + 1500))

    # converter-inlet O2 and SO2, dry and SO3-free
    y_o2 = analyzer("AT-302") / 100.0
    so2_tags = [k for k in ("AT-301A", "AT-301B") if k in c]
    so2 = {k: analyzer(k) / 100.0 for k in so2_tags}
    sel = steady & inw
    for k in so2_tags:
        sel &= np.isfinite(so2[k])

    def F_from_o2(f):
        return n_s * (1 + 0.5 * f - 1.5 * f * y_o2) / (O2_AIR - y_o2)

    # furnace SO3 share from O2 + SO2 where the pair is physical, else a typical 1 %
    f = 0.01
    ys = np.nanmean([np.median(so2[k][sel]) for k in so2_tags])
    yo = np.median(y_o2[sel])
    # y_s (F - 1.5 f n) = (1 - f) n and y_o (F - 1.5 f n) = O2_AIR F - n (1 + 0.5 f)
    best = None
    for ff in np.linspace(0.0, 0.03, 301):
        Fn = (1 + 0.5 * ff - 1.5 * ff * yo) / (O2_AIR - yo)
        ys_pred = (1 - ff) / (Fn - 1.5 * ff)
        err = abs(ys_pred - ys)
        if best is None or err < best[0]:
            best = (err, ff)
    if best[0] / ys < 0.004:
        f = best[1]
    F_o2 = F_from_o2(f)
    # carry the O2-based flow through transients with the venturi, scaled on steady data
    k = np.median(F_o2[sel] / F_ven[sel])
    F = k * F_ven

    # converter gas holdup for the sulfur that actually burns
    tau = 40.0
    a = 5.0 / (tau + 5.0)
    n_eff = np.empty_like(n_s)
    acc = n_s[0]
    for i, v in enumerate(n_s):
        acc += a * (v - acc)
        n_eff[i] = acc

    # interpass absorber inlet
    y_ipa = analyzer("AT-401") / 100.0
    x_guess = 0.95
    for _ in range(5):
        if scrubbed(inst["AT-401"]):
            q_ipa = y_ipa * (F - 1.5 * n_eff * x_guess)
        else:
            q_ipa = y_ipa * (F - 0.5 * n_eff * x_guess)
        x_guess = 1.0 - q_ipa[inw].sum() / n_eff[inw].sum()
    x_ipa = x_guess

    # stack
    y_st = analyzer("AT-501") * 1e-6
    yo_st = analyzer("AT-502") / 100.0
    if trips:
        floor = []
        for a0, b0 in trips:
            m = (t >= a0 + 300) & (t < b0 + 60) & (yo_st > 0.20)
            if m.sum() > 10:
                floor.append(np.median(y_st[m]))
        if floor:
            y_st = y_st - np.median(floor)
    x_tot = 0.998
    for _ in range(5):
        N_st = F - 1.5 * n_eff * x_tot
        o2_exp = (O2_AIR * F - n_eff * (1 + 0.5 * x_tot)) / N_st
        ms = steady & inw & (n_eff > 0.5 * np.median(n_eff[inw]))
        leak = (np.median(yo_st[ms]) - np.median(o2_exp[ms])) / (O2_AIR - np.median(o2_exp[ms]))
        leak = leak if leak > 0.01 else 0.0
        q_st = y_st / (1.0 - leak) * N_st
        x_tot = 1.0 - q_st[inw].sum() / n_eff[inw].sum()
    acid = (n_eff * x_tot)[inw].sum()
    so2_kg_t = 1000.0 * M_SO2 * q_st[inw].sum() / (M_H2SO4 * acid)
    return {"ipa_conversion_pct": round(100.0 * x_ipa, 4), "so2_kg_per_t": round(float(so2_kg_t), 5)}


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: acid_test.py PLANT_DIR OUT_JSON")
    res = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(res) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
