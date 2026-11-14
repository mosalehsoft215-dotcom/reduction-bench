"""A second reduction written apart from the reference.
- thrust zero: the pressurised hold after shutdown (tanks off the thrust frame), or, on a frame that carries the
  tanks, the hold before ignition plus the weight of propellant that has left, less pressurant that came in;
- fuel: the tank weight with the ullage pressurant taken out where the tank is weighed, else the flowmeter;
- chamber-pressure zero: two-point exponential through the post-shutdown excess;
- throat: At = A0 + (A1 - A0) s^q, with q chosen so that Pc At / mdot is steadiest over the mainstage;
- sums by the trapezoid rule."""
import csv
import json
import math
from pathlib import Path

import numpy as np

GM = {"He": 4.0026, "GN2": 28.0134}


def solve(folder):
    d = Path(folder)
    sh = json.loads((d / "stand_sheet.json").read_text(encoding="utf-8"))
    with open(d / "stand_data.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    keys = {k.split()[0]: k for k in rows[0]}
    get = lambda k: np.array([float(r[keys[k]]) for r in rows])
    t = get("time_s")
    ev = {}
    with open(d / "sequence.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            ev.setdefault(r["event"], float(r["time_s"]))
    press, ign, main = ev["pressurisation start"], ev["ignition command"], ev["mainstage"]
    off, vent = ev["shutdown command"], ev["vent open"]
    w = (t >= main) & (t < off)
    kind = {c["tag"]: c.get("kind") for c in sh["channels"]}
    dt = 0.02

    rho_o = 1141.0 - 4.6 * (get("TT-OX") - 90.2)
    rho_f = 810.0 - 0.75 * (get("TT-F") - 15.0)
    ox = get("FM-OX") * 1e-3 * rho_o
    fuel_meter = get("FM-F") * 1e-3 * rho_f if "FM-F" in keys else None

    lc = get("LC-1") * 1e3
    on_frame = "thrust frame" in sh.get("stand", "")
    if on_frame:
        base = lc[(t > press + 6) & (t < ign - 2)].mean()
        cum = np.zeros(t.size)
        for flow, rho_l, pt, tu, tank in ((ox, rho_o, "PT-OX", "TT-OXU", "LOX"), (fuel_meter, rho_f, "PT-F", "TT-FU", "RP-1")):
            rho_g = get(pt) * 1e3 * GM[sh["tanks"][tank]["pressurant"]] / (8314.462 * get(tu))
            net = flow * (1 - rho_g / rho_l)
            cum += np.concatenate([[0.0], np.cumsum(0.5 * (net[1:] + net[:-1]) * dt)])
        cum -= cum[np.searchsorted(t, ign)]
        force = lc - base - 9.80665 * np.where(t >= ign, cum, 0.0)
    else:
        force = lc - lc[(t > off + 5) & (t < vent - 2)].mean()

    amb = get("PB-1") * 1e3
    pc = get("PC-1") * 1e3 + (amb if kind["PC-1"] == "gauge" else 0.0)
    ex = pc - amb
    ma, mb = (t > off + 2) & (t < off + 4), (t > vent - 4) & (t < vent - 2)
    ea, eb = ex[ma].mean(), ex[mb].mean()
    if ea > 5e3 and eb > 1e3:
        k = math.log(ea / eb) / (t[mb].mean() - t[ma].mean())
        pc = pc - np.where(t >= ign, ea * math.exp(k * (t[ma].mean() - off)), 0.0)

    if fuel_meter is not None and "WT-F" not in keys:
        fuel_kg = np.trapezoid(fuel_meter[w], dx=dt)
        fuel_rate = fuel_meter
    elif "WT-F" in keys:
        edge = lambda x, at: x[(t >= at - 0.5) & (t < at + 0.5)].mean()
        rho_g = (get("PT-F") * 1e3 * GM[sh["tanks"]["RP-1"]["pressurant"]] / (8314.462 * get("TT-FU")))[w].mean()
        fuel_kg = (edge(get("WT-F"), main) - edge(get("WT-F"), off)) / (1 - rho_g / rho_f[w].mean())
        fuel_rate = fuel_meter if fuel_meter is not None else np.full(t.size, fuel_kg / (off - main))
    total = np.trapezoid(ox[w], dx=dt) + fuel_kg
    mdot = ox + fuel_rate

    e = sh["engine"]
    a0 = math.pi * (e["throat_diameter_mm_before"] / 2000) ** 2
    a1 = math.pi * (e["throat_diameter_mm_after"] / 2000) ** 2
    s = np.clip((t - ign) / (off - ign), 0, 1)
    if abs(a1 - a0) > 1e-9:
        best = None
        for q in np.linspace(0.2, 1.6, 141):
            at = a0 + (a1 - a0) * s ** q
            keep = w & (s > 0.05) & (mdot > 0)
            c = (pc * at)[keep] / mdot[keep]
            spread = np.std(c) / np.mean(c)
            if best is None or spread < best[0]:
                best = (spread, q)
        throat = a0 + (a1 - a0) * s ** best[1]
    else:
        throat = np.full(t.size, a0)
    return {"isp_s": float(np.trapezoid(force[w], dx=dt) / (9.80665 * total)),
            "cstar_m_s": float(np.trapezoid((pc * throat)[w], dx=dt) / total)}


reduce = solve
