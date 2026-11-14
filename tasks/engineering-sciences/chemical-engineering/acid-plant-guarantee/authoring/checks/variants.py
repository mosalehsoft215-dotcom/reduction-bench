"""The reference reduction with single steps swapped for the readings an analyst might take instead."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "solution"))
import acid_test as A

OPTS = {}


def reduce(plant_dir, **o):
    import json
    d = Path(plant_dir)
    plant = json.loads((d / "plant.json").read_text(encoding="utf-8"))
    t0, t, c = A.read_trend(d / "trend.csv")
    log = A.read_log(d / "log.txt", t0)
    cyl = A.read_cylinders(d / "cylinders.csv")
    inst = {i["tag"]: i for i in plant["instruments"]}
    cors = A.corrections(log, cyl)
    win = [e["t"] for e in log if e["words"][:1] == ["TEST"]]
    w0, w1 = min(win), max(win)
    inw = (t >= w0) & (t < w1)

    def analyzer(tag):
        x = cors[tag].true(c[tag], t) if (tag in cors and not o.get("no_cal")) else c[tag].copy()
        return A.fill(x, A.invalid_mask(t, log, {tag}))

    n_s = c["FT-101"] * plant["sulfur"]["S_mass_share"] / A.M_S
    ven = plant["air_venturi"]
    rho1 = c["PT-201"] * A.M_AIR / (A.R * (c["TT-201"] + 273.15))
    F_ven = ven["C_kg_h"] * np.sqrt(np.maximum(c["PDT-201"], 0.0) * rho1) / A.M_AIR
    stops = [e["t"] for e in log if "sulfur-feed-stopped" in e["words"]]
    starts = [e["t"] for e in log if "sulfur-feed-restarted" in e["words"]]
    trips = list(zip(stops, starts))
    steady = np.ones(t.size, bool)
    for a, b in trips:
        steady &= ~((t >= a - 120) & (t < b + 1200))
    for e in log:
        if "setpoint-change" in e["words"]:
            steady &= ~((t >= e["t"]) & (t < e["t"] + 1500))
    y_o2 = analyzer("AT-302") / 100.0
    so2_tags = [k for k in ("AT-301A", "AT-301B") if k in c]
    so2 = {k: analyzer(k) / 100.0 for k in so2_tags}
    sel = steady & inw
    for k in so2_tags:
        sel &= np.isfinite(so2[k])
    ys = np.nanmean([np.median(so2[k][sel]) for k in so2_tags])
    f = o.get("f", 0.01)
    F_o2 = n_s * (1 + 0.5 * f - 1.5 * f * y_o2) / (A.O2_AIR - y_o2)
    route = o.get("flow", "o2")
    if route == "o2":
        F = np.median(F_o2[sel] / F_ven[sel]) * F_ven
    elif route == "venturi":
        F = F_ven
    elif route == "composition":
        F_c = n_s * (1 - f) / ys + 1.5 * f * n_s
        F = np.median(F_c[sel] / F_ven[sel]) * F_ven
    elif route == "dcs":
        F = c["FI-201"] / 22.414
    tau = 40.0
    a = 5.0 / (tau + 5.0)
    n_eff = np.empty_like(n_s)
    acc = n_s[0]
    for i, v in enumerate(n_s):
        acc += a * (v - acc)
        n_eff[i] = acc
    y_ipa = analyzer("AT-401") / 100.0
    scr = A.scrubbed(inst["AT-401"])
    if o.get("basis") == "insitu":
        scr = False
    if o.get("basis") == "scrubbed":
        scr = True
    if o.get("ipa") == "ratio":
        # textbook formula from inlet and outlet concentrations
        b = np.median(y_ipa[sel]); aa = ys
        kf = 1.5 if scr else 0.5
        x_ipa = (aa - b) / (aa * (1 - kf * b))
    else:
        x_guess = 0.95
        for _ in range(5):
            q_ipa = y_ipa * (F - (1.5 if scr else 0.5) * n_eff * x_guess)
            x_guess = 1.0 - q_ipa[inw].sum() / n_eff[inw].sum()
        x_ipa = x_guess
    y_st = analyzer("AT-501") * 1e-6
    yo_st = analyzer("AT-502") / 100.0
    if trips and not o.get("no_floor"):
        floor = []
        for a0, b0 in trips:
            m = (t >= a0 + 300) & (t < b0 + 60) & (yo_st > 0.20)
            if m.sum() > 10:
                floor.append(np.median(y_st[m]))
        if floor:
            y_st = y_st - np.median(floor)
    keep = inw.copy()
    if o.get("cut_trip"):
        for a0, b0 in trips:
            keep &= ~((t >= a0 - 120) & (t < b0 + 1200))
    x_tot = 0.998
    for _ in range(5):
        N_st = F - 1.5 * n_eff * x_tot
        o2_exp = (A.O2_AIR * F - n_eff * (1 + 0.5 * x_tot)) / N_st
        ms = steady & inw & (n_eff > 0.5 * np.median(n_eff[inw]))
        leak = (np.median(yo_st[ms]) - np.median(o2_exp[ms])) / (A.O2_AIR - np.median(o2_exp[ms]))
        leak = leak if (leak > 0.01 and not o.get("no_leak")) else 0.0
        if o.get("stack") == "cems_flow":
            N_st = c["FT-501"] * 273.15 / (c["TT-501"] + 273.15) * c["PT-501"] / 101.325 / 22.414
        elif o.get("stack") == "inert":
            N_st = F * (1 - A.O2_AIR) / (1 - yo_st - y_st)
        q_st = y_st / (1.0 - leak) * N_st
        x_tot = 1.0 - q_st[keep].sum() / n_eff[keep].sum()
    acid = (n_eff * x_tot)[keep].sum()
    so2_kg_t = 1000.0 * A.M_SO2 * q_st[keep].sum() / (A.M_H2SO4 * acid)
    return {"ipa_conversion_pct": 100.0 * x_ipa, "so2_kg_per_t": float(so2_kg_t)}
