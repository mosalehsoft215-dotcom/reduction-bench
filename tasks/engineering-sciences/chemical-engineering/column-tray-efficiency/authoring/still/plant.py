"""Writes the records of one column test run: column_sheet.json, dcs_trend.csv, lab_results.csv, shift_log.txt.

Each run is a series of test points. The flows named in a point entry are held by flow controllers, each of which
keeps its own meter's indication at the set point, so the true flow is that indication divided by the meter's own
factor (reading = factor x true flow). The steady state of a point is the one solve_column gives at the true flows
with the run's true tray efficiencies; compositions and temperatures move from point to point as one first-order lag.
"""
import csv
import json
from datetime import datetime, timedelta

import numpy as np
from scipy.optimize import brentq

from column import Mixture, solve_column

SYSTEMS = {
    "methanol-water": {
        "light": {"name": "methanol", "molar_mass_g_mol": 32.042, "antoine": {"A": 8.08097, "B": 1582.271, "C": 239.726},
                  "cp_liquid_J_molK": 81.1, "cp_vapour_J_molK": 44.1, "dH_vap_nbp_J_mol": 35270.0,
                  "density_20C_kg_m3": 791.8, "expansion_1_K": 0.00119},
        "heavy": {"name": "water", "molar_mass_g_mol": 18.015, "antoine": {"A": 8.07131, "B": 1730.63, "C": 233.426},
                  "cp_liquid_J_molK": 75.3, "cp_vapour_J_molK": 33.6, "dH_vap_nbp_J_mol": 40650.0,
                  "density_20C_kg_m3": 998.2, "expansion_1_K": 0.00045},
        "nrtl": {"a12": 0.0, "b12_K": -127.76, "a21": 0.0, "b21_K": 425.33, "alpha": 0.2994}},
    "ethanol-water": {
        "light": {"name": "ethanol", "molar_mass_g_mol": 46.069, "antoine": {"A": 8.20417, "B": 1642.89, "C": 230.300},
                  "cp_liquid_J_molK": 112.4, "cp_vapour_J_molK": 65.6, "dH_vap_nbp_J_mol": 38560.0,
                  "density_20C_kg_m3": 789.3, "expansion_1_K": 0.00109},
        "heavy": {"name": "water", "molar_mass_g_mol": 18.015, "antoine": {"A": 8.07131, "B": 1730.63, "C": 233.426},
                  "cp_liquid_J_molK": 75.3, "cp_vapour_J_molK": 33.6, "dH_vap_nbp_J_mol": 40650.0,
                  "density_20C_kg_m3": 998.2, "expansion_1_K": 0.00045},
        "nrtl": {"a12": -0.9852, "b12_K": 302.2365, "a21": 3.7555, "b21_K": -676.0314, "alpha": 0.3}},
}
# the pressure column runs at 100-140 C, above the range of the low-temperature Antoine sets
SYSTEMS["methanol-water, high-temperature Antoine"] = dict(
    SYSTEMS["methanol-water"],
    light=dict(SYSTEMS["methanol-water"]["light"], antoine={"A": 7.97010, "B": 1521.23, "C": 233.970}),
    heavy=dict(SYSTEMS["methanol-water"]["heavy"], antoine={"A": 8.14019, "B": 1810.94, "C": 244.485}))

# heating steam: saturation temperature C, pressure bar a, latent heat J/kg
STEAM = {"4bar": (143.6, 4.0, 2.133e6), "6bar": (158.8, 6.0, 2.0858e6)}
FLOW_LAG_MIN = 14.0


def flashed_sample(mix, x, T, P_amb):
    """Tray liquid (bubble point at tray pressure) let down through an open valve to P_amb with no heat exchange:
    it flashes to equilibrium, the vapour escapes, and the liquid left is what the laboratory analyses."""
    Tb, _ = mix.bubble(np.array(x), np.array(P_amb))
    if T <= float(Tb):
        return x
    h0 = float(mix.hL(x, T))
    x_dew, _ = mix.dew_liquid(x, P_amb)

    def excess(xl):
        Tl, yl = mix.bubble(np.array(xl), np.array(P_amb))
        phi = (x - xl) / (float(yl) - xl)
        return (1 - phi) * float(mix.hL(xl, Tl)) + phi * float(mix.hV(yl, Tl)) - h0
    return brentq(excess, x_dew + 1e-12, x, xtol=1e-14, rtol=1e-13)


def density(mixd, x1, T):
    """Liquid density, ideal mixing of volumes, each pure liquid expanding linearly from 20 C."""
    v = 0.0
    for c, x in ((mixd["light"], x1), (mixd["heavy"], 1 - x1)):
        rho = c["density_20C_kg_m3"] * (1 - c["expansion_1_K"] * (T - 293.15))
        v = v + x * c["molar_mass_g_mol"] / 1000 / rho
    M = x1 * mixd["light"]["molar_mass_g_mol"] / 1000 + (1 - x1) * mixd["heavy"]["molar_mass_g_mol"] / 1000
    return M / v


def heat_per_kg_steam(spec):
    """J that one kg of heating steam gives to the kettle (its latent heat, plus the sensible heat of a condensate
    that leaves subcooled)."""
    Tsat, _, lam = STEAM[spec.get("steam", "4bar")]
    sub = spec.get("condensate_C")
    return lam + (0.0 if sub is None else spec["condensate_cp_J_kgK"] * (Tsat - sub))


def steady(mix, col, pt, E, guess=None, mixd=None):
    """One test point at its true flows (kg/h; reflux_m3_h for a reflux given as a volume flow), C, Pa.

    Iterates the molar conversions of the flows, which depend on the compositions of the streams."""
    zF = pt["zF"]
    F = pt["feed_kg_h"] / 3600 / mix.molar_mass(zF)
    c = dict(col, feed_tray=pt.get("feed_tray", col["feed_tray"]))
    xD, x0, xs = 0.9, 0.9, 0.5
    r = None
    for _ in range(25):
        TR = pt["TR_C"] + 273.15
        if "reflux_m3_h" in pt:
            R_kg_s = pt["reflux_m3_h"] * density(mixd, x0, TR) / 3600
        else:
            R_kg_s = pt["reflux_kg_h"] / 3600
        op = dict(F=F, zF=zF, hF=mix.hL(zF, pt["TF_C"] + 273.15), L0=R_kg_s / mix.molar_mass(x0), TR=TR,
                  D=pt["dist_kg_h"] / 3600 / mix.molar_mass(xD), P_top=pt["P_top_Pa"], T_amb=pt["T_amb_C"] + 273.15,
                  S=pt.get("side_kg_h", 0.0) / 3600 / mix.molar_mass(xs))
        r = solve_column(mix, c, op, *E, guess=guess if r is None else r["z"])
        if not r["ok"]:
            raise RuntimeError(f"column solve failed, residual {r['resid']:.2e}")
        new = (r["xD"], r["x0"], r["x"][c["side_tray"] - 1] if c.get("side_tray") else 0.5)
        if max(abs(new[0] - xD), abs(new[1] - x0), abs(new[2] - xs)) < 1e-12:
            break
        xD, x0, xs = new
    r["op"] = op
    r["R_kg_h"] = R_kg_s * 3600
    return r


def steady_on_steam(mix, col, pt, E, steam_true_kg_h, q_per_kg, loss_frac, mixd, guess=None):
    """Reflux on drum level, reboiler steam on flow control: the reflux is whatever the steam boils up."""
    target = steam_true_kg_h * q_per_kg / 3600 / (1 + loss_frac)       # W taken up by the process side

    def gap(R):
        return steady(mix, col, dict(pt, reflux_kg_h=R), E, guess=guess, mixd=mixd)["Q_reb"] - target
    R = brentq(gap, 0.2 * pt["dist_kg_h"], 6.0 * pt["dist_kg_h"], xtol=1e-10, rtol=1e-13)
    return steady(mix, col, dict(pt, reflux_kg_h=R), E, guess=guess, mixd=mixd)


def preheated(mix, col, pt, E, ph, mixd, guess=None):
    """Feed heated against the bottoms in a counter-current exchanger of effectiveness ph['eff'] (no losses): the
    feed temperature at the nozzle and the bottoms outlet follow from the exchanger and the column together."""
    TF = pt["TF_C"]
    for _ in range(80):
        r = steady(mix, col, dict(pt, TF_C=TF), E, guess=guess, mixd=mixd)
        cF = r["op"]["F"] * (pt["zF"] * mix.cpl[0] + (1 - pt["zF"]) * mix.cpl[1])
        cB = r["B"] * (r["xB"] * mix.cpl[0] + (1 - r["xB"]) * mix.cpl[1])
        TB = r["T"][-1] - 273.15
        Q = ph["eff"] * min(cF, cB) * (TB - ph["T_in_C"])
        TF_new = ph["T_in_C"] + Q / cF
        if abs(TF_new - TF) < 1e-10:
            break
        TF = TF_new
    r["TF_C"] = TF
    r["T_bot_out_C"] = TB - Q / cB
    return r


def write_run(folder, spec, rng):
    folder.mkdir(parents=True, exist_ok=True)
    mixd = SYSTEMS[spec["system"]]
    mix = Mixture(mixd)
    col = spec["column"]
    col = dict(col, heat_loss_W_K=np.array(col["heat_loss_W_K"], dtype=float))
    E = (spec["E_rect"], spec["E_strip"])
    inst = spec["instruments"]
    t0 = datetime.fromisoformat(spec["start"])
    streams = ("feed", "reflux", "dist", "bottoms", "side", "steam")
    factor = {k: 1 + rng.normal(0, 0.001) for k in streams}       # each meter: reading = factor x true flow
    for k, f in spec.get("as_found", {}).items():                  # a meter found off calibration after the run
        factor[k] *= f
    controlled = set(spec.get("controlled", ("feed", "reflux", "dist", "side")))
    out_of_service = set(spec.get("meters_out", ()))
    vol_reflux = inst["reflux_meter"] == "volumetric"
    lineout = spec.get("lineout_min", 70)
    window = 120
    sample_min = spec.get("sample_minutes", (30, 90))
    tau_c = spec.get("tau_comp_min", 14.0)
    tc_out = spec.get("tc_out", False)
    ph = spec.get("preheater")
    tc_bias = {j: rng.normal(0, 0.08) for j in list(col["tray_TIs"]) + ["reb", "fin", "bout"]}
    q_kg = heat_per_kg_steam(spec)
    loss = spec["reboiler_loss_frac"]

    states = []
    guess = None
    for pt in spec["points"]:
        sp = dict(pt["setpoints"])
        true = {key: sp[key] / factor[key] for key in ("feed", "dist", "side") if key in sp}
        tp = dict(zF=pt["zF"], feed_kg_h=true["feed"], dist_kg_h=true["dist"], TF_C=pt.get("TF_C", 0.0),
                  TR_C=pt["TR_C"], P_top_Pa=pt["P_top_Pa"], T_amb_C=pt["T_amb_C"], feed_tray=pt["feed_tray"])
        if "side" in true:
            tp["side_kg_h"] = true["side"]
        if "steam" in sp:
            r = steady_on_steam(mix, col, tp, E, sp["steam"] / factor["steam"], q_kg, loss, mixd, guess)
        else:
            if vol_reflux:
                tp["reflux_m3_h"] = sp["reflux"] / factor["reflux"]
            else:
                tp["reflux_kg_h"] = sp["reflux"] / factor["reflux"]
            if ph:
                r = preheated(mix, col, tp, E, ph, mixd, guess)
                tp["TF_C"] = r["TF_C"]
            else:
                r = steady(mix, col, tp, E, guess, mixd)
        guess = r["z"]
        true_flow = {"feed": tp["feed_kg_h"], "dist": tp["dist_kg_h"], "side": tp.get("side_kg_h", 0.0),
                     "bottoms": r["B"] * mix.molar_mass(r["xB"]) * 3600,
                     "steam": r["Q_reb"] * (1 + loss) / q_kg * 3600,
                     "reflux": r["R_kg_h"] / density(mixd, r["x0"], tp["TR_C"] + 273.15) if vol_reflux else r["R_kg_h"]}
        shown = {key: (sp[key] if (key in sp and key in controlled) else true_flow[key] * factor[key]) for key in streams}
        temps = {j: r["T"][j - 1] - 273.15 for j in col["tray_TIs"]}
        temps["reb"] = r["T"][-1] - 273.15
        raw = {"distillate": (float(r["xD"]), None), "bottoms": (float(r["xB"]), None)}
        raw.update({f"tray {j}": (float(r["x"][j - 1]), float(r["T"][j - 1])) for j in col["sample_trays"]})
        if col.get("side_tray"):
            raw["side draw"] = (float(r["x"][col["side_tray"] - 1]), None)
        states.append(dict(shown=shown, temps=temps, TF=tp["TF_C"], TR=r["TR"] - 273.15,
                           P=pt["P_top_Pa"] - (col["dP_condenser_Pa"] if inst["pressure_tap"] == "reflux drum" else 0.0),
                           Tamb=pt["T_amb_C"], r=r, tp=tp, sp=sp, true_flow=true_flow, raw=raw,
                           Tbout=r.get("T_bot_out_C"), Tfin=ph["T_in_C"] if ph else None))

    # one first-order response: from its value at a move, each quantity approaches the new point's steady value
    period = lineout + window

    def lagged(values, since, tau):
        i = min(int(since // period), len(values) - 1)
        if i == 0:
            return values[0]
        s = values[0]
        for j in range(1, i):
            s = values[j] + (s - values[j]) * np.exp(-period / tau)
        return values[i] + (s - values[i]) * np.exp(-(since - i * period) / tau)

    trend, lab, log = [], [], []
    log.append((t0, f"Test run started by {spec['crew']}. Column running at the settings of point 1 since the night shift."))
    sample_no = 1
    tags = inst["tags"]
    flow_keys = [k for k in streams if k not in out_of_service and (k != "side" or col.get("side_tray"))]
    for k, cur in enumerate(states, start=1):
        sp = cur["sp"]
        if "steam" in sp:
            mid = f"steam {sp['steam']:.1f} kg/h"
        elif vol_reflux:
            mid = f"reflux {sp['reflux']:.4f} m3/h"
        else:
            mid = f"reflux {sp['reflux']:.1f} kg/h"
        note = f"Point {k}: feed {sp['feed']:.0f} kg/h on tray {cur['tp']['feed_tray']}, {mid}, distillate {sp['dist']:.1f} kg/h"
        if col.get("side_tray"):
            note += f", side draw {sp['side']:.1f} kg/h"
        move = t0 + timedelta(minutes=period * (k - 1))
        log.append((move, note + "."))
        for m in range(period):
            since = (k - 1) * period + m
            row = {"time": move + timedelta(minutes=m)}
            for key in flow_keys:
                val = lagged([s["shown"][key] for s in states], since, FLOW_LAG_MIN)
                row[key] = val * (1 + rng.normal(0, 0.004)) if val > 0 else 0.0
            for j in list(col["tray_TIs"]) + ["reb"]:
                row[f"T{j}"] = lagged([s["temps"][j] for s in states], since, tau_c) + tc_bias[j] + rng.normal(0, 0.05)
            row["TF"] = lagged([s["TF"] for s in states], since, tau_c if ph else FLOW_LAG_MIN) + rng.normal(0, 0.05)
            row["TR"] = lagged([s["TR"] for s in states], since, tau_c) + rng.normal(0, 0.05) + spec.get("tr_offset_K", 0.0)
            row["P"] = lagged([s["P"] for s in states], since, FLOW_LAG_MIN) + rng.normal(0, 30.0)
            row["Tamb"] = lagged([s["Tamb"] for s in states], since, FLOW_LAG_MIN) + rng.normal(0, 0.1)
            if ph:
                row["Tfin"] = cur["Tfin"] + tc_bias["fin"] + rng.normal(0, 0.05)
                row["Tbout"] = lagged([s["Tbout"] for s in states], since, tau_c) + tc_bias["bout"] + rng.normal(0, 0.05)
            trend.append(row)
        opened = move + timedelta(minutes=lineout)
        log.append((opened, f"Point {k} lined out; test window opened."))
        for minute in sample_min:
            ts = opened + timedelta(minutes=minute)
            since = (k - 1) * period + lineout + minute
            points = [("feed", states[0]["tp"]["zF"])]
            for name, (_, Tv0) in cur["raw"].items():
                xv = lagged([s["raw"][name][0] for s in states], since, tau_c)
                if Tv0 is not None and spec.get("open_tray_valves"):
                    Tv = lagged([s["raw"][name][1] for s in states], since, tau_c)
                    xv = flashed_sample(mix, xv, Tv, spec["barometric_Pa"])
                points.append((name, xv))
            for name, x in points:
                xm = float(np.clip(x + rng.normal(0, 0.0003 + 0.004 * min(x, 1 - x) + 0.001 * x), 1e-6, 1 - 1e-6))
                if inst["lab_basis"] == "mass":
                    M1, M2 = mixd["light"]["molar_mass_g_mol"], mixd["heavy"]["molar_mass_g_mol"]
                    val = 100 * xm * M1 / (xm * M1 + (1 - xm) * M2)
                else:
                    val = 100 * xm
                lab.append((f"S{sample_no:03d}", ts, name, val, ts + timedelta(hours=3 + int(rng.integers(0, 5)))))
                sample_no += 1
        log.append((opened + timedelta(minutes=window), f"Point {k} test window closed."))
    for extra in spec.get("log_extra", []):
        log.append((t0 + timedelta(minutes=extra[0]), extra[1]))
    for chk in spec.get("meter_checks", []):
        # collection test after the run, column held at one point's settings: the meter's steady indication against
        # what a weighed tank or pot took up (the true mass flow over the stated minutes)
        st = states[chk["point"] - 1]
        key = chk["key"]
        mass_h = st["r"]["R_kg_h"] if key == "reflux" else st["true_flow"][key]
        log.append((t0 + timedelta(minutes=chk["at_min"]),
                    chk["text"].format(ind=st["shown"][key], mass=mass_h * chk["minutes"] / 60, T=st["TR"])))
    log.sort(key=lambda e: e[0])

    head, keys = ["timestamp", tags["feed"]], ["feed"]
    if not ph:
        head.append(tags["TF"]); keys.append("TF")
    if "reflux" in flow_keys:
        head.append(tags["reflux"]); keys.append("reflux")
    head.append(tags["TR"]); keys.append("TR")
    for key in ("dist", "bottoms", "side", "steam"):
        if key in flow_keys:
            head.append(tags[key]); keys.append(key)
    head.append(tags["P"]); keys.append("P")
    if not tc_out:
        head += [f"TI{100 + j}_tray{j}_C" for j in col["tray_TIs"]] + ["TI199_reboiler_C"]
        keys += [f"T{j}" for j in col["tray_TIs"]] + ["Treb"]
    if ph:
        head += [ph["tag_in"], ph["tag_out"]]; keys += ["Tfin", "Tbout"]
    head.append("TI001_ambient_C"); keys.append("Tamb")
    with open(folder / "dcs_trend.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(head)
        for row in trend:
            line = [row["time"].strftime("%Y-%m-%d %H:%M")]
            for key in keys:
                if key == "P":
                    line.append(f"{row['P'] / 1000 if inst['pressure_unit'] == 'kPa' else row['P'] / 100:.2f}")
                elif key == "reflux" and vol_reflux:
                    line.append(f"{row['reflux']:.4f}")
                else:
                    line.append(f"{row[key]:.2f}")
            w.writerow(line)
    with open(folder / "lab_results.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        unit = "wt_pct" if inst["lab_basis"] == "mass" else "mol_pct"
        w.writerow(["sample", "taken", "point", f"{mixd['light']['name']}_{unit}", "reported"])
        for s in lab:
            w.writerow([s[0], s[1].strftime("%Y-%m-%d %H:%M"), s[2], f"{s[3]:.3f}", s[4].strftime("%Y-%m-%d %H:%M")])
    with open(folder / "shift_log.txt", "w", newline="\n") as fh:
        for tm, text in log:
            fh.write(f"{tm.strftime('%d.%m.%Y %H:%M')} | {text}\n")
    sheet = spec["sheet"](mixd, col)
    (folder / "column_sheet.json").write_text(json.dumps(sheet, indent=2) + "\n", encoding="ascii", newline="\n")
    return {"factor": factor,
            "points": [{"true_flow": s["true_flow"], "TF_C": s["tp"]["TF_C"], "R_kg_h": float(s["r"]["R_kg_h"]),
                        "xD": float(s["r"]["xD"]), "xB": float(s["r"]["xB"])} for s in states]}
