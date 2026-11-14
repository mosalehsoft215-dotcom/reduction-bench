"""Static-source position error of one flight-test campaign, and the total-temperature probe's recovery factor.

    python3 static_source.py <campaign folder> <answer json>

Steps:
1. Records. Transducer volts become pressures through their calibrations; the static transducer's reading is
   corrected for its zero shift with cabin temperature. The fuel totalizer is scaled back by its stated error to
   give the fuel used, and the weight follows from the zero-fuel weight and the fuel at engine start.
2. Pressure lag. Each line acts like a first-order lag. Its time constant scales from the ground test as the
   inverse of the line's pressure (Poiseuille flow: tau ~ mu L V / (d^4 p)); the ground test was made on the
   runway (pressure from the tower barometer) with the cabin at its flight temperature, so the viscosity does not
   change. The true pressure is the reading plus tau times its rate of change.
3. Recovery factor. At one pressure altitude the static air temperature is the same for every run, so
   TAT = T + r V^2 kv / (2 cp) with V the true airspeed (TAT = T (1 + 0.2 r M^2) and M = V / a, the speed of
   sound of the moist air a^2 = gamma R T / kv, kv = 1 - 0.378 e / p, e from the recorded dew point); r is fitted over the four-leg runs (one T per pressure
   altitude), or over the level accelerations when a campaign has no four-leg run.
4. Free-stream static pressure:
   - fly-bys: the photo grid is graduated for the vertical plane of the runway centreline, so the static-port
     mark's height above the lens is the grid reading times d / D, d the camera's actual horizontal distance
     from the aircraft (square to the runway, from the GPS position at the passage) and D the distance to the
     centreline; the tower barometer carries the pressure up to that height (hypsometric, virtual temperature);
   - four-leg runs: true airspeed and wind from the four legs' GPS velocities (constant wind, one true airspeed);
   - level accelerations: the wind comes from the four-leg runs at that pressure altitude if there are any, else
     from the pair of runs on reciprocal headings; then V = |ground velocity - wind|.
   For GPS points, T = TAT - r V^2 / (2 cp), M = V / a(T) and p = total pressure / (1 + 0.2 M^2)^3.5.
5. Error model. (ps_indicated - p) / (pt - p) = a_flap + b_flap CL + c M^2, by least squares, with CL from the
   weight and the true dynamic pressure.
6. The requested corrections follow by iteration at the indicated condition: dVpc = Vc - Vi and
   dHpc = Hc - Hi, Hc and Hi both ISA pressure altitudes.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

R_AIR, G0, GAMMA = 287.053, 9.80665, 1.4
CP = GAMMA * R_AIR / (GAMMA - 1.0)
P_SL, T_SL, LAPSE = 101325.0, 288.15, 0.0065
A_SL = np.sqrt(GAMMA * R_AIR * T_SL)
KNOT, FOOT = 0.514444, 0.3048

OPTIONS = dict(lag=True, lag_pressure=True, lateral=True, zero_shift=True, totalizer=True, lift=True,
               mach_term=True, humid=True, antenna=True, humid_sound=True, stratosphere=True, pitch_arm=True)
H_TROP = 11000.0
T_TROP = T_SL - LAPSE * H_TROP
P_TROP = P_SL * (T_TROP / T_SL) ** (G0 / (R_AIR * LAPSE))
STRATOSPHERE = [True]          # the ISA relations above 11 km (switched off only by the misreading study)


def hp_of(p):
    """ISA pressure altitude (m): the troposphere, and the isothermal layer above 11 km"""
    if STRATOSPHERE[0] and p < P_TROP:
        return H_TROP - R_AIR * T_TROP / G0 * np.log(p / P_TROP)
    return T_SL / LAPSE * (1.0 - (p / P_SL) ** (R_AIR * LAPSE / G0))


def p_of(h):
    if STRATOSPHERE[0] and h > H_TROP:
        return P_TROP * np.exp(-G0 * (h - H_TROP) / (R_AIR * T_TROP))
    return P_SL * (1.0 - LAPSE * h / T_SL) ** (G0 / (R_AIR * LAPSE))


def cas(qc):
    return A_SL * np.sqrt(5.0 * ((qc / P_SL + 1.0) ** (2.0 / 7.0) - 1.0))


def qc_of(v):
    return P_SL * ((1.0 + 0.2 * (v / A_SL) ** 2) ** 3.5 - 1.0)


def mach(pt, ps):
    return np.sqrt(5.0 * ((pt / ps) ** (2.0 / 7.0) - 1.0))


def poly(c, v):
    return c[0] + c[1] * v + c[2] * v * v


# ---------------------------------------------------------------------------------------------- records
def load(folder, o):
    folder = Path(folder)
    ac = json.loads((folder / "aircraft.json").read_text())
    rows = list(csv.DictReader(open(folder / "flight.csv")))
    col = {k: np.array([float(r[k]) for r in rows]) for k in rows[0] if k != "card"}
    col["card"] = np.array([r["card"] for r in rows])
    tr = ac["transducers"]
    col["pt"] = poly(tr["total"]["kpa_from_volts"], col["total_v"]) * 1e3
    ps = poly(tr["static"]["kpa_from_volts"], col["static_v"])
    if o["zero_shift"]:
        ps = ps - tr["static"]["zero_shift_kpa_per_k"] * (col["cabin_c"] - tr["static"]["calibrated_at_c"])
    col["ps"] = ps * 1e3
    # moist air: speed of sound a^2 = gamma R T_v, T_v = T / kv, kv = 1 - 0.378 e / p (e from the dew point)
    e = 610.94 * np.exp(17.625 * col["dewpoint_c"] / (col["dewpoint_c"] + 243.04))
    col["kv"] = 1.0 - 0.378 * e / col["ps"] if o["humid_sound"] else np.ones_like(col["ps"])
    fuel = col["fuel_used_kg"] / (1.0 + ac["fuel_totalizer"]["reads_high_percent"] / 100.0 if o["totalizer"] else 1.0)
    col["weight"] = ac["weights"]["zero_fuel_kg"] + ac["weights"]["fuel_at_engine_start_kg"] - fuel
    cards = list(csv.DictReader(open(folder / "cards.csv")))
    tower = list(csv.DictReader(open(folder / "tower.csv"))) if (folder / "tower.csv").exists() else []
    if o["lag"]:
        unlag(col, ac, ground_pressure(ac, tower), o)
    return ac, col, cards, tower


def ground_pressure(ac, tower):
    """runway-level pressure (2 m up, where the lines were on the ground) from the first tower reading"""
    site = ac["site"]
    if not tower:
        return p_of(site["runway_elevation_m"] + 2.0)
    t = float(tower[0]["temp_c"]) + 273.15
    dz = site["runway_elevation_m"] + 2.0 - site["tower_barometer_height_m"]
    return float(tower[0]["baro_hpa"]) * 100.0 * np.exp(-G0 * dz / (R_AIR * t))


def unlag(col, ac, p_ground, o):
    """x = y + tau(y) dy/dt within each recorded block (a gap in time starts a new block)"""
    t = col["time_s"]
    edges = np.concatenate([[0], np.nonzero(np.diff(t) > 1.0)[0] + 1, [len(t)]])
    lines = ac["pitot_static"]["lines"]
    for key, line in (("pt", "total"), ("ps", "static")):
        y = col[key].copy()
        for a, b in zip(edges[:-1], edges[1:]):
            yy, tt = y[a:b], t[a:b]
            rate = smooth_rate(yy, tt)
            tau = lines[line]["ground_lag_s"] * np.ones_like(yy)
            if o["lag_pressure"]:
                tau = tau * p_ground / yy
            col[key][a:b] = yy + tau * rate


def smooth_rate(y, t, half=6):
    """slope of a local straight line over +-half samples (3 s each side at 2 Hz)"""
    out = np.zeros_like(y)
    n = len(y)
    for i in range(n):
        a, b = max(0, i - half), min(n, i + half + 1)
        if b - a >= 3:
            out[i] = np.polyfit(t[a:b] - t[i], y[a:b], 1)[0]
    return out


# ---------------------------------------------------------------------------------------------- methods
def window(col, label):
    return col["card"] == label


def flyby_points(ac, col, tower, o):
    site = ac["site"]
    hdg = np.radians(site["runway_heading_deg"])
    along, right = np.array([np.sin(hdg), np.cos(hdg)]), np.array([np.cos(hdg), -np.sin(hdg)])
    cam = site["camera"]
    pts = []
    for tw in tower:
        t0 = float(tw["time_s"])
        near = np.abs(col["time_s"] - t0) <= 3.0
        i = int(np.argmin(np.abs(col["time_s"] - t0)))
        if tw["grid_m"].strip():
            lateral = col["east_m"][i] * right[0] + col["north_m"][i] * right[1]
            d = cam["left_of_centreline_m"] + lateral if o["lateral"] else cam["left_of_centreline_m"]
            z = cam["lens_height_m"] + float(tw["grid_m"]) * d / cam["left_of_centreline_m"]
        else:
            # camera failed: the ports sit below the GPS antenna by its body-axis offset turned through the pitch
            g = ac["gps"]
            th = np.radians(col["pitch_deg"][near].mean()) if o["pitch_arm"] else 0.0
            dz = g["antenna_above_static_ports_m"] * np.cos(th) - g.get("antenna_aft_of_static_ports_m", 0.0) * np.sin(th)
            z = col["gps_height_m"][near].mean() - (dz if o["antenna"] else 0.0)
        tk = float(tw["temp_c"]) + 273.15
        if o["humid"]:
            e = float(tw["rh_pct"]) / 100.0 * 610.94 * np.exp(17.625 * (tk - 273.15) / (tk - 273.15 + 243.04))
            tk = tk / (1.0 - 0.378 * e / (float(tw["baro_hpa"]) * 100.0))
        p = float(tw["baro_hpa"]) * 100.0 * np.exp(-G0 * (z - site["tower_barometer_height_m"]) / (R_AIR * tk))
        pts.append(point(col, near, p, flap=col["flap_deg"][i], kind="flyby"))
    return pts


def point(col, sel, p, flap, kind, m=None):
    pt, ps = col["pt"][sel].mean(), col["ps"][sel].mean()
    m = mach(pt, p) if m is None else m
    w = col["weight"][sel].mean()
    return dict(kind=kind, flap=int(flap), p=p, pt=pt, ps=ps, m=m, w=w, f=(ps - p) / (pt - p))


def leg_mean(col, sel):
    return dict(ve=col["v_east"][sel].mean(), vn=col["v_north"][sel].mean(), tat=col["tat_c"][sel].mean() + 273.15,
                hdg=np.radians(np.degrees(np.arctan2(np.sin(np.radians(col["heading_deg"][sel])).mean(),
                                                      np.cos(np.radians(col["heading_deg"][sel])).mean()))))


def fourleg_solve(legs):
    """TAS and wind from ground velocities of legs flown at one airspeed: |Vg - W| = V for each leg"""
    g = np.array([[lg["ve"], lg["vn"]] for lg in legs])
    w = g.mean(0)
    v = np.linalg.norm(g - w, axis=1).mean()
    for _ in range(50):
        d = g - w
        r = np.linalg.norm(d, axis=1)
        j = np.column_stack([-d[:, 0] / r, -d[:, 1] / r, -np.ones(len(g))])
        step = np.linalg.lstsq(j, -(r - v), rcond=None)[0]
        w = w + step[:2]
        v = v + step[2]
        if np.abs(step).max() < 1e-9:
            break
    return v, w


def gps_runs(ac, col, cards, o):
    """four-leg runs and level accelerations -> lists of temperature samples and stable points (r pending)"""
    groups = {}
    for cd in cards:
        if cd["card"].startswith("L"):
            groups.setdefault(cd["card"].split("-")[0], []).append(cd)
    fourlegs = []
    for key, cds in groups.items():
        legs = [leg_mean(col, window(col, cd["card"])) for cd in cds]
        sels = [window(col, cd["card"]) for cd in cds]
        v, w = fourleg_solve(legs)
        sel = np.any(sels, axis=0)
        fourlegs.append(dict(key=key, v=v, wind=w, tat=np.mean([lg["tat"] for lg in legs]), sel=sel,
                             kv=col["kv"][sel].mean(),
                             hp=float(cds[0]["pressure_altitude_ft"]), flap=int(cds[0]["flap_deg"]), legs=sels))
    accels = [cd for cd in cards if cd["card"].startswith("A")]
    runs = []
    for cd in accels:
        sel = window(col, cd["card"])
        hp = float(cd["pressure_altitude_ft"])
        runs.append(dict(card=cd["card"], sel=sel, hp=hp, flap=int(cd["flap_deg"])))
    for run in runs:
        same = [fl for fl in fourlegs if abs(fl["hp"] - run["hp"]) < 300]
        if same:
            run["wind"] = np.mean([fl["wind"] for fl in same], axis=0)
        else:
            run["wind"] = pair_wind(col, run, reciprocal(run, runs))
        g = np.column_stack([col["v_east"][run["sel"]], col["v_north"][run["sel"]]])
        run["v"] = np.linalg.norm(g - run["wind"], axis=1)
        run["tat"] = col["tat_c"][run["sel"]] + 273.15
        run["kv"] = col["kv"][run["sel"]]
    return fourlegs, runs


def reciprocal(run, runs):
    for other in runs:
        if other is not run and abs(other["hp"] - run["hp"]) < 300:
            return other
    return None


def pair_wind(col, a, b):
    """wind from two runs on reciprocal headings: at the same indicated airspeed the true airspeed is the same,
    so the mean of the two ground velocities is the wind"""
    def ias(sel):
        return cas(col["pt"][sel] - col["ps"][sel])
    va, vb = ias(a["sel"]), ias(b["sel"])
    ga = np.column_stack([col["v_east"][a["sel"]], col["v_north"][a["sel"]]])
    gb = np.column_stack([col["v_east"][b["sel"]], col["v_north"][b["sel"]]])
    lo, hi = max(va.min(), vb.min()), min(va.max(), vb.max())
    grid = np.linspace(lo, hi, 40)[2:-2]
    oa, ob = np.argsort(va), np.argsort(vb)
    wa = np.column_stack([np.interp(grid, va[oa], ga[oa, k]) for k in (0, 1)])
    wb = np.column_stack([np.interp(grid, vb[ob], gb[ob, k]) for k in (0, 1)])
    return (0.5 * (wa + wb)).mean(0)


def recovery(fourlegs, runs):
    """TAT = T_run + r V^2 / (2 cp): least squares, one static temperature per run or per pressure altitude"""
    rows, y, groups = [], [], []
    for fl in fourlegs:
        groups.append(("hp", round(fl["hp"])))
    for run in runs:
        groups.append(("hp", round(run["hp"])))
    names = sorted(set(groups), key=str)
    idx = {g: i for i, g in enumerate(names)}
    for fl in fourlegs:
        r = np.zeros(len(names) + 1)
        r[idx[("hp", round(fl["hp"]))]] = 1.0
        r[-1] = fl["v"] ** 2 * fl["kv"] / (2 * CP)
        rows.append(r)
        y.append(fl["tat"])
    for run in runs:
        for v, tt, kv in zip(run["v"][::4], run["tat"][::4], run["kv"][::4]):
            r = np.zeros(len(names) + 1)
            r[idx[("hp", round(run["hp"]))]] = 1.0
            r[-1] = v * v * kv / (2 * CP)
            rows.append(r)
            y.append(tt)
    sol = np.linalg.lstsq(np.array(rows), np.array(y), rcond=None)[0]
    return float(sol[-1])


def gps_points(col, fourlegs, runs, r):
    pts = []
    for fl in fourlegs:
        for sel in fl["legs"]:
            t = fl["tat"] - r * fl["v"] ** 2 * fl["kv"] / (2 * CP)
            m = fl["v"] / np.sqrt(GAMMA * R_AIR * t / fl["kv"])
            pt = col["pt"][sel].mean()
            p = pt / (1 + 0.2 * m * m) ** 3.5
            pts.append(point(col, sel, p, fl["flap"], "fourleg", m=m))
    for run in runs:
        idx = np.nonzero(run["sel"])[0]
        for k in range(0, len(idx) - 8, 8):
            seg = idx[k:k + 8]
            sel = np.zeros_like(run["sel"])
            sel[seg] = True
            v = run["v"][k:k + 8].mean()
            kv = run["kv"][k:k + 8].mean()
            t = run["tat"][k:k + 8].mean() - r * v * v * kv / (2 * CP)
            m = v / np.sqrt(GAMMA * R_AIR * t / kv)
            pt = col["pt"][sel].mean()
            p = pt / (1 + 0.2 * m * m) ** 3.5
            pts.append(point(col, sel, p, run["flap"], "accel", m=m))
    return pts


# ---------------------------------------------------------------------------------------------- model
def lift_coefficient(pt_or_none, p, m, w, wing):
    return w * G0 / (0.5 * GAMMA * p * m * m * wing)


def fit(points, ac, o):
    wing = ac["aircraft"]["wing_area_m2"]
    flaps = sorted({pp["flap"] for pp in points})
    rows, y = [], []
    for pp in points:
        x = lift_coefficient(None, pp["p"], pp["m"], pp["w"], wing) if o["lift"] else cas(pp["pt"] - pp["ps"]) / KNOT / 100.0
        r = np.zeros(2 * len(flaps) + 1)
        k = flaps.index(pp["flap"])
        r[2 * k], r[2 * k + 1] = 1.0, x
        r[-1] = pp["m"] ** 2 if o["mach_term"] else 0.0
        rows.append(r)
        y.append(pp["f"])
    sol = np.linalg.lstsq(np.array(rows), np.array(y), rcond=None)[0]
    return {fl: (sol[2 * k], sol[2 * k + 1], sol[-1]) for k, fl in enumerate(flaps)}


def correction(req, coef, ac, o):
    wing = ac["aircraft"]["wing_area_m2"]
    a, b, c = coef[int(req["flap_deg"])]
    ps_i = p_of(req["pressure_altitude_ft"] * FOOT)
    pt = ps_i + qc_of(req["ias_kt"] * KNOT)
    dp = 0.0
    for _ in range(200):
        p = ps_i - dp
        m = mach(pt, p)
        x = lift_coefficient(None, p, m, req["weight_kg"], wing) if o["lift"] else req["ias_kt"] / 100.0
        new = (a + b * x + c * m * m) * (pt - p)
        if abs(new - dp) < 1e-9:
            break
        dp = new
    p = ps_i - dp
    if req["key"] == "dvpc_kt":
        return float(cas(pt - p) / KNOT - req["ias_kt"])
    return float((hp_of(p) - req["pressure_altitude_ft"] * FOOT) / FOOT)


def solve(folder, **changes):
    o = dict(OPTIONS, **changes)
    STRATOSPHERE[0] = o["stratosphere"]
    ac, col, cards, tower = load(folder, o)
    fourlegs, runs = gps_runs(ac, col, cards, o)
    # four-leg runs average the gusts over whole legs: they alone set r when there are any
    r = recovery(fourlegs, []) if fourlegs else recovery([], runs)
    pts = flyby_points(ac, col, tower, o) + gps_points(col, fourlegs, runs, r)
    coef = fit(pts, ac, o)
    out = dict(recovery_factor=r)
    for req in ac["requested"]:
        out[req["key"]] = correction(req, coef, ac, o)
    return out, dict(coef={str(k): [float(x) for x in v] for k, v in coef.items()}, points=len(pts))


if __name__ == "__main__":
    result, info = solve(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result), json.dumps(info))
