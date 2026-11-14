"""Wrong routes through the closing records, each one step away from the reference, scored like the judge (pick-up 2 %, release 0.60 ms).

python3 wrong_routes.py <draws folder> [draw ...]      (draw folders made by ../model/build_data.py <dir> <draw>)

Per closing the reference's own steps are reused (time base, coil fit, pole touches); a route changes one decision.
Writes wrong_routes.csv here: route, draw, unit, error of each key in ms, inside both limits or not.
"""
import csv
import json
import math
import os
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "solution"))
import close_timing as R  # noqa: E402

ALPHA = 0.00393
LIMIT_PCT, LIMIT_MS = 2.0, 0.60


def per_closing(folder, rescale=True):
    folder, sheet, rows = R.load(folder)
    net = sheet["network"]
    l_loop = net["generator_subtransient_inductance_h"] + net["transformer_leakage_inductance_h_lv"]
    obs = []
    for row in rows:
        rec = R.read_dfr(folder / "dfr" / (row["closing"] + ".csv"))
        t = rec["t_ms"] / 1000.0
        if rescale:
            pre = t < -0.005
            k = R.sine_frequency(t[pre], rec["VB_A"][pre]) / float(row["sync_bus_hz"])
            if abs(k - 1) > 0.003:
                t = t * k
        t_e, i_inf, tau = R.coil_fit(t, rec["ICOIL"])
        t_e_pure = R.coil_fit(t, rec["ICOIL"], allow_lag=False)[0]
        sigma = np.std(rec["IG_A"][t < 0])
        sm = np.convolve(rec["ICOIL"], np.ones(3) / 3, mode="same")
        on = int(np.flatnonzero(rec["ICOIL"] > 0.15)[0])
        dip = int(np.flatnonzero((np.diff(sm) < -0.02) & (np.arange(sm.size - 1) > on + 20))[0])
        onset = {p: float(t[np.flatnonzero(np.abs(rec["IG_" + p]) > 10 * sigma)[0]]) for p in "ABC"}
        obs.append(dict(t_dip=float(t[dip]), i_dip=float(sm[dip]), t_e=t_e, t_e_pure=t_e_pure, i_inf=i_inf, tau=tau, touch=R.pole_touches(rec, t, l_loop), onset=onset,
                        temp=float(row["cubicle_temp_c"]), vbat=float(row["battery_v"]), vdc=float(np.mean(rec["VDC"]))))
    return sheet, obs


def safe_log(x):
    return math.log(max(x, 1e-9))


def law_fit(obs, latch_of, start_of=lambda o: o["t_e"]):
    """I_pu and the three pole offsets from touch - start = latch(I_pu, closing) + c_p over all closings."""
    def resid(x):
        out = []
        for o in obs:
            lt = latch_of(x[0], o)
            out += [lt + c - (o["touch"][p] - start_of(o)) for p, c in zip("ABC", x[1:])]
        return np.array(out)
    return least_squares(resid, [5.0, 0.035, 0.035, 0.035], bounds=([0.5, -0.1, -0.1, -0.1], [20, 0.2, 0.2, 0.2]),
                         loss="soft_l1", f_scale=1e-4).x


def keys(i_pu, c):
    """The two graded values: pick-up current (A) and mean release-to-touch time (ms)."""
    return float(i_pu), 1000 * float(np.mean(c))


def loop_from_tau(sheet, obs, source="vdc", start_of=lambda o: o["t_e"]):
    """The reference's decisions: battery voltage, loop resistance from the coil time constant."""
    r20, l = sheet["closing_coil"]["resistance_ohm_at_20c"], sheet["closing_coil"]["inductance_h"]
    rest = np.mean([l / o["tau"] - r20 * (1 + ALPHA * (o["temp"] - 20)) for o in obs])

    def latch(ip, o):
        r = r20 * (1 + ALPHA * (o["temp"] - 20)) + rest
        return -(l / r) * safe_log(1 - ip * r / o[source])
    x = law_fit(obs, latch, start_of)
    return keys(x[0], x[1:])


def clamp_current(sheet, obs):
    """Latch at a fixed current read straight off the clamp: t = -tau ln(1 - I_pu / I_inf)."""
    x = law_fit(obs, lambda ip, o: -o["tau"] * safe_log(1 - ip / o["i_inf"]))
    return keys(x[0], x[1:])


def documented_loop(sheet, obs, source="vbat"):
    """Battery voltage over the documented resistances: coil copper at the cubicle temperature plus the cable."""
    r20, l = sheet["closing_coil"]["resistance_ohm_at_20c"], sheet["closing_coil"]["inductance_h"]
    rc = sheet["control_cable_loop_resistance_ohm"]

    def latch(ip, o):
        r = r20 * (1 + ALPHA * (o["temp"] - 20)) + rc
        return -(l / r) * safe_log(1 - ip * r / o[source])
    x = law_fit(obs, latch)
    return keys(x[0], x[1:])


def certificate_pickup(sheet, obs):
    """Pick-up current from the factory certificate's voltage dependence; the field closings give only the offsets."""
    r20, l = sheet["closing_coil"]["resistance_ohm_at_20c"], sheet["closing_coil"]["inductance_h"]
    un = sheet["rated_control_voltage_v"]
    pts = [(row["coil_voltage_fraction"] * un, [row["pole_%s_ms" % p] / 1000 for p in "ABC"])
           for row in sheet["factory_timing_certificate"]["closing_ms"]]

    def resid(x):
        return np.array([-(l / r20) * safe_log(1 - x[0] * r20 / v) + x[1 + k] - ts[k] for v, ts in pts for k in range(3)])
    ip = least_squares(resid, [5.0, 0.03, 0.03, 0.03]).x[0]
    c = [np.median([o["touch"][p] - o["t_e"] + o["tau"] * safe_log(1 - ip / o["i_inf"]) for o in obs]) for p in "ABC"]
    return keys(ip, c)


def armature_dip(sheet, obs):
    """Latch release read off the clamp as the start of the armature dip: the clamp current there is the pick-up
    current, and each pole's release-to-touch time is its touch less that instant."""
    return keys(np.mean([o["i_dip"] for o in obs]),
                [np.mean([o["touch"][p] - o["t_dip"] for o in obs]) for p in "ABC"])


def low_effort(sheet, obs):
    """A quick first program: the armature dip as release, and each pole's touch where its current first stands out
    of the noise."""
    return keys(np.mean([o["i_dip"] for o in obs]),
                [np.mean([o["onset"][p] - o["t_dip"] for o in obs]) for p in "ABC"])


ROUTES = [
    ("reference_decisions", "battery bus voltage, loop resistance from tau, field pick-up current",
     lambda s, o, u: loop_from_tau(s, o)),
    ("command_as_energisation", "release referred to BI_CLOSE instead of the start of coil current",
     lambda s, o, u: loop_from_tau(s, o, start_of=lambda q: 0.0)),
    ("time_stamps_trusted", "recorder time stamps taken as they are",
     lambda s, o, u: loop_from_tau(s, u)),
    ("clamp_without_lag", "coil current fitted as a clean exponential, the clamp's own lag not allowed for",
     lambda s, o, u: loop_from_tau(s, o, start_of=lambda q: q["t_e_pure"])),
    ("clamp_current_latch", "latch at a fixed clamp current, I_inf from the clamp",
     lambda s, o, u: clamp_current(s, o)),
    ("documented_loop_battery", "battery_v over coil copper plus documented cable resistance",
     lambda s, o, u: documented_loop(s, o, "vbat")),
    ("documented_loop_vdc", "VDC over coil copper plus documented cable resistance",
     lambda s, o, u: documented_loop(s, o, "vdc")),
    ("certificate_pickup", "pick-up current from the factory certificate, offsets from the field",
     lambda s, o, u: certificate_pickup(s, o)),
    ("armature_dip_release", "latch release taken at the start of the armature dip on the clamp, pick-up current read there",
     lambda s, o, u: armature_dip(s, o)),
    ("low_effort_replica", "armature dip as release, current onsets as touches",
     lambda s, o, u: low_effort(s, o)),
]


def detect(job):
    """Per-closing results of one unit of one draw, with and without the time-base check (worker process)."""
    folder = job
    s, o = per_closing(folder)
    _, u = per_closing(folder, rescale=False)
    return s, o, u


def main(top, draws):
    rows = []
    jobs = [(d, uid) for d in draws for uid in json.load(open(os.path.join(top, "draw%d" % d, "truth.json")))]
    todo = [(d, uid) for d, uid in jobs if not os.path.exists(os.path.join(top, "draw%d" % d, "routes_cache.pkl"))]
    with ProcessPoolExecutor() as pool:
        found = dict(zip(todo, pool.map(detect, [os.path.join(top, "draw%d" % d, uid) for d, uid in todo])))
    for draw in draws:
        gen = os.path.join(top, "draw%d" % draw)
        truth = json.load(open(os.path.join(gen, "truth.json")))
        cache = os.path.join(gen, "routes_cache.pkl")
        if os.path.exists(cache):
            det = pickle.load(open(cache, "rb"))
        else:
            det = {uid: found[(draw, uid)] for uid in truth}
            pickle.dump(det, open(cache, "wb"))
        for name, what, f in ROUTES:
            line = "draw %d %-24s" % (draw, name)
            for uid in truth:
                s, o, u = det[uid]
                ip, rt = f(s, o, u)
                ea = 100 * (ip / truth[uid]["pickup_current_a"] - 1)
                er = rt - truth[uid]["release_to_touch_ms"]
                inside = abs(ea) <= LIMIT_PCT and abs(er) <= LIMIT_MS
                rows.append([name, what, draw, uid, "%.3f" % ea, "%.3f" % er, int(inside)])
                line += " %s%s%+.1f%%/%+.2f" % (uid[-2:], "." if inside else "X", ea, er)
            print(line, flush=True)
    with open(os.path.join(HERE, "wrong_routes.csv"), "w", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["route", "decision", "draw", "unit", "pickup_error_pct", "release_error_ms", "inside_limits"])
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1], [int(d) for d in sys.argv[2:]] or [0, 1, 2, 3, 4])
