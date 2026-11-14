"""Temperature rises of a power transformer from one short-circuit temperature-rise test record.

python3 heat_run_rises.py <test record folder> <answer.json>

1. Readings: certificate corrections; 1-min log; switch-off instant and test current from the record.
2. Ambient: the oil cups, leaving out any cup that warms relative to the coolest cup once the losses flow
   (a cup in the cooler exhaust), checked against the cup positions and the weather station.
3. Top oil: the two cover pockets; if they disagree, the one that is consistent with the cooler-inlet check
   thermocouple and the station indicator is kept (top oil cannot read below the cooler inlet oil).
4. End of the total-loss period (last 60 min): top-oil rise and average-oil rise (top and bottom mean),
   corrected to rated losses with (P_N/P)^x.
5. Windings: resistance after switch-off -> temperature with the conductor's constant (235 copper, 225
   aluminium) and the cold resistance at the cold average oil temperature; while the measuring current
   is still settling the reading carries an inductive excess proportional to the current deficit, fitted as
   a factor on the resistance; the winding-to-average-oil difference is fitted as b exp(-t/tau) against the
   logged average oil (which keeps cooling) and extrapolated to the switch-off second. Gradient corrected to
   rated current with (I_N/I)^y; average winding rise = corrected average-oil rise + corrected gradient.
"""
import csv
import json
import re
import sys

import numpy as np
from scipy.optimize import least_squares

EXPONENTS = {"ONAN": (0.8, 1.6), "ONAF": (0.9, 1.6), "ODAF": (1.0, 2.0)}
KMAT = {"copper": 235.0, "aluminium": 225.0}


def clock(s):
    parts = [int(x) for x in s.strip().split(":")]
    return parts[0] * 60 + parts[1] + (parts[2] / 60.0 if len(parts) > 2 else 0.0)


def load(folder):
    sheet = json.load(open(folder + "/transformer.json"))
    with open(folder + "/heat_run.csv") as fh:
        rd = list(csv.DictReader(fh))
    log = {k: np.array([float(r[k]) for r in rd]) for k in rd[0] if k != "time"}
    t = np.array([clock(r["time"]) for r in rd]) + 0.5  # 1-min averages centred on the half minute
    with open(folder + "/winding_resistance.csv") as fh:
        res = list(csv.DictReader(fh))
    rec = open(folder + "/test_record.txt").read()
    return sheet, t, log, res, rec


def analyse(folder):
    sheet, t, log, res, rec = load(folder)
    sens = {s["id"]: s for s in sheet["sensors"]}

    def corr(tag):
        return log[tag] + sens[tag].get("certificate_correction_K", 0.0)

    x, y = EXPONENTS[sheet["cooling"]]
    m = re.search(r"^(\d\d:\d\d)\s+supply on", rec, re.M)
    t_on = clock(m.group(1))
    e1 = clock(re.search(r"^(\d\d:\d\d)\s+end of total-loss", rec, re.M).group(1))
    t_off = clock(re.search(r"^(\d\d:\d\d:\d\d)\s+supply off", rec, re.M).group(1))
    m = re.search(r"top oil \(TOT-1\) reads ([\d.]+) C, bottom oil \(BOT\) reads ([\d.]+) C", rec)
    theta1 = 0.5 * (float(m.group(1)) + sens["TOT-1"]["certificate_correction_K"]
                    + float(m.group(2)) + sens["BOT"]["certificate_correction_K"])
    m = re.search(r"HV \(line-to-line, \S+\) ([\d.]+) mOhm, LV \(line-to-line, \S+\) ([\d.]+) mOhm", rec)
    R1 = {"HV": float(m.group(1)), "LV": float(m.group(2))}

    end1 = (t > e1 - 60) & (t < e1)
    pre = t < t_on
    # ambient cups: exclude those that warm relative to the coolest cup once losses flow
    cups = [s for s in sheet["sensors"] if s["measures"] == "ambient air" and s["grade"] == "test"]
    rise = {}
    for s in cups:
        rise[s["id"]] = None
    base = {s["id"]: np.mean(corr(s["id"])[pre]) for s in cups}
    endv = {s["id"]: np.mean(corr(s["id"])[end1]) for s in cups}
    coolest = min(cups, key=lambda s: endv[s["id"]] - base[s["id"]])["id"]
    keep = [s["id"] for s in cups
            if (endv[s["id"]] - base[s["id"]]) - (endv[coolest] - base[coolest]) < 1.0]
    amb = np.mean([corr(k) for k in keep], axis=0)
    # top oil pockets
    p1, p2 = corr("TOT-1"), corr("TOT-2")
    d = abs(np.mean(p1[end1]) - np.mean(p2[end1]))
    if d > 0.5:
        ref = np.mean(log["OTI"][end1])
        tci = np.mean(log["TCI"][end1])
        cand = [(abs(np.mean(p[end1]) - ref) + (5.0 if np.mean(p[end1]) < tci else 0.0), i)
                for i, p in enumerate((p1, p2))]
        top = (p1, p2)[min(cand)[1]]
    else:
        top = 0.5 * (p1 + p2)
    bot = corr("BOT")
    avg_oil = 0.5 * (top + bot)
    P = np.mean(log["P_kW"][end1])
    PN = sheet["guaranteed_losses_kW"]["total"]
    fP = (PN / P) ** x
    top_rise = (np.mean(top[end1]) - np.mean(amb[end1])) * fP
    om_rise = (np.mean(avg_oil[end1]) - np.mean(amb[end1])) * fP
    # test current in the gradient hour
    last = (t > t_off - 30) & (t < t_off - 1)
    i_pu = np.mean(log["I_A"][last]) / sheet["rated_current_HV_A"]
    out = {"top_oil_rise_K": float(top_rise)}
    diag = {"ambient_cups": keep, "P_test": float(P), "i_pu": float(i_pu)}
    # oil after switch-off (minute log, interpolated)
    for w in ("HV", "LV"):
        rows = [r for r in res if r["winding"] == w]
        ts = np.array([clock(r["time"]) for r in rows])
        I = np.array([float(r["current_A"]) for r in rows])
        R = np.array([float(r["resistance_mOhm"]) for r in rows])
        I_final = np.median(I[-10:])
        dev = 1.0 - I / I_final
        sig = max(np.std(dev[-20:]), 1e-6)
        smooth = np.array([np.median(dev[max(0, i - 2):i + 3]) for i in range(len(dev))])
        high = np.where(smooth > 2e-3)[0]
        first = int(high.max()) + 3 if len(high) else 0  # a short margin after the last unsettled reading
        ok = np.zeros(len(dev), bool)
        ok[first:] = True
        Kc = KMAT[sheet["windings"][w]["conductor"]]
        thw = R / R1[w] * (Kc + theta1) - Kc
        tt = (ts - t_off)
        oil = np.interp(ts, t, avg_oil)
        dd = thw - oil

        # the reading still carries the inductive term while the current settles; it scales with the
        # current deficit, so it is fitted as a factor (1 + k * deficit) on the resistance
        def resid(v):
            th_c = R[ok] / (1.0 + v[2] * dev[ok]) / R1[w] * (Kc + theta1) - Kc
            return v[0] * np.exp(-tt[ok] / v[1]) - (th_c - oil[ok])
        f = least_squares(resid, [15.0, 7.0, 0.0], x_scale=[1.0, 1.0, 1.0])
        b, tau, kdc = f.x
        oil_off = float(np.interp(t_off, t, avg_oil))
        g = b  # winding minus average oil at switch-off
        gN = g * (1.0 / i_pu) ** y
        out["%s_winding_rise_K" % w.lower()] = float(om_rise + gN)
        diag[w] = dict(theta_w_off=oil_off + b, g=g, tau_min=tau, used=int(ok.sum()), of=len(ok),
                       rms_K=float(np.sqrt(np.mean(f.fun ** 2))))
    return {k: round(v, 4) for k, v in out.items()}, diag


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: heat_run_rises.py <test record folder> <answer file>")
    rises, notes = analyse(sys.argv[1])
    with open(sys.argv[2], "w", encoding="utf-8") as answer:
        answer.write(json.dumps(rises) + "\n")
    sys.stdout.write(json.dumps({"rises": rises, "notes": notes}, default=str) + "\n")
