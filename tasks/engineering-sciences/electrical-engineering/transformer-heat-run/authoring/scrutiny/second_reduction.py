"""Second, separately coded reduction of a temperature-rise test record.

How it departs from solution/heat_run_rises.py:
  * resistance recomputed from voltage and current; readings kept from two readings after the five-point mean
    current deficit first falls below 3e-5, with no model of the settling term;
  * the logged average oil after switch-off smoothed by its own fitted exponential a + c exp(-t/To); the
    winding fitted as that oil curve plus b exp(-t/T) (the reference interpolates the logged oil and models
    the settling term instead of dropping it);
  * ambient cups judged against the weather station: a cup whose offset to the station grows by more than
    2 K between the hour before the supply and the last hour of the total-loss period is left out;
  * top-oil pocket judged by physics: a pocket reading below the cooler-inlet check thermocouple plus 0.3 K
    is left out (the oil cools on its way to the cooler);
  * steady values from a straight-line fit over the last 60 minutes evaluated at the period end.
python second_reduction.py <test folder>   writes one record's rises to stdout as JSON
python second_reduction.py                 all seven records of the bundle -> second_reduction.csv
"""
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


def pjoin(*parts):
    return str(Path(*parts))

TABLE = {"ONAN": (0.8, 1.6), "ONAF": (0.9, 1.6), "ODAF": (1.0, 2.0)}
CONST = {"copper": 235.0, "aluminium": 225.0}


def minutes(s):
    h, m, *sec = s.split(":")
    return int(h) * 60 + int(m) + (int(sec[0]) / 60.0 if sec else 0.0)


def reduce(folder):
    info = json.load(open(folder + "/transformer.json"))
    cal = {s["id"]: s.get("certificate_correction_K", 0.0) for s in info["sensors"]}
    log = list(csv.reader(open(folder + "/heat_run.csv")))
    head = log[0]
    data = np.array([[float(v) for v in r[1:]] for r in log[1:]])
    tm = np.array([minutes(r[0]) + 0.5 for r in log[1:]])
    col = {h: data[:, i - 1] + cal.get(h, 0.0) for i, h in enumerate(head) if i > 0}
    text = open(folder + "/test_record.txt").read()
    on = minutes(re.search(r"(\d\d:\d\d)\s+supply on", text).group(1))
    endtl = minutes(re.search(r"(\d\d:\d\d)\s+end of total-loss", text).group(1))
    off = minutes(re.search(r"(\d\d:\d\d:\d\d)\s+supply off", text).group(1))
    vals = [float(v) for v in re.findall(r"reads ([\d.]+) C", text)]
    th_cold = 0.5 * (vals[0] + cal["TOT-1"] + vals[1] + cal["BOT"])
    r_cold = dict(zip(("HV", "LV"), [float(v) for v in re.findall(r"\) ([\d.]+) mOhm", text)]))
    x, y = TABLE[info["cooling"]]

    win = (tm > endtl - 60) & (tm < endtl)
    before = tm < on

    def at_end(v):
        c = np.polyfit(tm[win] - endtl, v[win], 1)
        return c[1]

    cups = [s["id"] for s in info["sensors"] if s["measures"] == "ambient air" and s["grade"] == "test"]
    good = []
    for cid in cups:
        grow = np.mean((col[cid] - col["WX"])[win]) - np.mean((col[cid] - col["WX"])[before])
        if grow < 2.0:
            good.append(cid)
    ambient = np.mean([col[c] for c in good], axis=0)
    pockets = [p for p in ("TOT-1", "TOT-2") if np.mean(col[p][win]) > np.mean(col["TCI"][win]) + 0.3]
    top = np.mean([col[p] for p in pockets], axis=0)
    mid = 0.5 * (top + col["BOT"])
    power = np.mean(col["P_kW"][win])
    scale = (info["guaranteed_losses_kW"]["total"] / power) ** x
    d_top = (at_end(top) - at_end(ambient)) * scale
    d_mid = (at_end(mid) - at_end(ambient)) * scale
    hour = (tm > off - 40) & (tm < off - 2)
    k_cur = np.mean(col["I_A"][hour]) / info["rated_current_HV_A"]
    rows = list(csv.DictReader(open(folder + "/winding_resistance.csv")))
    out = {"top_oil_rise_K": float(d_top)}
    for wd in ("HV", "LV"):
        rr = [r for r in rows if r["winding"] == wd]
        t = np.array([minutes(r["time"]) for r in rr]) - off
        cur = np.array([float(r["current_A"]) for r in rr])
        res = np.array([float(r["voltage_mV"]) for r in rr]) / cur
        deficit = 1.0 - cur / np.median(cur[-8:])
        run5 = np.convolve(deficit, np.ones(5) / 5, mode="same")
        below = np.where(run5 < 3e-5)[0]
        keep = np.arange(len(t)) >= (below.min() + 2 if len(below) else 0)
        kc = CONST[info["windings"][wd]["conductor"]]
        temp = res / r_cold[wd] * (kc + th_cold) - kc

        # the average oil after switch-off is smoothed by its own exponential, the winding then rides on it
        after = (tm > off) & (tm < off + 25)
        to, yo = tm[after] - off, mid[after]
        po, _ = curve_fit(lambda tt, a, c, T: a + c * np.exp(-tt / T), to, yo,
                          p0=[yo[-1] - 5.0, 5.0, 100.0], method="trf",
                          bounds=([-50.0, -100.0, 5.0], [200.0, 200.0, 2000.0]), maxfev=20000)
        oil_fit = po[0] + po[1] * np.exp(-t / po[2])
        tk, yk = t[keep], (temp - oil_fit)[keep]
        p, _ = curve_fit(lambda tt, b, T: b * np.exp(-tt / T), tk, yk, p0=[max(yk[0], 1.0), 7.0],
                         method="trf", bounds=([0.0, 1.0], [80.0, 30.0]), maxfev=20000)
        w_off = po[0] + po[1] + p[0]
        grad = (w_off - (po[0] + po[1])) * (1.0 / k_cur) ** y
        out["%s_winding_rise_K" % wd.lower()] = float(d_mid + grad)
    return out


def table():
    here = Path(__file__).resolve().parent
    bundle = here.parents[1]
    record = json.load(open(pjoin(bundle, "authoring", "thermal", "model_record.json")))
    keys = ("top_oil_rise_K", "hv_winding_rise_K", "lv_winding_rise_K")
    with open(pjoin(here, "second_reduction.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["record"] + list(keys) + ["miss_" + k for k in keys])
        for rid, v in record.items():
            r = reduce(pjoin(bundle, v["folder"]))
            w.writerow([rid] + ["%.3f" % r[k] for k in keys] + ["%+.3f" % (r[k] - v[k]) for k in keys])
            print(rid, ["%+.3f" % (r[k] - v[k]) for k in keys], flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(json.dumps(reduce(sys.argv[1])))
    else:
        table()
