"""A second reduction written apart from the reference: 2-minute block means, analyzer delays found by
cross-correlation, the air flow from an O2 balance block by block, everything else from the same physics."""
import csv
import json
import re
from datetime import datetime
from pathlib import Path

import warnings

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

AIR_O2 = 20.946e-2
MS, MSO2, MACID, MAIR = 32.06, 64.066, 98.079, 28.9647
BLOCK = 24  # samples of 5 s


def load(plant_dir):
    d = Path(plant_dir)
    meta = json.loads((d / "plant.json").read_text(encoding="utf-8"))
    with open(d / "trend.csv", encoding="utf-8") as fh:
        r = csv.reader(fh)
        head = next(r)
        data = [row for row in r]
    start = datetime.fromisoformat(data[0][0])
    sec = np.array([(datetime.fromisoformat(x[0]) - start).total_seconds() for x in data])
    arr = np.array([[float(v) if v else np.nan for v in x[1:]] for x in data])
    trend = {h: arr[:, i] for i, h in enumerate(head[1:])}
    events = []
    pat = re.compile(r"^(\S+) \| (.*)$")
    for line in (d / "log.txt").read_text(encoding="utf-8").splitlines():
        m = pat.match(line)
        if not m:
            continue
        when = (datetime.fromisoformat(m.group(1)) - start).total_seconds()
        fields = m.group(2).split(" | ")
        kv = dict(f.split("=", 1) for f in fields if "=" in f)
        events.append((when, [f for f in fields if "=" not in f], kv))
    cyl = {}
    with open(d / "cylinders.csv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            cyl[row["cylinder"]] = float(row["certified"])
    return meta, sec, trend, events, cyl


def calibrate(tag, x, sec, events, cyl):
    zt, zv, gt, gv = [], [], [], []
    for when, words, kv in events:
        if words[:2] != ["CAL", tag]:
            continue
        if kv["kind"] == "zero-adjust":
            zt.append(when); zv.append(0.0)
        elif kv["kind"] == "span-adjust":
            gt.append(when); gv.append(1.0)
        elif kv["kind"] == "zero-check":
            zt.append(when); zv.append(float(kv["as_found"]))
        elif kv["kind"] == "span-check":
            ref = 100 * AIR_O2 if kv["gas"] == "ambient-air" else cyl[kv["gas"]]
            z = np.interp(when, zt, zv) if zt else 0.0
            gt.append(when); gv.append((float(kv["as_found"]) - z) / ref)
    z = np.interp(sec, zt, zv) if zt else 0.0
    g = np.interp(sec, gt, gv) if gt else 1.0
    y = (x - z) / g
    # remove calibration-gas stretches
    bad = np.zeros(sec.size, bool)
    open_at = None
    for when, words, kv in events:
        if words[0] == "CAL" and words[1] == tag:
            bad |= (sec > when - 90) & (sec < when + 270)
        if words[0] == "CAL-PERIOD" and tag in words[1]:
            if words[2] == "begin":
                open_at = when
            else:
                bad |= (sec > open_at - 90) & (sec < when + 200)
    y[bad] = np.nan
    return y


def delay_of(y, ref, max_lag=40):
    """Lag (samples) at which y best follows ref."""
    a = np.nan_to_num(y - np.nanmean(y))
    b = np.nan_to_num(ref - np.nanmean(ref))
    best, lag = -np.inf, 0
    for k in range(0, max_lag):
        c = np.dot(a[k:], b[:len(b) - k])
        if c > best:
            best, lag = c, k
    return lag


def shift(y, k):
    return np.concatenate([y[k:], np.full(k, np.nan)]) if k else y


def blocks(x, mask):
    n = x.size // BLOCK
    xx = x[:n * BLOCK].reshape(n, BLOCK)
    mm = mask[:n * BLOCK].reshape(n, BLOCK)
    with np.errstate(invalid="ignore"):
        out = np.where(mm.all(axis=1), np.nanmean(np.where(mm, xx, np.nan), axis=1), np.nan)
    return out


def solve(plant_dir):
    meta, sec, tr, ev, cyl = load(plant_dir)
    inst = {i["tag"]: i for i in meta["instruments"]}
    t_start = [w for w, words, _ in ev if words == ["TEST", "start"]][0]
    t_end = [w for w, words, _ in ev if words == ["TEST", "end"]][0]
    win = (sec >= t_start) & (sec < t_end)
    sulfur = tr["FT-101"] * meta["sulfur"]["S_mass_share"] / MS
    # sulfur that burns: 40 s gas hold-up of furnace and boiler
    burn = sulfur.copy()
    for i in range(1, burn.size):
        burn[i] = burn[i - 1] + (5.0 / 45.0) * (sulfur[i] - burn[i - 1])
    rho = tr["PT-201"] * MAIR / (8.314462 * (tr["TT-201"] + 273.15))
    air_v = meta["air_venturi"]["C_kg_h"] * np.sqrt(np.clip(tr["PDT-201"], 0, None) * rho) / MAIR
    o2 = calibrate("AT-302", tr["AT-302"], sec, ev, cyl) / 100
    o2 = shift(o2, delay_of(-o2, burn))
    trip_spans = []
    st = [w for w, words, _ in ev if "sulfur-feed-stopped" in words]
    rs = [w for w, words, _ in ev if "sulfur-feed-restarted" in words]
    trip_spans = list(zip(st, rs))
    quiet = np.ones(sec.size, bool)
    for a, b in trip_spans:
        quiet &= ~((sec > a - 200) & (sec < b + 1500))
    for w, words, _ in ev:
        if "setpoint-change" in words:
            quiet &= ~((sec > w) & (sec < w + 1800))
    use = quiet & win & np.isfinite(o2)
    # air from the O2 balance with the furnace SO3 taken as 1 % of the sulfur (it moves the result by 0.2 %)
    f = 0.01
    air_o2 = sulfur * (1 + 0.5 * f - 1.5 * f * o2) / (AIR_O2 - o2)
    ratio = np.nanmedian(blocks(air_o2 / air_v, use))
    air = ratio * air_v

    def basis_factor(tag, x_share, absorbed):
        """moles of gas the analyzer refers its reading to, per mole of air flow terms"""
        s = inst[tag].get("sampling", "")
        so3_out = burn * (x_share - absorbed)
        n_tot = air - 0.5 * burn * x_share - burn * absorbed
        return n_tot - so3_out if "acid scrubber" in s else n_tot

    ipa = calibrate("AT-401", tr["AT-401"], sec, ev, cyl) / 100
    ipa = shift(ipa, delay_of(ipa, -burn / air, 60) if np.nanstd(ipa) > 0 else 0)
    x = 0.95
    for _ in range(8):
        q = ipa * basis_factor("AT-401", x, 0.0)
        m = win & np.isfinite(q)
        x = 1 - np.nansum(q[m]) / np.nansum(burn[m])
    x_ipa = x
    ppm = calibrate("AT-501", tr["AT-501"], sec, ev, cyl) * 1e-6
    so2_st_o2 = calibrate("AT-502", tr["AT-502"], sec, ev, cyl) / 100
    for a, b in trip_spans:
        m = (sec > a + 240) & (sec < b + 60) & (so2_st_o2 > 0.205)
        if m.sum() > 8:
            ppm = ppm - np.nanmedian(ppm[m])
    xt = 0.998
    for _ in range(8):
        n_st = air - 1.5 * burn * xt
        o2_expected = (AIR_O2 * air - burn * (1 + 0.5 * xt)) / n_st
        m = use & np.isfinite(so2_st_o2)
        dil = (np.nanmedian(blocks(so2_st_o2, m)) - np.nanmedian(blocks(o2_expected, m))) / \
              (AIR_O2 - np.nanmedian(blocks(o2_expected, m)))
        dil = max(dil, 0.0) if dil > 0.012 else 0.0
        q = ppm / (1 - dil) * n_st
        # fill the calibration gaps by interpolation of the molar flow
        good = np.isfinite(q)
        q = np.interp(np.arange(q.size), np.flatnonzero(good), q[good])
        xt = 1 - q[win].sum() / burn[win].sum()
    kg_t = 1000 * MSO2 * q[win].sum() / (MACID * (burn * xt)[win].sum())
    return {"ipa_conversion_pct": 100 * x_ipa, "so2_kg_per_t": float(kg_t)}


reduce = solve
