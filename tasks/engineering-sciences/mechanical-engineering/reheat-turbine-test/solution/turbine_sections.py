"""HP and IP section efficiencies from one enthalpy-drop test campaign.

python3 hpip_test.py <campaign folder> <answer.json>

1. Window means of every tag (test windows from test_log.txt, both ends included).
2. Instrument corrections: certificates, post-test checks, barometric pressure for gauge transmitters,
   condensate legs at the turbine-hall air temperature.
3. Duplicate thermocouples are screened: a pair that separates during the campaign is resolved with
   the station crossover thermocouple, and the drifting element is dropped.
4. IP section: hot reheat steam and the mid-span packing leak (HP inlet steam, or control-stage shell
   steam on partial-arc units) mix before the IP blading. The leak ratio and the IP efficiency are
   fitted together to the IP exhaust enthalpies of all test points.
5. HP section: from the stop-valve state to the HP exhaust state. When the cold reheat thermowells sit
   downstream of the reheat desuperheater, the HP exhaust enthalpy is recovered from a mixing balance
   whose steam flow is the feedwater flow less the packing leak found in step 4.
"""
import json
import math
import re
import sys

import numpy as np
from iapws import IAPWS97
from scipy.optimize import least_squares

G = 9.80665


def h_pT(p_bar, T_C):
    return IAPWS97(P=p_bar / 10.0, T=T_C + 273.15).h


def s_ph(p_bar, h):
    return IAPWS97(P=p_bar / 10.0, h=h).s


def h_ps(p_bar, s):
    return IAPWS97(P=p_bar / 10.0, s=s).h


def rho_liq(p_bar, T_C):
    return IAPWS97(P=p_bar / 10.0, T=T_C + 273.15).rho


def read_trend(path):
    with open(path) as fh:
        head = fh.readline().strip().split(",")
        rows = [ln.strip().split(",") for ln in fh if ln.strip()]
    minutes = np.array([int(r[0][:2]) * 60 + int(r[0][3:5]) for r in rows])
    cols = {}
    for j, tag in enumerate(head[1:], 1):
        vals = []
        for r in rows:
            v = r[j] if j < len(r) else ""
            vals.append(float(v) if v not in ("", "nan", "NaN") else np.nan)
        cols[tag] = np.array(vals)
    return minutes, cols


def read_log(path):
    text = open(path).read()
    windows = []
    for m in re.finditer(r"^Test\s+(\d+)\s+(\d\d):(\d\d)-(\d\d):(\d\d)", text, re.M):
        windows.append((int(m.group(2)) * 60 + int(m.group(3)), int(m.group(4)) * 60 + int(m.group(5))))
    post = {}
    for m in re.finditer(r"^(\S+) post-test check: reads ([0-9.]+) (bar|K) (high|low)", text, re.M):
        err = float(m.group(2)) * (1 if m.group(4) == "high" else -1)
        post[m.group(1)] = err
    dead = set(re.findall(r"^(\S+) out of service", text, re.M))
    return windows, post, dead


class Campaign:
    def __init__(self, folder):
        self.sheet = json.load(open(folder + "/unit.json"))
        self.minutes, self.raw = read_trend(folder + "/trend.csv")
        self.windows, self.post, self.dead = read_log(folder + "/test_log.txt")
        self.inst = {i["tag"]: i for i in self.sheet["instruments"]}
        self.dropped = {}

    def mask(self, w):
        return (self.minutes >= w[0]) & (self.minutes <= w[1])

    def mean(self, tag, w):
        v = self.raw[tag][self.mask(w)]
        return float(np.nanmean(v))

    def temperature(self, tag, w):
        return self.mean(tag, w) + self.inst[tag].get("certificate_correction", 0.0) - self.post.get(tag, 0.0)

    def pressure(self, tag, w):
        i = self.inst[tag]
        p = self.mean(tag, w) + i.get("certificate_correction", 0.0) - self.post.get(tag, 0.0)
        if i["unit"] == "bar(g)":
            p += self.mean("PB", w) / 100.0
        line = i.get("sensing_line", {})
        if line.get("filled") == "condensate":
            H = line["transmitter_below_tap_m"]
            amb = self.mean("TE-AMB", w)
            rho = rho_liq(p, amb)
            for _ in range(3):
                p_tap = p - rho * G * H / 1e5
                rho = rho_liq(max(p_tap, 0.05), amb)
            p = p - rho * G * H / 1e5
        return p

    def screen_pair(self, a, b, referee=None):
        """Usable tags of a duplicate pair for each window. Calibrated test thermocouples agree within a few
        tenths of a kelvin; a pair that disagrees by more, or separates during the campaign, is resolved with
        the station thermocouple and the element that departs from it is dropped."""
        live = [t for t in (a, b) if t not in self.dead]
        use = {w: [t for t in live if np.isfinite(self.mean(t, w))] for w in self.windows}
        if len(live) < 2 or referee is None:
            return use
        d = np.array([self.temperature(a, w) - self.temperature(b, w) for w in self.windows])
        if np.max(np.abs(d)) < 0.6 and np.ptp(d) < 0.3:
            return use
        score = {}
        for t in (a, b):
            off = np.array([self.temperature(t, w) - self.mean(referee, w) for w in self.windows])
            score[t] = abs(float(np.mean(off))) + float(np.ptp(off))
        bad = max(score, key=score.get)
        self.dropped[bad] = "all windows"
        return {w: [t for t in use[w] if t != bad] for w in self.windows}

    def enthalpy_two(self, tags, p, w):
        return float(np.mean([h_pT(p, self.temperature(t, w)) for t in tags]))


def analyse(folder):
    c = Campaign(folder)
    turbine = c.sheet["turbine"]
    partial = "control-stage" in turbine["mid_span_packing"]
    downstream = turbine["cold_reheat_thermowells"].startswith("downstream")
    ipx_use = c.screen_pair("TE-IPX-A", "TE-IPX-B", referee="TI-IPX")
    crh_use = c.screen_pair("TE-CRH-A", "TE-CRH-B")
    pts = []
    for w in c.windows:
        p_ms = c.pressure("PT-MS", w)
        h_ms = c.enthalpy_two(["TE-MS-A", "TE-MS-B"], p_ms, w)
        p_crh = c.pressure("PT-CRH", w)
        h_crh_meas = c.enthalpy_two(crh_use[w], p_crh, w)
        p_hrh = c.pressure("PT-HRH", w)
        h_hrh = c.enthalpy_two(["TE-HRH-A", "TE-HRH-B"], p_hrh, w)
        p_x = c.pressure("PT-IPX", w)
        h_x = c.enthalpy_two(ipx_use[w], p_x, w)
        if partial:
            p_fs = c.pressure("PT-FS", w)
            h_leak = h_pT(p_fs, c.temperature("TE-FS", w))
        else:
            h_leak = h_ms
        m_fw = c.mean("FT-FW", w)
        m_sp = c.mean("FT-RHS", w)
        h_w = h_pT(c.mean("PT-RHS", w), c.mean("TE-RHS", w))
        pts.append(dict(p_ms=p_ms, h_ms=h_ms, p_crh=p_crh, h_crh_meas=h_crh_meas, p_hrh=p_hrh, h_hrh=h_hrh,
                        p_x=p_x, h_x=h_x, h_leak=h_leak, m_fw=m_fw, m_sp=m_sp, h_w=h_w))

    def ip_residuals(v):
        eta, r = v
        out = []
        for q in pts:
            h_in = (q["h_hrh"] + r * q["h_leak"]) / (1 + r)
            hs = h_ps(q["p_x"], s_ph(q["p_hrh"], h_in))
            out.append(h_in - eta * (h_in - hs) - q["h_x"])
        return out

    fit = least_squares(ip_residuals, [0.9, 0.02], x_scale=[0.01, 0.01], diff_step=1e-6)
    eta_ip, r = fit.x

    effs = []
    for q in pts:
        h_crh = q["h_crh_meas"]
        if downstream and q["m_sp"] > 0:
            m_leak = r * (q["m_fw"] + q["m_sp"]) / (1 + r)
            m_steam = q["m_fw"] - m_leak
            h_crh = q["h_crh_meas"] + q["m_sp"] / m_steam * (q["h_crh_meas"] - q["h_w"])
        hs = h_ps(q["p_crh"], s_ph(q["p_ms"], q["h_ms"]))
        effs.append((q["h_ms"] - h_crh) / (q["h_ms"] - hs))
    eta_hp = float(np.mean(effs))
    diag = dict(leak_ratio=float(r), ip_rms_kJ=float(np.sqrt(np.mean(np.square(fit.fun)))),
                hp_by_point=[float(e) for e in effs], dropped=c.dropped)
    return dict(hp_efficiency=round(eta_hp, 6), ip_efficiency=round(float(eta_ip), 6)), diag


def main():
    ans, diag = analyse(sys.argv[1])
    with open(sys.argv[2], "w") as fh:
        json.dump(ans, fh)
    print(json.dumps(dict(ans, **diag)))


if __name__ == "__main__":
    main()
