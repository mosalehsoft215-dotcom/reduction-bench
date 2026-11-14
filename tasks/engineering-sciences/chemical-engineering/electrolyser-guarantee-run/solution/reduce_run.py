#!/usr/bin/env python3
"""Reduce one membrane electrolyser guarantee run to caustic current efficiency and DC energy.

    python3 reduce_run.py <line folder> <answer.json>

Current: the DCCT channels, checked against two independent physical routes, the AC revenue meter
less the factory rectifier losses (divided by the rectifier output voltage) and the hydrogen meter
through Faraday's law (all cathode current makes hydrogen), its water vapour re-evaluated at the
measured gas temperature. Voltage: VT-E, checked against the rectifier output voltage less the
busbar drop computed from the busbar geometry. NaOH: product mass (FT-310 volume x DT-311 density)
times the NaOH fraction from the density, plus the change of NaOH held in the catholyte loop; the
product meter is proved against the tank whenever a product stop lets the level rise freely.
"""
import csv
import datetime as dt
import json
import math
import sys
from pathlib import Path

import numpy as np

FARADAY = 96485.33212
M_NAOH = 39.99711e-3
V_NORMAL = 22.41397e-3

CORRECT = dict(current="checked", voltage="checked", fraction="density", prove_meter=True,
               inventory="loop", h2_vapour="measured")


def caustic_density(w, t_c):
    return 981.4 + 1255.5 * w - 335.0 * w * w - (0.36 + 0.75 * w) * (t_c - 20.0)


def caustic_fraction(rho, t_c):
    w = np.full(np.shape(rho), 0.32)
    for _ in range(25):
        w = w - (caustic_density(w, t_c) - rho) / (1255.5 - 670.0 * w - 0.75 * (t_c - 20.0))
    return w


def water_vapour_kpa(t_c):
    return 0.1333224 * 10.0 ** (8.07131 - 1730.63 / (233.426 + t_c))


def minutes_of(text):
    t = dt.datetime.strptime(text, "%Y-%m-%dT%H:%M:%SZ")
    return (t - dt.datetime(2000, 1, 1)).total_seconds() / 60.0


def read_table(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def column(rows, name):
    return np.array([float(r[name]) if r.get(name, "") not in ("", None) else np.nan for r in rows])


class Line:
    def __init__(self, folder):
        f = Path(folder)
        self.sheet = json.loads((f / "sheet.json").read_text(encoding="utf-8"))
        self.plan = json.loads((f / "run_plan.json").read_text(encoding="utf-8"))
        hist = read_table(f / "historian.csv")
        self.t = np.array([minutes_of(r["time_utc"]) for r in hist])
        self.h = {k: column(hist, k) for k in hist[0] if k != "time_utc"}
        self.lab = read_table(f / "lab.csv")
        self.log = read_table(f / "oplog.csv")
        self.checks = read_table(f / "checks.csv")
        self.t0 = minutes_of(self.plan["window_start_utc"])
        self.t1 = minutes_of(self.plan["window_end_utc"])
        self.win = (self.t >= self.t0) & (self.t < self.t1)
        el = self.sheet["electrolyser"]
        self.n_el = el["elements"]
        self.holdup = el["catholyte_holdup_outside_tank_m3"]
        tank = self.sheet["catholyte_tank"]
        self.tank_area = math.pi * tank["inside_diameter_m"] ** 2 / 4.0
        self.rho_cal = tank["LT-330"]["calibration_density_kg_m3"]
        text = self.sheet["instruments"]["FT-410"]
        self.t_cfg = float(text.split("saturated at ")[1].split(" degC")[0])

    def has(self, tag):
        return np.isfinite(self.h[tag][self.win]).mean() > 0.9

    # ---------------------------------------------------------------- current
    def loss_kw(self, i_ka):
        pts = self.sheet["rectifier"]["factory_loss_test"]
        x = np.array([p["dc_current_kA"] for p in pts])
        y = np.array([p["losses_kW"] for p in pts])
        c = np.polyfit(x, y, 2)
        return np.polyval(c, i_ka)

    def current_ac(self):
        p_kw = self.h["JT_AC_MW"] * 1e3
        v = self.h["VT_R_V"]
        i = p_kw / v
        for _ in range(4):
            i = (p_kw - self.loss_kw(i)) / v
        return i

    def current_h2(self, vapour):
        q = self.h["FT410_Nm3h"]
        p = self.h["PT411_kPa"]
        if vapour == "measured":
            q = q * (p - water_vapour_kpa(self.h["TT412_C"])) / (p - water_vapour_kpa(self.t_cfg))
        mol_s = q / 3600.0 / V_NORMAL
        return 2.0 * FARADAY * mol_s / self.n_el / 1e3

    def current(self, mode, vapour):
        w = self.win
        chans = [self.h[k] for k in ("DCCT_A_kA", "DCCT_B_kA") if self.has(k)]
        if mode.startswith("dcct"):
            return chans[0] if mode == "dcct_a" else np.nanmean(np.vstack(chans), axis=0)
        dcct = np.nanmean(np.vstack(chans), axis=0)
        if mode == "h2_first" and self.has("FT410_Nm3h"):
            return self.current_h2(vapour)
        if mode == "ac_first" and self.has("JT_AC_MW"):
            return self.current_ac()
        refs = []
        if self.has("JT_AC_MW"):
            refs.append(np.nanmean(self.current_ac()[w]) / np.nanmean(dcct[w]))
        if self.has("FT410_Nm3h"):
            refs.append(np.nanmean(self.current_h2(vapour)[w]) / np.nanmean(dcct[w]))
        if mode == "vote":
            # majority of the individual sources (each DCCT channel counts once)
            votes = [1.0] * len(chans) + refs
            r = float(np.median(votes))
            return dcct * (r if abs(r - 1.0) > 0.004 else 1.0)
        if not refs or all(abs(r - 1.0) < 0.004 for r in refs):
            return dcct
        r = float(np.mean(refs))
        if len(chans) == 2:
            single = [np.nanmean(c[w]) / np.nanmean(dcct[w]) for c in chans]
            best = int(np.argmin([abs(s - r) for s in single]))
            if abs(single[best] - r) < 0.004:
                return chans[best]
        return dcct * r

    # ---------------------------------------------------------------- voltage
    def voltage(self, mode, i_ka):
        vt_e = self.h["VT_E_V"]
        vt_r = self.h["VT_R_V"]
        if mode == "vt_e":
            return vt_e
        if mode == "vt_r":
            return vt_r
        if mode == "fitted_bus":
            drop = vt_r - vt_e
            r = np.nanmean(drop[self.win]) / np.nanmean(i_ka[self.win])
            return vt_r - r * i_ka
        # The busbar between VT-R and VT-E is a resistance: its drop vanishes with the current. With
        # VT-E = g U_E and VT-R = U_E + R(T) I, the drop VT-R - VT-E = (1 - g) VT-R + g R(T) I, so
        # two clearly different currents give g; one current level leaves g at 1.
        alpha = self.sheet["busbar"]["temperature_coefficient_per_K"]
        x = i_ka * 1e3 * (1.0 + alpha * (self.h["TT_BUS_C"] - 20.0))
        ok = np.isfinite(vt_e) & np.isfinite(vt_r) & np.isfinite(x)
        low = ok & (i_ka < 0.85 * np.nanmedian(i_ka[self.win]))
        if low.sum() < 30:
            return vt_e
        a, _b = np.linalg.lstsq(np.column_stack([vt_r[ok], x[ok]]), (vt_r - vt_e)[ok], rcond=None)[0]
        return vt_e / (1.0 - a) if abs(a) > 0.004 else vt_e

    # ---------------------------------------------------------------- caustic
    def fraction(self, mode):
        if mode == "density":
            return caustic_fraction(self.h["DT311_kgm3"], self.h["TT311_C"])
        if mode == "analyser":
            return self.h["AT311_pct"] / 100.0
        rows = [r for r in self.lab if r["status"] == "VALID" and r["point"] == "SP-310"]
        tl = np.array([minutes_of(r["time_utc"]) for r in rows])
        wl = np.array([float(r["naoh_pct"]) for r in rows]) / 100.0
        if mode == "lab_controlled":         # each lot's titrations scaled by certified / found on the control standard
            ctrl = {r["titrant_lot"]: float(r["naoh_pct"]) for r in self.lab if r["point"] == "CTRL"}
            cert = self.sheet["control_standard"]["naoh_pct"]
            wl = wl * np.array([cert / ctrl.get(r["titrant_lot"], cert) for r in rows])
            return np.interp(self.t, tl, wl)
        if mode == "lab_mean":
            sel = (tl >= self.t0) & (tl <= self.t1)
            return np.full(self.t.shape, wl[sel].mean())
        return np.interp(self.t, tl, wl)

    def stops(self):
        q = self.h["FT310_m3h"]
        low = q < 0.05 * np.nanmedian(q[self.win])
        out, k = [], 0
        while k < len(q):
            if low[k]:
                e = k
                while e < len(q) and low[e]:
                    e += 1
                if e - k >= 10:
                    out.append((k, e))
                k = e
            else:
                k += 1
        return out

    def meter_scale(self):
        m_tank = self.h["LT330_mm"] / 1000.0 * self.tank_area * self.rho_cal
        mass_rate = self.h["FT310_m3h"] * self.h["DT311_kgm3"] / 3600.0
        dm = self.h["FT320_kgh"] / 3600.0
        scales = []
        for a, b in self.stops():
            if a < 100 or b - a < 15:
                continue
            ins = np.arange(a + 2, b - 2)
            pre = np.arange(a - 95, a - 5)
            rise = np.polyfit(self.t[ins] * 60.0, m_tank[ins], 1)[0]
            drift = np.polyfit(self.t[pre] * 60.0, m_tank[pre], 1)[0]
            true_pre = rise - (dm[ins].mean() - dm[pre].mean()) - drift
            scales.append(mass_rate[pre].mean() / true_pre)
        if not scales:
            return 1.0
        s = float(np.mean(scales))
        return s if abs(s - 1.0) > 0.005 else 1.0

    def naoh_kg(self, frac_mode, prove, inventory):
        w = self.fraction(frac_mode)
        rate = self.h["FT310_m3h"] * self.h["DT311_kgm3"] / 3600.0
        if prove:
            rate = rate / self.meter_scale()
        win = self.win
        product = np.nansum(rate[win] * w[win]) * 60.0
        if inventory == "none":
            return product
        m_tank = self.h["LT330_mm"] / 1000.0 * self.tank_area * self.rho_cal
        if inventory == "tank":
            loop = m_tank * w
        else:
            hold = self.holdup * caustic_density(w, self.h["TT_CATH_C"])
            loop = (m_tank + hold) * w
        idx = np.flatnonzero(win)
        start = np.nanmean(loop[max(idx[0] - 2, 0): idx[0] + 3])
        end = np.nanmean(loop[idx[-1] - 2: idx[-1] + 3])
        return product + (end - start)

    def reduce(self, routes=None):
        r = dict(CORRECT)
        r.update(routes or {})
        i_ka = self.current(r["current"], r["h2_vapour"])
        naoh = self.naoh_kg(r["fraction"], r["prove_meter"], r["inventory"])
        charge = self.n_el * np.nansum(i_ka[self.win]) * 1e3 * 60.0 / FARADAY
        u = self.voltage(r["voltage"], i_ka)
        energy = np.nansum(u[self.win] * i_ka[self.win] * 1e3) * 60.0
        return {"current_efficiency_pct": round(100.0 * naoh / M_NAOH / charge, 4),
                "dc_energy_kwh_per_t": round(energy / naoh / 3600.0, 3)}


def main():
    folder, answer = sys.argv[1], sys.argv[2]
    result = Line(folder).reduce()
    Path(answer).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
