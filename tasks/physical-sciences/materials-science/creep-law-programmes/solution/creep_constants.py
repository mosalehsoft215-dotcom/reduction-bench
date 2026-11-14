#!/usr/bin/env python3
"""Reduce one constant-load creep programme to the steady-state creep law of its material.

usage: python3 creep_constants.py <programme folder> <result.json>

Writes stress_exponent, activation_energy_kJ_mol and design_creep_rate_per_s (true strain rate under the design
true stress at the design temperature).
"""
import csv
import datetime as dt
import json
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

R = 8.314462618
WINDOW = 16          # readings per rate window
T_STEADY = 0.4       # K: a window is used only while every thermocouple stays this close to its test median
LATE = 1.8           # windows beyond this multiple of the strain at the observed minimum rate


# ----------------------------------------------------------------------------------------------------------- input
def load_programme(folder):
    with open(os.path.join(folder, "programme.json"), encoding="utf-8") as fh:
        sheet = json.load(fh)
    with open(os.path.join(folder, "tests.csv"), encoding="utf-8", newline="") as fh:
        tests = list(csv.DictReader(fh))
    log = []
    with open(os.path.join(folder, "log.txt"), encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r"\s*(\d{2}\.\d{2}\.\d{4} \d{2}:\d{2})\s+(.*\S)", line)
            if m:
                log.append((dt.datetime.strptime(m.group(1), "%d.%m.%Y %H:%M"), m.group(2)))
    records = {}
    for t in tests:
        with open(os.path.join(folder, "records", t["test"] + ".csv"), encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh))
        num = lambda v: float(v) if v not in ("", None) else np.nan
        records[t["test"]] = {k: np.array([num(r[k]) for r in rows]) for k in rows[0]}
    return sheet, tests, log, records


def read_log(log, tests):
    """Turn the operator's entries into corrections and exclusions."""
    loaded = {t["test"]: dt.datetime.strptime(t["loaded"], "%d.%m.%Y %H:%M") for t in tests}
    hours = lambda tid, when: (when - loaded[tid]).total_seconds() / 3600.0
    out = dict(tc_error={}, pan_delta={}, ratio={}, swaps=[], blackouts=[], steps=[], ended_early=set())
    for when, text in log:
        named = re.findall(r"\bT\d{2}\b", text)
        for tc, val, sense in re.findall(r"(K-\d+)[^.]*?reads (\d+(?:\.\d+)?) K (high|low)", text):
            for tid in named:
                out["tc_error"][(tid, tc)] = float(val) if sense == "high" else -float(val)
        m = re.search(r"(\d+(?:\.\d+)?) kg (more|less) than the test sheet", text)
        if m:
            for tid in named:
                out["pan_delta"][tid] = float(m.group(1)) * (1 if m.group(2) == "more" else -1)
        m = re.search(r"[Ll]ever ratio of (\S+?)[^:]*:\s*(\d+(?:\.\d+)?)", text)
        if m:
            for tid in named:
                out["ratio"][tid] = float(m.group(2))
        m = re.search(r"(T\d{2}): (K-\d+) replaced by (K-\d+)", text)
        if m:
            out["swaps"].append((m.group(1), m.group(2), m.group(3), hours(m.group(1), when)))
        if re.search(r"mains failure", text, re.I):
            out["blackouts"].append([when, None])
        if re.search(r"mains restored", text, re.I) and out["blackouts"]:
            out["blackouts"][-1][1] = when
        m = re.search(r"(T\d{2}): extensometer .*re-?seated", text)
        if m:
            out["steps"].append((m.group(1), hours(m.group(1), when)))
        m = re.search(r"(T\d{2}): .*(?:test ended|not restarted)", text)
        if m:
            out["ended_early"].add(m.group(1))
    return out


# ------------------------------------------------------------------------------------------------------- specimen
def slices(d, Lp, rho, D, span, n_par=160, n_fil=80, n_sh=6):
    rg, rs = d / 2, D / 2
    xf = np.sqrt(rho ** 2 - (rho - (rs - rg)) ** 2)
    sh = max((span - Lp) / 2 - xf, 1e-6)
    edges = np.concatenate([np.linspace(0, Lp / 2, n_par // 2 + 1), Lp / 2 + np.linspace(0, xf, n_fil + 1)[1:],
                            Lp / 2 + xf + np.linspace(0, sh, n_sh + 1)[1:]])
    mid, du = 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)
    w = mid - Lp / 2
    r = np.where(w <= 0, rg, np.where(w < xf, rg + rho - np.sqrt(np.clip(rho ** 2 - w ** 2, 0, None)), rs))
    return np.concatenate([-mid[::-1], mid]), np.concatenate([du[::-1], du]), np.pi * np.concatenate([r[::-1], r]) ** 2


# --------------------------------------------------------------------------------------------------- one test
def windows_of(test, sheet, notes, rec):
    """Rate windows of one test: (extension rate mm/h, creep extension mm, three thermocouple temperatures C)."""
    tid = test["test"]
    t, ext = rec["time_h"], rec["extension_mm"].copy()
    certs = sheet["thermocouples"]
    names = [test["tc_lower"], test["tc_middle"], test["tc_upper"]]
    temps = []
    for col, name in zip(("tc_lower_C", "tc_middle_C", "tc_upper_C"), names):
        raw = rec[col]
        which = np.full(len(t), name, dtype=object)
        for sid, old, new, h in notes["swaps"]:
            if sid == tid and old == name:
                which[t >= h] = new
        true = np.full(len(t), np.nan)
        for tc in set(which):
            c = certs[tc]
            sel = which == tc
            true[sel] = raw[sel] + np.interp(raw[sel], c["reading_C"], c["correction_K"]) - notes["tc_error"].get((tid, tc), 0.0)
        temps.append(true)
    temps = np.array(temps)
    # extensometer steps reported in the log: remove them from the record
    for sid, h in notes["steps"]:
        if sid != tid:
            continue
        i1 = int(np.searchsorted(t, h, side="right"))
        before, after = slice(max(i1 - 9, 0), i1), slice(i1, min(i1 + 9, len(t)))
        sb = np.polyfit(t[before], ext[before], 1)
        sa = np.polyfit(t[after], ext[after], 1)
        jump = np.polyval(sa, h) - np.polyval(sb, h)
        ext[i1:] -= jump
    ok = np.all(np.isfinite(temps), axis=0)
    med = np.nanmedian(temps[:, t > 0.2 * t[-1]], axis=1)
    ok &= np.all(np.abs(temps - med[:, None]) < T_STEADY, axis=0)
    step = np.median(np.diff(t))
    cut = np.zeros(len(t), bool)
    cut[1:] = np.diff(t) > 1.5 * step
    for sid, h in notes["steps"]:
        if sid == tid:
            cut[int(np.searchsorted(t, h, side="right"))] = True
    out = []
    i = 1
    while i + WINDOW <= len(t):
        idx = np.arange(i, i + WINDOW)
        if ok[idx].all() and not cut[idx[1:]].any():
            rate = np.polyfit(t[idx], ext[idx], 1)[0]
            out.append((t[idx].mean(), rate, ext[idx].mean() - ext[0], temps[:, idx].mean(axis=1)))
            i += WINDOW
        else:
            i += 1
    return out


class Test:
    def __init__(self, row, sheet, notes, rec):
        self.id = row["test"]
        dr = sheet["drawings"][row["drawing"]]
        self.Lp = float(row["parallel_length_mm"])
        self.x, self.dx, self.a = slices(float(row["diameter_mm"]), self.Lp, dr["fillet_radius_mm"],
                                         dr["shoulder_diameter_mm"], dr["ridge_spacing_mm"])
        pan = float(row["pan_mass_kg"]) + notes["pan_delta"].get(self.id, 0.0)
        ratio = notes["ratio"].get(self.id, sheet["machines"][row["machine"]]["lever_ratio"])
        self.F = pan * sheet["local_gravity_m_s2"] * ratio
        exp_tab = sheet["material"]["expansion_mean_from_20C"]
        self.alpha = lambda T: np.interp(T, exp_tab["temperature_C"], exp_tab["per_K"])
        self.windows = windows_of(row, sheet, notes, rec)
        self.early = self.id in notes["ended_early"]

    def profile(self, tc):
        lo, mi, up = tc
        h = self.Lp / 2
        return mi + (up - lo) / (2 * h) * self.x + (lo + up - 2 * mi) / (2 * h * h) * self.x ** 2

    def prepare(self, n, Q):
        """Per window: slice temperatures, hot slice lengths and sections, slice strains."""
        self.rows = []
        for tm, rate, creep, tc in self.windows:
            T = self.profile(tc)
            k = 1 + self.alpha(T) * (T - 20.0)
            dx, a = self.dx * k, self.a * k ** 2
            TK = T + 273.15
            i0 = np.argmin(np.abs(self.x))
            w = (a[i0] / a) ** n * np.exp(-Q / R * (1 / TK - 1 / TK[i0]))
            lo, hi = 0.0, 0.5
            for _ in range(60):
                eg = 0.5 * (lo + hi)
                if np.sum(dx * (np.exp(w * eg) - 1)) > creep:
                    hi = eg
                else:
                    lo = eg
            self.rows.append(dict(t=tm, rate=rate, dx=dx, a=a, TK=TK, e=w * eg, eg=eg))


def model_rate(row, F, lnr0, n, Q, s_ref, T_ref):
    s = F * np.exp(row["e"]) / row["a"]
    de = np.exp(lnr0 + n * np.log(s / s_ref) - Q / R * (1 / row["TK"] - 1 / T_ref))
    return np.sum(row["dx"] * np.exp(row["e"]) * de) * 3600.0


def select_late(test, n, Q, s_ref, T_ref):
    rows = test.rows
    if len(rows) < 4:
        return []
    raw = np.array([r["rate"] for r in rows])
    k_min = int(np.argmin(raw))
    if test.early or k_min >= 0.85 * (len(rows) - 1):
        return []
    e_min = rows[k_min]["eg"]
    return [r for r in rows if r["eg"] >= LATE * e_min]


def solve(folder):
    sheet, tests, log, records = load_programme(folder)
    notes = read_log(log, tests)
    objs = [Test(row, sheet, notes, records[row["test"]]) for row in tests]
    s_ref = float(np.median([o.F / o.a.min() for o in objs]))
    T_ref = float(np.median([float(r["setpoint_C"]) for r in tests])) + 273.15
    n, Q, lnr0 = 5.0, 300e3, np.log(1e-8)
    for _ in range(4):
        used = []
        for o in objs:
            o.prepare(n, Q)
            used += [(o, r) for r in select_late(o, n, Q, s_ref, T_ref)]
        y = np.log([r["rate"] for _, r in used])

        def resid(p):
            return np.array([np.log(model_rate(r, o.F, p[0], p[1], p[2] * 1e3, s_ref, T_ref)) for o, r in used]) - y

        fit = least_squares(resid, [lnr0, n, Q / 1e3], x_scale=[1.0, 0.5, 20.0])
        lnr0, n, Q = fit.x[0], fit.x[1], fit.x[2] * 1e3
    Td = sheet["design_point"]["temperature_C"] + 273.15
    sd = sheet["design_point"]["stress_MPa"]
    rate_d = float(np.exp(lnr0 + n * np.log(sd / s_ref) - Q / R * (1 / Td - 1 / T_ref)))
    rms = float(np.sqrt(np.mean(fit.fun ** 2)))
    return {"stress_exponent": float(n), "activation_energy_kJ_mol": float(Q / 1e3), "design_creep_rate_per_s": rate_d,
            "windows_used": len(used), "rms_log_residual": rms}


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: creep_constants.py <programme folder> <result.json>")
    result = solve(sys.argv[1])
    with open(sys.argv[2], "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1)


if __name__ == "__main__":
    main()
