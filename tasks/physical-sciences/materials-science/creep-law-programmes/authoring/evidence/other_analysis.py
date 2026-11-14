"""A second, separately written reduction of a creep programme (shares no code with solution/creep_constants.py).

Differences from the reference: records read with numpy.genfromtxt; log read by keyword tests sentence by sentence;
the specimen described by Gauss-Legendre nodes on the parallel length and on each fillet arc (shoulder ignored);
local rates from quadratic fits over 21 readings centred on every 7th reading; the steady part chosen in time
(after twice the time of the smoothed minimum rate); the law fitted by Gauss-Newton with analytic derivatives.

usage: python3 other_analysis.py <programme folder> [<result.json>]"""
import json
import os
import sys
from datetime import datetime

import numpy as np

RGAS = 8.314462618
NODES, WEIGHTS = np.polynomial.legendre.leggauss(24)


def clock(text):
    return datetime.strptime(text, "%d.%m.%Y %H:%M")


def numbers_after(text, word):
    i = text.find(word)
    return [] if i < 0 else [float(v) for v in __import__("re").findall(r"\d+(?:\.\d+)?", text[i:])]


class Programme:
    def __init__(self, folder):
        self.folder = folder
        self.sheet = json.load(open(os.path.join(folder, "programme.json"), encoding="utf-8"))
        rows = np.genfromtxt(os.path.join(folder, "tests.csv"), delimiter=",", dtype=str, encoding="utf-8")
        head = list(rows[0])
        self.tests = [dict(zip(head, r)) for r in rows[1:]]
        self.entries = []
        for line in open(os.path.join(folder, "log.txt"), encoding="utf-8"):
            if len(line) > 18 and line[2] == "." and line[16:18] == "  ":
                self.entries.append((clock(line[:16]), line[18:].strip()))
        self.read_entries()

    def read_entries(self):
        import re
        self.tc_bias, self.extra_kg, self.lever, self.replaced, self.knocks, self.stopped = {}, {}, {}, [], [], set()
        start = {t["test"]: clock(t["loaded"]) for t in self.tests}
        for when, text in self.entries:
            for sentence in re.split(r"(?<=\.)\s+", text):
                pass
            tests = set(re.findall(r"T\d\d", text))
            low = text.lower()
            if " reads " in low and " k " in low + " ":
                m = re.search(r"(K-\d+).*? reads (\S+) K (high|low)", text)
                for t in tests:
                    self.tc_bias[(t, m.group(1))] = float(m.group(2)) * (1 if m.group(3) == "high" else -1)
            if "than the test sheet" in low:
                m = re.search(r"(\S+) kg (more|less) than", text)
                for t in tests:
                    self.extra_kg[t] = float(m.group(1)) * (1 if m.group(2) == "more" else -1)
            if low.startswith("lever ratio"):
                value = numbers_after(text, ":")[0]
                for t in tests:
                    self.lever[t] = value
            if " replaced by " in low:
                t = text.split(":")[0]
                old, new = re.search(r"(K-\d+) replaced by (K-\d+)", text).groups()
                self.replaced.append((t, old, new, (when - start[t]).total_seconds() / 3600))
            if "extensometer" in low and "seated" in low:
                t = text.split(":")[0]
                self.knocks.append((t, (when - start[t]).total_seconds() / 3600))
            if "burnt out" in low or "not restarted" in low:
                self.stopped.add(text.split(":")[0])


def corrected(prog, test, col_name, tc, t, raw):
    cert_ids = np.array([tc] * len(t), dtype=object)
    for tid, old, new, h in prog.replaced:
        if tid == test and old == tc:
            cert_ids[t >= h] = new
    out = raw.copy()
    for cid in set(cert_ids):
        c = prog.sheet["thermocouples"][cid]
        s = cert_ids == cid
        out[s] = raw[s] + np.interp(raw[s], c["reading_C"], c["correction_K"]) - prog.tc_bias.get((test, cid), 0.0)
    return out


def specimen_nodes(d, Lp, rho, D):
    """Positions, axial weights and cold radii at quadrature nodes: parallel length and the two fillet arcs."""
    rg, rs = d / 2, D / 2
    xf = np.sqrt(rho ** 2 - (rho - (rs - rg)) ** 2)
    xp = 0.5 * Lp * NODES
    wp = 0.5 * Lp * WEIGHTS
    u = 0.5 * xf * (NODES + 1)
    wu = 0.5 * xf * WEIGHTS
    rf = rg + rho - np.sqrt(rho ** 2 - u ** 2)
    x = np.concatenate([xp, Lp / 2 + u, -Lp / 2 - u])
    w = np.concatenate([wp, wu, wu])
    r = np.concatenate([np.full_like(xp, rg), rf, rf])
    return x, w, r


def local_rates(t, y, half=10, stride=7):
    out = []
    for c in range(half, len(t) - half, stride):
        sl = slice(c - half, c + half + 1)
        tt = t[sl] - t[c]
        if np.ptp(tt) > (2 * half + 0.5) * np.median(np.diff(t[sl])):
            continue
        p = np.polyfit(tt, y[sl], 2)
        out.append((c, p[1]))
    return out


def reduce(folder):
    prog = Programme(folder)
    sh = prog.sheet
    a_tab = sh["material"]["expansion_mean_from_20C"]
    alpha = lambda T: np.interp(T, a_tab["temperature_C"], a_tab["per_K"])
    obs = []
    for test in prog.tests:
        tid = test["test"]
        if tid in prog.stopped:
            continue
        data = np.genfromtxt(os.path.join(folder, "records", tid + ".csv"), delimiter=",", names=True,
                             encoding="utf-8", missing_values="", filling_values=np.nan)
        t, ext = data["time_h"], data["extension_mm"].copy()
        for kt, h in prog.knocks:
            if kt == tid:
                j = np.searchsorted(t, h, side="right")
                left = np.polyfit(t[j - 6:j], ext[j - 6:j], 1)
                right = np.polyfit(t[j:j + 6], ext[j:j + 6], 1)
                ext[j:] -= np.polyval(right, h) - np.polyval(left, h)
                t_mark = h
        tcs = [corrected(prog, tid, c, test[k], t, data[c]) for c, k in
               (("tc_lower_C", "tc_lower"), ("tc_middle_C", "tc_middle"), ("tc_upper_C", "tc_upper"))]
        tcs = np.vstack(tcs)
        late_med = np.nanmedian(tcs[:, len(t) // 3:], axis=1)
        calm = np.all(np.abs(tcs - late_med[:, None]) < 0.5, axis=0) & np.all(np.isfinite(tcs), axis=0)
        rates = local_rates(t, ext)
        good = []
        for c, r in rates:
            lo, hi = max(c - 10, 0), min(c + 11, len(t))
            if not calm[lo:hi].all():
                continue
            if any(kt == tid and t[lo] <= h <= t[hi - 1] for kt, h in prog.knocks):
                continue
            good.append((c, r))
        if len(good) < 5:
            continue
        rr = np.array([r for _, r in good])
        sm = np.convolve(rr, np.ones(5) / 5, mode="valid")
        k = int(np.argmin(sm)) + 2
        if k > 0.85 * len(good):
            continue
        t_min = t[good[k][0]]
        dr = sh["drawings"][test["drawing"]]
        x, w, r0 = specimen_nodes(float(test["diameter_mm"]), float(test["parallel_length_mm"]),
                                  dr["fillet_radius_mm"], dr["shoulder_diameter_mm"])
        pan = float(test["pan_mass_kg"]) + prog.extra_kg.get(tid, 0.0)
        lever = prog.lever.get(tid, sh["machines"][test["machine"]]["lever_ratio"])
        F = pan * sh["local_gravity_m_s2"] * lever
        h = float(test["parallel_length_mm"]) / 2
        for c, r in good:
            if t[c] < 2.0 * t_min:
                continue
            lo, mi, up = tcs[:, max(c - 10, 0):c + 11].mean(axis=1)
            # Lagrange form of the parabola through (-h, lo), (0, mi), (h, up)
            T = lo * x * (x - h) / (2 * h * h) - mi * (x + h) * (x - h) / (h * h) + up * x * (x + h) / (2 * h * h)
            grow = 1 + alpha(T) * (T - 20.0)
            obs.append(dict(r=r, F=F, w=w * grow, area=np.pi * (r0 * grow) ** 2, TK=T + 273.15,
                            creep=ext[c] - ext[0], x=x))
    return prog, obs


def fit(prog, obs):
    s_ref = np.median([o["F"] / o["area"].min() for o in obs])
    T_ref = np.median([o["TK"][np.argmin(np.abs(o["x"]))] for o in obs])
    p = np.array([np.log(np.median([o["r"] for o in obs]) / 40.0 / 3600.0), 5.0, 300.0])
    for outer in range(5):
        lnr, n, Qk = p
        for o in obs:   # strains of the nodes for the current law: proportional to the local steady rate
            i0 = np.argmin(np.abs(o["x"]))
            rel = (o["area"][i0] / o["area"]) ** n * np.exp(-Qk * 1e3 / RGAS * (1 / o["TK"] - 1 / o["TK"][i0]))
            lo_, hi_ = 0.0, 0.5
            for _ in range(50):
                m = 0.5 * (lo_ + hi_)
                lo_, hi_ = (m, hi_) if np.sum(o["w"] * (np.exp(rel * m) - 1)) < o["creep"] else (lo_, m)
            o["e"] = rel * m
        for it in range(40):
            J, res = [], []
            for o in obs:
                s = o["F"] * np.exp(o["e"]) / o["area"]
                g = np.exp(lnr + n * np.log(s / s_ref) - Qk * 1e3 / RGAS * (1 / o["TK"] - 1 / T_ref))
                terms = o["w"] * np.exp(o["e"]) * g * 3600.0
                tot = terms.sum()
                res.append(np.log(tot) - np.log(o["r"]))
                J.append([1.0, np.sum(terms * np.log(s / s_ref)) / tot,
                          np.sum(terms * (-1e3 / RGAS) * (1 / o["TK"] - 1 / T_ref)) / tot])
            J, res = np.array(J), np.array(res)
            step = np.linalg.solve(J.T @ J, -J.T @ res)
            p = p + step
            lnr, n, Qk = p
            if np.max(np.abs(step / np.array([1, 1, 100]))) < 1e-10:
                break
    Td = prog.sheet["design_point"]["temperature_C"] + 273.15
    sd = prog.sheet["design_point"]["stress_MPa"]
    lnr, n, Qk = p
    rate = float(np.exp(lnr + n * np.log(sd / s_ref) - Qk * 1e3 / RGAS * (1 / Td - 1 / T_ref)))
    return {"stress_exponent": float(n), "activation_energy_kJ_mol": float(Qk), "design_creep_rate_per_s": rate,
            "points": len(obs), "rms": float(np.sqrt(np.mean(res ** 2)))}


def analyse(folder):
    return fit(*reduce(folder))


if __name__ == "__main__":
    result = analyse(sys.argv[1])
    if len(sys.argv) > 2:
        json.dump(result, open(sys.argv[2], "w"), indent=1)
    print(json.dumps(result))
