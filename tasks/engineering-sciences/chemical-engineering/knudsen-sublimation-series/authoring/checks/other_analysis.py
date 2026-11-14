"""A separately written analysis of one effusion series, used only to check that the limits are fair.

    python3 other_analysis.py <series folder> <answer.json>

Written from formats.txt alone; nothing is shared with solution/. Own reading of every file and of the log (word
rules), temperatures resampled to a 2 s grid and integrated with Simpson's rule, the waits taken at the record's
first and last reading, a crystal in a liner cooled by the heat its sublimation takes and by the black-body
radiation of the liner's orifice to the 77 K shroud (Crank-Nicolson on the grid), the certificate interpolated by a quadratic through its three nearest rows in (1/T, ln p),
and every unknown (dH298, ln p298, ln C of each foil, the offset of each uncalibrated spare thermometer, an unknown
lead resistance) found by a hand-written Levenberg-Marquardt iteration on ln(mass loss) with finite-difference
derivatives.
"""
import csv
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.signal import lfilter

GAS = 8.314462618
STEFAN = 5.670374419e-8
SHROUD_K = 77.0
TREF = 298.15


def form1(T, h, lp, dc):
    return lp - h / GAS * (1 / T - 1 / TREF) + dc / GAS * (TREF / T - 1 + np.log(T / TREF))


def form2_shift(T, Tt, h_t, dcp):
    # ln p(II) - ln p(I): minus the integral of the transition enthalpy over R T^2 from Tt to T
    return -((h_t - dcp * Tt) * (1 / Tt - 1 / T) + dcp * np.log(T / Tt)) / GAS


def ids_in(words, ids):
    return [w.strip(".,;:()") for w in words if w.strip(".,;:()") in ids]


def span(named, words, ids):
    """'R09 to R11' -> R09..R11; 'R03 and R07' -> both; otherwise the ids as written."""
    text = " ".join(words)
    if len(named) == 2 and f"{named[0]} to {named[1]}" in text:
        return ids[ids.index(named[0]):ids.index(named[1]) + 1]
    return named


class Series:
    def __init__(self, folder):
        self.folder = folder = Path(folder)
        self.cellj = json.loads((folder / "cell.json").read_text())
        comp = json.loads((folder / "compounds.json").read_text())
        self.runs = list(csv.DictReader(open(folder / "runs.csv", newline="")))
        self.ids = [r["run"] for r in self.runs]
        self.kind = {r["run"]: r["compound"] for r in self.runs}
        self.wg = {r["run"]: (r["before_mg"], r["after_mg"]) for r in csv.DictReader(open(folder / "weighings.csv"))}
        th = self.cellj["thermometer"]
        self.res = th["column"] == "resistance_ohm"
        self.G = self.cellj.get("liner", {}).get("conductance_W_per_K")
        self.CL = self.cellj.get("liner", {}).get("heat_capacity_J_per_K")
        if th.get("location") == "liner":      # the controller holds the liner (and the crystals) on the program
            self.G = self.CL = None
        self.empty = self.cellj.get("empty_mass_mg")
        # the orifice of the closed liner sends black-body radiation to the shroud
        self.hole = STEFAN * math.pi * (0.5e-3 * self.cellj["orifice_diameter_mm"]) ** 2
        self.rough = False      # True while starting values are looked for (no liner cooling then)
        self.two_wire = self.res and int(th.get("wires", 4)) == 2
        self.lead_known = th.get("lead_resistance_ohm") if self.two_wire else 0.0
        self.raw = {}
        for r in self.runs:
            rows = np.array([[float(a), float(b)] for a, b in list(csv.reader(open(folder / r["file"])))[1:]])
            self.raw[r["run"]] = rows
        ref = comp["reference"]
        self.Mr = ref["molar_mass_g_mol"] * 1e-3
        tab = np.array([[float(a), float(b)] for a, b in list(csv.reader(open(folder / ref["vapour_pressure_table"])))[1:]])
        self.ix, self.iy = 1.0 / (tab[:, 0] + 273.15), np.log(tab[:, 1])
        self.cp_ref = ref.get("Cp_crystal_J_K_mol")
        area = self.cellj.get("bed_area_mm2", 0.0) * 1e-6
        self.evap_cond = {"reference": ref.get("evaporation_coefficient"),
                          "sample": comp["sample"].get("evaporation_coefficient")}
        self.evap_cond = {k: (None if v is None else v * area) for k, v in self.evap_cond.items()}
        smp = comp["sample"]
        self.cp_smp = smp["Cp_crystal_J_K_mol"]
        self.Ms = smp["molar_mass_g_mol"] * 1e-3
        self.dc = smp["Cp_gas_J_K_mol"] - smp["Cp_crystal_J_K_mol"]
        tr = smp.get("transition")
        self.trans = None if tr is None else (tr["temperature_C"] + 273.15, tr["dH_kJ_mol"] * 1000.0,
                                              tr["Cp_crystal_J_K_mol"] - smp["Cp_crystal_J_K_mol"])
        self.read_log()

    def read_log(self):
        ids = self.ids
        self.corr = {i: 0.0 for i in ids}
        self.covered = {i: None for i in ids}
        self.icecheck = None
        self.spare = {i: None for i in ids}
        self.counts = {i: True for i in ids}
        self.foil = {i: 0 for i in ids}
        self.start = {}
        self.stays = []
        nfoil, nspare = 0, 0
        self.part = {}
        for line in (self.folder / "run_log.txt").read_text().splitlines():
            stamp, _, text = line.partition(" | ")
            when = datetime.strptime(stamp, "%d.%m.%Y %H:%M:%S")
            words = text.split()
            if len(words) >= 2 and words[0] in ids and words[1] in ("under", "started", "vented.", "weighed."):
                rid = words[0]
                if words[1] == "under":
                    self.stays.append({"from": when, "runs": []})
                elif words[1] == "started":
                    self.stays[-1]["runs"].append(rid)
                    self.start[rid] = when
                    self.foil[rid] = nfoil
                elif words[1] == "vented.":
                    self.stays[-1]["to"] = when
                continue
            low = text.lower()
            named = ids_in(words, ids)
            if "spare" in low and "onwards" in low:
                first = ids_in(low.split("fitted for")[-1].upper().split(), ids)[0]
                if "no certificate" in low:
                    for i in ids[ids.index(first):]:
                        self.corr[i], self.spare[i] = 0.0, nspare
                    nspare += 1
                else:
                    cert = low.split("certificate:")[1].split(")")[0].split()
                    if cert[0] == "correction":
                        c = float(cert[1])
                    else:
                        c = float(cert[1]) if cert[3] == "low" else -float(cert[1])
                    k = cert.index("between")
                    rng = (float(cert[k + 1]), float(cert[k + 3]))
                    for i in ids[ids.index(first):]:
                        self.corr[i], self.spare[i], self.covered[i] = c, None, rng
            if "not for evaluation" in low:
                for i in span(named, words, ids):
                    self.counts[i] = False
            if "repeated as" in low:
                head = text.split("repeated as")[0].split()
                for i in span(ids_in(head, ids), head, ids):
                    self.counts[i] = False
            if "ohm at the ice point and" in low:
                w = low.split()
                self.icecheck = (float(w[w.index("ohm") - 1]), float(w[len(w) - 1 - w[::-1].index("ohm") - 1]))
            if "new foil of the same nominal size fitted" in low:
                nfoil += 1
            if len(words) > 2 and words[0] in ids and words[1] == "weighed" and words[2] == "with" and "mg on its own" in low:
                before_mg = low.split("mg on its own")[0].split()
                self.part[words[0]] = float(before_mg[-1])
        self.nfoil, self.nspare = nfoil + 1, nspare
        for s in self.stays:
            s["mass"] = (float(self.wg[s["runs"][0]][0]) - float(self.wg[s["runs"][-1]][1])
                         + self.part.get(s["runs"][-1], 0.0))
            s["counts"] = all(self.counts[i] for i in s["runs"])
            for i in s["runs"]:
                if self.covered[i] is not None and not self.res:
                    v = self.raw[i][:, 1] + self.corr[i]
                    if v.min() < self.covered[i][0] or v.max() > self.covered[i][1]:
                        s["counts"] = False     # outside the certificate's range: the spare's error is unknown

    def through(self, c, kind):
        """Orifice and crystal surface in series (the surface passes vapour as evaporation coefficient x bed area)."""
        e = self.evap_cond[kind]
        return c if e is None else c * e / (c + e)

    def celsius(self, rid, v, lead):
        if not self.res:
            return v
        th = self.cellj["thermometer"]
        a, b, r0 = th["A_per_K"], th["B_per_K2"], th["R0_ohm"]
        if self.icecheck:
            r0 = self.icecheck[0] - lead
            a = ((self.icecheck[1] - lead) / r0 - 1 - 1e4 * b) / 100
        y = (v - lead) / r0 - 1
        return 2 * y / (a + np.sqrt(a * a + 4 * b * y))

    def cooled(self, T, lnp_f, dH_f, c, M):
        """Crystal temperature in a liner: the heat its sublimation takes comes through liner.conductance."""
        if self.G is None or self.rough:
            return T
        x = T.copy()
        for _ in range(60):
            q = c * np.exp(lnp_f(x)) * np.sqrt(M / (2 * math.pi * GAS * x)) * dH_f(x) / M
            q = q + self.hole * (x ** 4 - SHROUD_K ** 4)
            x_new = T - q / self.G
            done = np.max(np.abs(x_new - x)) < 1e-8
            x = x_new
            if done:
                break
        return x

    def ref_dH(self, T):
        """Sublimation enthalpy of the reference: -R d(ln p)/d(1/T) of the local quadratic."""
        x = 1.0 / np.asarray(T)
        k = np.clip(np.searchsorted(-self.ix, -x), 1, len(self.ix) - 2)
        x0, x1, x2 = self.ix[k - 1], self.ix[k], self.ix[k + 1]
        y0, y1, y2 = self.iy[k - 1], self.iy[k], self.iy[k + 1]
        d = (y0 * ((x - x1) + (x - x2)) / ((x0 - x1) * (x0 - x2)) + y1 * ((x - x0) + (x - x2)) / ((x1 - x0) * (x1 - x2))
             + y2 * ((x - x0) + (x - x1)) / ((x2 - x0) * (x2 - x1)))
        return -GAS * d

    def ref_lnp(self, T):
        x = 1.0 / np.asarray(T)
        k = np.clip(np.searchsorted(-self.ix, -x), 1, len(self.ix) - 2)
        x0, x1, x2 = self.ix[k - 1], self.ix[k], self.ix[k + 1]
        y0, y1, y2 = self.iy[k - 1], self.iy[k], self.iy[k + 1]
        return (y0 * (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2)) + y1 * (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
                + y2 * (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)))

    def losses(self, p):
        nf = self.nfoil
        h, lp, lnC = p[0], p[1], p[2:2 + nf]
        dT = p[2 + nf:2 + nf + self.nspare]
        lead = p[-1] if (self.two_wire and self.lead_known is None) else (self.lead_known or 0.0)
        out = []
        form2 = False
        for s in self.stays:
            if self.CL is not None and not self.rough:
                if self.kind[s["runs"][0]] == "reference":
                    form2 = False
                out.append(self.stay_lag(s, h, lp, lnC, dT, lead, form2))
                if (self.kind[s["runs"][0]] == "sample" and self.trans is not None and not form2
                        and self.peak >= self.trans[0]):
                    form2 = True
                continue
            m, t_prev_end = 0.0, s["from"]
            for k, rid in enumerate(s["runs"]):
                rows = self.raw[rid]
                Tc = self.celsius(rid, rows[:, 1], lead) + self.corr[rid]
                if self.spare[rid] is not None:
                    Tc = Tc - dT[self.spare[rid]]
                grid = np.arange(rows[0, 0], rows[-1, 0] + 1e-9, 2.0)
                T = np.interp(grid, rows[:, 0], Tc) + 273.15
                c = self.through(math.exp(lnC[self.foil[rid]]), self.kind[rid])
                M = self.Mr if self.kind[rid] == "reference" else self.Ms
                if self.kind[rid] == "reference":
                    form2 = False
                    Tx = self.cooled(T, self.ref_lnp, self.ref_dH, c, M)
                    lnp = self.ref_lnp(Tx)
                else:
                    def lnp_f(u, two):
                        y = form1(u, h, lp, self.dc)
                        return y + form2_shift(u, *self.trans) if two else y

                    def dH_f(u, two):
                        y = h + self.dc * (u - TREF)
                        if two:
                            Tt, ht, dcp = self.trans
                            y = y - (ht + dcp * (u - Tt))
                        return y
                    Tx = self.cooled(T, lambda u: lnp_f(u, form2), lambda u: dH_f(u, form2), c, M)
                    lnp = lnp_f(Tx, form2)
                    if self.trans is not None and not form2 and Tx.max() >= self.trans[0]:
                        first = int(np.argmax(Tx >= self.trans[0]))
                        tail = self.cooled(T[first:], lambda u: lnp_f(u, True), lambda u: dH_f(u, True), c, M)
                        Tx = np.concatenate([Tx[:first], tail])
                        lnp = np.concatenate([lnp[:first], lnp_f(tail, True)])
                        form2 = True
                rate = c * np.exp(lnp) * np.sqrt(M / (2 * math.pi * GAS * Tx))
                # wait before the program, at its first reading (before any change of form in this program)
                w = (self.start[rid] - t_prev_end).total_seconds()
                r0 = rate[0]
                m += r0 * w + simpson(rate, x=grid)
                t_prev_end = self.start[rid].timestamp() + rows[-1, 0]
                t_prev_end = datetime.fromtimestamp(t_prev_end)
                last_rate = rate[-1]
            m += last_rate * (s["to"] - t_prev_end).total_seconds()
            out.append(m)
        return np.array(out)

    def stay_lag(self, s, h, lp, lnC, dT, lead, two=False):
        """Loss (kg) of a stay in a liner cell with heat capacity: the whole stay on a 2 s grid (waits at the first
        and last reading of the neighbouring program), the crystal temperature from
        (C_liner + n Cp) dT/dt = G (T_block - T) - q by the trapezoidal (Crank-Nicolson) rule, q re-evaluated until
        the trajectory settles."""
        grid_T, grid_c = [], []
        t_prev_end = s["from"]
        for rid in s["runs"]:
            rows = self.raw[rid]
            Tc = self.celsius(rid, rows[:, 1], lead) + self.corr[rid]
            if self.spare[rid] is not None:
                Tc = Tc - dT[self.spare[rid]]
            c = self.through(math.exp(lnC[self.foil[rid]]), self.kind[rid])
            w = (self.start[rid] - t_prev_end).total_seconds()
            grid_T.append(np.full(int(round(w / 2.0)), Tc[0] + 273.15))
            g = np.arange(rows[0, 0], rows[-1, 0] + 1e-9, 2.0)
            grid_T.append(np.interp(g, rows[:, 0], Tc)[1:] + 273.15)
            grid_c.append(np.full(len(grid_T[-2]) + len(grid_T[-1]), c))
            t_prev_end = datetime.fromtimestamp(self.start[rid].timestamp() + rows[-1, 0])
            last = Tc[-1] + 273.15
        n_end = int(round((s["to"] - t_prev_end).total_seconds() / 2.0))
        grid_T.append(np.full(n_end, last))
        grid_c.append(np.full(n_end, c))
        Tb, cv = np.concatenate(grid_T), np.concatenate(grid_c)
        if self.kind[s["runs"][0]] == "reference":
            lnp_f, dH_f, M, cp = self.ref_lnp, self.ref_dH, self.Mr, self.cp_ref
        else:
            if two:
                Tt, ht, dcp = self.trans
                lnp_f = lambda u: form1(u, h, lp, self.dc) + form2_shift(u, Tt, ht, dcp)
                dH_f = lambda u: h + self.dc * (u - TREF) - (ht + dcp * (u - Tt))
            else:
                lnp_f = lambda u: form1(u, h, lp, self.dc)
                dH_f = lambda u: h + self.dc * (u - TREF)
            M, cp = self.Ms, self.cp_smp
        charge = (float(self.wg[s["runs"][0]][0]) - self.empty) * 1e-6
        T = Tb.copy()
        lost = np.zeros(len(Tb))     # kg gone at each grid point: the crystals' heat capacity falls with it
        for _ in range(14):
            flux = cv * np.exp(lnp_f(T)) * np.sqrt(M / (2 * math.pi * GAS * T))
            q = flux * dH_f(T) / M + self.hole * (T ** 4 - SHROUD_K ** 4)
            u = Tb - q / self.G
            new = np.empty(len(Tb))
            prev_T, prev_u = Tb[0], u[0]
            for i0 in range(0, len(Tb), 120):       # 4 min blocks of the 2 s grid
                i1 = min(len(Tb), i0 + 120)
                tau = (self.CL + (charge - lost[(i0 + i1) // 2]) / M * cp) / self.G
                k = 2.0 / (2 * tau)
                b, a = [k / (1 + k), k / (1 + k)], [1.0, -(1 - k) / (1 + k)]
                new[i0:i1] = lfilter(b, a, u[i0:i1], zi=[(-a[1]) * prev_T + b[1] * prev_u])[0]
                prev_T, prev_u = new[i1 - 1], u[i1 - 1]
            lost = np.concatenate([[0.0], np.cumsum(0.5 * (flux[1:] + flux[:-1]) * 2.0)])
            done = np.max(np.abs(new - T)) < 1e-7
            T = new
            if done:
                break
        self.peak = float(T.max())
        rate = cv * np.exp(lnp_f(T)) * np.sqrt(M / (2 * math.pi * GAS * T))
        return float(simpson(rate, dx=2.0))

    def solve(self):
        nf, ns = self.nfoil, self.nspare
        free_lead = self.two_wire and self.lead_known is None
        p = np.array([90e3, math.log(1e-2)] + [math.log(5e-7)] * nf + [0.0] * ns + [0.0] * free_lead)
        keep = np.array([s["counts"] for s in self.stays])
        obs = np.array([s["mass"] for s in self.stays]) * 1e-6
        # starting conductance from the reference stays, then a coarse search for the sample
        refs = [j for j, s in enumerate(self.stays) if keep[j] and self.kind[s["runs"][0]] == "reference"]
        self.rough = True
        q = p.copy()
        q[2:2 + nf] = 0.0
        pred = self.losses(q)
        for f in range(nf):
            js = [j for j in refs if self.foil[self.stays[j]["runs"][0]] == f]
            p[2 + f] = math.log(np.mean(obs[js] / pred[js])) if js else np.nan
        good = p[2:2 + nf][np.isfinite(p[2:2 + nf])]
        p[2:2 + nf] = np.where(np.isfinite(p[2:2 + nf]), p[2:2 + nf], good.mean())
        best = None
        for h in np.arange(50e3, 160e3, 3e3):
            for lp in np.arange(-14.0, 4.0, 0.5):
                q = p.copy()
                q[0], q[1] = h, lp
                with np.errstate(over="ignore", invalid="ignore"):
                    r = np.log(self.losses(q)[keep] / obs[keep])
                if np.all(np.isfinite(r)) and (best is None or r @ r < best[0]):
                    best = (r @ r, h, lp)
        p[0], p[1] = best[1], best[2]
        self.rough = False
        scale = np.array([1e3, 1.0] + [1.0] * nf + [0.1] * ns + [0.05] * free_lead)

        def resid(p):
            mg = obs[keep] * 1e6
            return np.log(self.losses(p)[keep] / obs[keep]) * mg / np.sqrt(0.003 ** 2 + (3e-4 * mg) ** 2)
        lam = 1e-3
        r = resid(p)
        for _ in range(200):
            J = np.empty((len(r), len(p)))
            for k in range(len(p)):
                dp = np.zeros_like(p)
                dp[k] = 1e-4 * scale[k]
                J[:, k] = (resid(p + dp) - resid(p - dp)) / (2 * dp[k])
            A, g = J.T @ J, J.T @ r
            while True:
                step = np.linalg.solve(A + lam * np.diag(np.diag(A)), -g)
                r_new = resid(p + step)
                if r_new @ r_new < r @ r:
                    p, r, lam = p + step, r_new, lam / 3
                    break
                lam *= 4
                if lam > 1e12:
                    break
            if np.all(np.abs(step) < 1e-9 * scale) or lam > 1e12:
                break
        return {"dH_sub_298_kJ_mol": p[0] / 1000, "p_298_Pa": math.exp(p[1]), "chi2_per_stay": float(r @ r / len(r))}


def main():
    Path(sys.argv[2]).write_text(json.dumps(Series(sys.argv[1]).solve(), indent=1) + "\n")


if __name__ == "__main__":
    main()
