"""Sublimation enthalpy and vapour pressure at 298.15 K from one Knudsen effusion series.

    python3 analyse_series.py <series folder> <answer.json>

1. Read the sheets, the weighings, every run's temperature record and the operator's log. Pt100 resistances are
   converted with the stated quadratic after the leads' resistance is taken off (a 2-wire Pt100 whose leads were
   not measured gets that resistance as a fitted unknown).
2. From the log: which runs share a stay in the chamber, when each stay reached vacuum, when each program started
   and when the chamber was vented; which runs' losses are not used; where an orifice foil was changed; where a
   thermometer was replaced (with a certificate: a known offset; without one: an offset to be fitted). A
   certificate covers only its stated range: a stay in which a run's temperatures leave that range is not used,
   since nothing is known of the spare's error there.
3. Predict the mass each stay loses: the effusion rate integrated over the whole stay - the wait at the program's
   starting temperature, the program itself (its record), any wait between two programs and the wait at the last
   program's end temperature until venting. The sample's vapour pressure follows the stated heat capacities; a
   sample with a stated transition changes form where its temperature first reaches the transition temperature,
   in whatever run that happens, and keeps the new form. In a liner cell with the thermometer in the block the
   crystal temperature is the liner's: (C_liner + n Cp) dT/dt = G (T_block - T) - q_sublimation - q_radiation,
   q_radiation = sigma A_orifice (T^4 - T_shroud^4) because the orifice of the closed liner looks at the 77 K
   shroud and a small hole in a closed cup radiates as a black body.
4. Fit everything at once by weighted least squares on the weighed losses of the stays that count: dH298 and
   ln p298 of the sample as loaded, the effusion conductance of every orifice foil (a foil without reference runs
   is fixed by the sample runs made with it), the offset of every spare thermometer without a certificate and an
   unmeasured lead resistance.
"""
import csv
import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
from scipy.signal import lfilter

R, T0 = 8.314462618, 298.15
SIGMA_MG = 0.003            # scatter of a weighed loss (two weighings)
SIGMA_REL = 3e-4            # and the part that grows with the loss (thermometer noise in the waits)
STEFAN = 5.670374419e-8     # W/(m^2 K^4)
T_SHROUD = 77.0             # K, liquid-nitrogen shroud facing the orifice


def end_temperature(t, T, first):
    """The block temperature at the start (first=True) or end of a program: a straight line through the
    record's first or last 7 lines (all on the ramp), taken at the end."""
    i = slice(0, 7) if first else slice(-7, None)
    k = np.polyfit(t[i], T[i], 1)
    return float(np.polyval(k, t[0] if first else t[-1]))


def ln_p_a(T, dH298, lnp298, dCp):
    return lnp298 - dH298 / R * (1 / T - 1 / T0) + dCp / R * ((T0 / T - 1) + np.log(T / T0))


def ln_p_b_minus_a(T, Tt, h, c):
    return -((h - c * Tt) * (1 / Tt - 1 / T) + c * np.log(T / Tt)) / R


def named(text, ids):
    """Run ids a remark names: '<a> to <b>' and '<a> and <b>' groups, else the single ids."""
    out = []
    for a, sep, b in re.findall(r"R(\d+)\s+(to|and)\s+R(\d+)", text):
        out += [f"R{i:02d}" for i in range(int(a), int(b) + 1)] if sep == "to" else [f"R{int(a):02d}", f"R{int(b):02d}"]
    if not out:
        out = re.findall(r"R\d+", text)
    return [x for x in out if x in ids]


def reference_law(path):
    """ln p of the reference and its sublimation enthalpy (Clausius-Clapeyron slope), from the certificate table as
    a cubic spline of ln p against 1/T."""
    d = np.genfromtxt(path, delimiter=",", names=True)
    x = 1.0 / (d["temperature_C"] + 273.15)
    order = np.argsort(x)
    spline = CubicSpline(x[order], np.log(d["p_Pa"])[order])
    slope = spline.derivative()
    return (lambda T: spline(1.0 / np.asarray(T))), (lambda T: -R * slope(1.0 / np.asarray(T)))


def read_series(folder):
    folder = Path(folder)
    cell = json.loads((folder / "cell.json").read_text())
    comp = json.loads((folder / "compounds.json").read_text())
    with open(folder / "runs.csv", newline="") as fh:
        runs = list(csv.DictReader(fh))
    with open(folder / "weighings.csv", newline="") as fh:
        weigh = {r["run"]: (r["before_mg"].strip(), r["after_mg"].strip()) for r in csv.DictReader(fh)}
    th = cell["thermometer"]
    for r in runs:
        d = np.genfromtxt(folder / r["file"], delimiter=",", names=True)
        r["t"], r["read"] = d["time_s"], d[th["column"]]
        r["use"], r["foil"], r["known_offset"], r["spare"], r["cert"] = True, 0, 0.0, None, None
    ids = [r["run"] for r in runs]
    by = {r["run"]: r for r in runs}
    stays, cur = [], None
    foil, known, spare, n_spare, cert = 0, 0.0, None, 0, None
    extra = {}                  # weighings made with a part still on the cell: its mass
    check = None                # a Pt100's resistances read at 0 and 100 C (leads included)
    for line in (folder / "run_log.txt").read_text().splitlines():
        stamp, text = line.split(" | ", 1)
        when = datetime.strptime(stamp.strip(), "%d.%m.%Y %H:%M:%S")
        m = re.match(r"(R\d+) (under vacuum|started|vented|weighed)(?=\.$| \()", text)
        if m:
            rid, what = m.groups()
            if what == "under vacuum":
                cur = {"vacuum": when, "runs": []}
                stays.append(cur)
            elif what == "started":
                by[rid].update(start=when, foil=foil, known_offset=known, spare=spare, cert=cert)
                cur["runs"].append(rid)
            elif what == "vented":
                cur["vented"] = when
            continue
        low = text.lower()
        # a certificate gives the error ("reads x K low/high": true = reading +/- x) or the correction ("correction
        # c K": true = reading + c), each over the range t1 to t2 it covers
        m = re.search(r"certificate: reads ([\d.]+) k (low|high) between (-?[\d.]+) and (-?[\d.]+) c", low)
        mc = re.search(r"certificate: correction ([+-]?[\d.]+) k between (-?[\d.]+) and (-?[\d.]+) c", low)
        k = re.search(r"fitted for (R\d+) onwards", text)
        if k and "spare" in low:
            first = ids.index(k.group(1))
            if m or mc:
                if m:
                    delta = float(m.group(1)) * (1 if m.group(2) == "low" else -1)
                    rng = (float(m.group(3)), float(m.group(4)))
                else:
                    delta, rng = float(mc.group(1)), (float(mc.group(2)), float(mc.group(3)))
                for i in ids[first:]:
                    by[i]["known_offset"], by[i]["spare"], by[i]["cert"] = delta, None, rng
                known, spare, cert = delta, None, rng
            elif "no certificate" in low:
                for i in ids[first:]:
                    by[i]["known_offset"], by[i]["spare"], by[i]["cert"] = 0.0, n_spare, None
                known, spare, n_spare, cert = 0.0, n_spare, n_spare + 1, None
        if "not for evaluation" in low:
            for i in named(text.split("not for evaluation")[0], ids):
                by[i]["use"] = False
        m = re.search(r"repeated as R\d+", text)
        if m:
            for i in named(text[:m.start()], ids):
                by[i]["use"] = False
        if "new foil of the same nominal size fitted" in low:
            foil += 1
        m = re.search(r"read ([\d.]+) ohm at the ice point and ([\d.]+) ohm at 100 c", low)
        if m:
            check = (float(m.group(1)), float(m.group(2)))
        m = re.match(r"(R\d+) weighed with .*\(.*?([\d.]+) mg on its own\)", text)
        if m:
            extra[m.group(1)] = float(m.group(2))
    for s in stays:
        first, last = s["runs"][0], s["runs"][-1]
        s["before"] = float(weigh[first][0])
        s["extra_mg"] = extra.get(last, 0.0)
        s["after"] = float(weigh[last][1]) - s["extra_mg"]
        s["use"] = all(by[i]["use"] for i in s["runs"])
    cell["thermometer"]["check_ohm"] = check
    return cell, comp, runs, by, stays


def pt100_celsius(th, Rv, lead):
    """Temperature of a Pt100 from the resistance read; a check at 0 and 100 C gives the element's own R0 and A."""
    A, B, R0 = th["A_per_K"], th["B_per_K2"], th["R0_ohm"]
    if th.get("check_ohm"):
        r1, r2 = th["check_ohm"]
        R0 = r1 - lead
        A = ((r2 - lead) / R0 - 1 - 1e4 * B) / 100
    ratio = (Rv - lead) / R0
    return (-A + np.sqrt(A * A - 4 * B * (1 - ratio))) / (2 * B)


class Model:
    def __init__(self, folder):
        self.cell, comp, runs, self.by, self.stays = read_series(folder)
        th = self.th = self.cell["thermometer"]
        self.pt100 = th["column"] == "resistance_ohm"
        # a 2-wire Pt100 reads its leads too; when their resistance was not measured it is fitted
        self.lead_free = self.pt100 and th.get("wires", 4) == 2 and "lead_resistance_ohm" not in th
        self.lead = th.get("lead_resistance_ohm", 0.0) if self.pt100 and th.get("wires", 4) == 2 else 0.0
        ref, smp = comp["reference"], comp["sample"]
        self.Mr, self.Ms = ref["molar_mass_g_mol"] / 1000, smp["molar_mass_g_mol"] / 1000
        self.dCp = smp["Cp_gas_J_K_mol"] - smp["Cp_crystal_J_K_mol"]
        self.ref_lnp, self.ref_dH = reference_law(Path(folder) / ref["vapour_pressure_table"])
        # a crystal in a liner is cooled by the heat its sublimation takes, which reaches it through the liner
        self.G = self.cell.get("liner", {}).get("conductance_W_per_K")
        # with a heat capacity the liner and its crystals lag the block: (C_liner + n Cp) dT/dt = G (T_block - T) - q
        self.Cl = self.cell.get("liner", {}).get("heat_capacity_J_per_K")
        if th.get("location") == "liner":
            # the controller holds the thermometer in the liner, and so the crystals, on the program
            self.G, self.Cl = None, None
        self.empty_mg = self.cell.get("empty_mass_mg")
        # the liner's orifice radiates to the shroud as a black body of the orifice's area
        self.rad = STEFAN * np.pi * (self.cell["orifice_diameter_mm"] * 5e-4) ** 2
        self.T_shroud = T_SHROUD
        self.form_from_block = False    # misreading switch: decide the form on the block temperature
        self.Cp_ref = ref.get("Cp_crystal_J_K_mol")
        # crystals with an evaporation coefficient below one pass vapour like a second orifice (coefficient times
        # the bed area) in series with the cell's orifice
        bed = self.cell.get("bed_area_mm2", 0.0) * 1e-6
        self.surf = {k: (d["evaporation_coefficient"] * bed if "evaporation_coefficient" in d else None)
                     for k, d in (("reference", ref), ("sample", smp))}
        self.Cp_smp = smp["Cp_crystal_J_K_mol"]
        self.coarse = False         # True only while the starting values are searched
        trs = smp.get("transition")
        self.trs = None if trs is None else (trs["temperature_C"] + 273.15, trs["dH_kJ_mol"] * 1e3,
                                             trs["Cp_crystal_J_K_mol"] - smp["Cp_crystal_J_K_mol"])
        self.layout()

    def layout(self):
        by, stays = self.by, self.stays
        for s in stays:
            for rid in s["runs"]:
                r = by[rid]
                if r["cert"] is not None:
                    t = r["read"] + r["known_offset"]
                    if t.min() < r["cert"][0] or t.max() > r["cert"][1]:
                        s["use"] = False
        self.use = np.array([s["use"] for s in stays])
        self.foils = sorted({by[i]["foil"] for s in stays for i in s["runs"]})
        self.n_spare = 1 + max([by[i]["spare"] for s in stays for i in s["runs"] if by[i]["spare"] is not None],
                               default=-1)
        # every stay as a list of segments: (kind, run id, record times or None, temperatures read, wait in s)
        self.segments = []
        for s in stays:
            seg, end_prev = [], None
            for k, rid in enumerate(s["runs"]):
                r = by[rid]
                wait = ((r["start"] - s["vacuum"]).total_seconds() if k == 0
                        else (r["start"] - end_prev).total_seconds())
                seg.append(("wait", rid, None, np.array([end_temperature(r["t"], r["read"], True)]), wait))
                seg.append(("program", rid, r["t"], r["read"], None))
                end_prev = r["start"] + timedelta(seconds=float(r["t"][-1]))
            last = by[s["runs"][-1]]
            seg.append(("wait", s["runs"][-1], None, np.array([end_temperature(last["t"], last["read"], False)]),
                        (s["vented"] - end_prev).total_seconds()))
            self.segments.append(seg)
        self.obs = np.array([s["before"] - s["after"] for s in stays])

    def escape(self, C, compound):
        a = self.surf[compound]
        return C if a is None else 1.0 / (1.0 / C + 1.0 / a)

    def kelvin(self, rid, read, x):
        """Cell temperature (K) of readings of run rid for parameters x."""
        r = self.by[rid]
        if self.pt100:
            read = pt100_celsius(self.th, read, x[-1] if self.lead_free else self.lead)
        T = read + r["known_offset"] + 273.15
        if r["spare"] is not None:
            T = T - x[2 + len(self.foils) + r["spare"]]
        return T

    def lnp_sample(self, T, x, form_b):
        la = ln_p_a(T, x[0] * 1e3, x[1], self.dCp)
        return la + ln_p_b_minus_a(T, *self.trs) if form_b else la

    def dH_sample(self, T, x, form_b):
        h = x[0] * 1e3 + self.dCp * (T - T0)
        if form_b:
            Tt, ht, c = self.trs
            h = h - (ht + c * (T - Tt))
        return h

    def crystal(self, T, lnp_fun, dH_fun, C, M):
        """Crystal temperature: the block temperature less the heat flow its sublimation takes over the liner's
        conductance (solved by fixed-point iteration); the block temperature itself without a liner."""
        if self.G is None or self.coarse:
            return T
        Tc = T.copy()
        for _ in range(30):
            q = C * np.exp(lnp_fun(Tc)) * np.sqrt(M / (2 * np.pi * R * Tc)) * dH_fun(Tc) / M
            q = q + self.rad * (Tc ** 4 - self.T_shroud ** 4)
            new = T - q / self.G
            if np.max(np.abs(new - Tc)) < 1e-7:
                return new
            Tc = new
        return Tc

    def stay_dynamic(self, seg, x, lnC, charge_kg, form_b=False):
        """Loss (mg) of one stay of a liner cell whose liner and crystals have heat capacity: the crystal temperature
        follows (C_liner + n Cp) dT/dt = G (T_block - T) - rate dH / M - sigma A (T^4 - T_shroud^4) from the block
        temperature at 'under vacuum', solved exactly for a block temperature held over each 10 s step (first-order
        filter, the crystals' heat capacity following the mass still in the cell), with the sublimation heat and
        the orifice's radiation iterated to convergence. Keeps the highest
        crystal (and block) temperature of the stay for the sample's form."""
        Tb, Cv = [], []
        for kind, rid, t, read, wait in seg:
            r = self.by[rid]
            T = self.kelvin(rid, read, x)
            C = self.escape(np.exp(lnC[self.foils.index(r["foil"])]) * 1e-7, r["compound"])
            if kind == "wait":
                n = int(round(wait / 10.0))
                Tb.append(np.full(n, T[0]))
                Cv.append(np.full(n, C))
            else:
                Tb.append(0.5 * (T[1:] + T[:-1]))
                Cv.append(np.full(len(T) - 1, C))
            is_ref = r["compound"] == "reference"
        Tb, Cv = np.concatenate(Tb), np.concatenate(Cv)
        if is_ref:
            lnp_f, dH_f, M, Cp = self.ref_lnp, self.ref_dH, self.Mr, self.Cp_ref
        else:
            lnp_f, dH_f = (lambda u: self.lnp_sample(u, x, form_b)), (lambda u: self.dH_sample(u, x, form_b))
            M, Cp = self.Ms, self.Cp_smp
        T0 = Tb[0]
        Tend = Tb.copy()
        gone = np.zeros(len(Tb))         # kg lost before each step (the crystals' heat capacity shrinks with it)
        for _ in range(10):
            Tstart = np.concatenate([[T0], Tend[:-1]])
            Tm = 0.5 * (Tstart + Tend)
            flow = Cv * np.exp(lnp_f(Tm)) * np.sqrt(M / (2 * np.pi * R * Tm))
            q = flow * dH_f(Tm) / M + self.rad * (Tm ** 4 - self.T_shroud ** 4)
            u = Tb - q / self.G
            new, state = np.empty(len(Tb)), T0
            for c0 in range(0, len(Tb), 30):        # 5 min pieces, each with its own heat capacity
                c1 = min(len(Tb), c0 + 30)
                a = np.exp(-10.0 * self.G / (self.Cl + (charge_kg - gone[(c0 + c1) // 2]) / M * Cp))
                new[c0:c1] = lfilter([1 - a], [1, -a], u[c0:c1], zi=[a * state])[0]
                state = new[c1 - 1]
            gone = np.concatenate([[0.0], np.cumsum(flow * 10.0)[:-1]])
            done = np.max(np.abs(new - Tend)) < 1e-7
            Tend = new
            if done:
                break
        Tstart = np.concatenate([[T0], Tend[:-1]])
        Tm = 0.5 * (Tstart + Tend)
        self.top = float(Tb.max()) if self.form_from_block else float(Tend.max())
        return float(np.sum(Cv * np.exp(lnp_f(Tm)) * np.sqrt(M / (2 * np.pi * R * Tm)))) * 10.0 * 1e6

    def predict(self, x):
        lnC = x[2:2 + len(self.foils)]
        out = np.zeros(len(self.stays))
        form_b = False
        for j, seg in enumerate(self.segments):
            if self.Cl is not None and not self.coarse:
                is_ref = self.by[self.stays[j]["runs"][0]]["compound"] == "reference"
                if is_ref:
                    form_b = False          # the cell was emptied of the sample
                out[j] = self.stay_dynamic(seg, x, lnC, (self.stays[j]["before"] - self.empty_mg) * 1e-6,
                                           form_b and not is_ref)
                if not is_ref and self.trs is not None and not form_b and self.top >= self.trs[0]:
                    form_b = True           # the crystal reached the transition in this stay
                continue
            total = 0.0
            for kind, rid, t, read, wait in seg:
                r = self.by[rid]
                T = self.kelvin(rid, read, x)
                C = self.escape(np.exp(lnC[self.foils.index(r["foil"])]) * 1e-7, r["compound"])
                if r["compound"] == "reference":
                    form_b = False          # the cell was emptied of the sample
                    Tc = self.crystal(T, self.ref_lnp, self.ref_dH, C, self.Mr)
                    f = np.exp(self.ref_lnp(Tc)) * np.sqrt(self.Mr / (2 * np.pi * R * Tc))
                else:
                    Tc = self.crystal(T, lambda u: self.lnp_sample(u, x, form_b),
                                      lambda u: self.dH_sample(u, x, form_b), C, self.Ms)
                    lnp = self.lnp_sample(Tc, x, form_b)
                    if self.trs is not None and not form_b and kind == "program" and Tc.max() >= self.trs[0]:
                        # the crystal changes form where it first reaches the transition temperature
                        i = int(np.argmax(Tc >= self.trs[0]))
                        Tb = self.crystal(T[i:], lambda u: self.lnp_sample(u, x, True),
                                          lambda u: self.dH_sample(u, x, True), C, self.Ms)
                        Tc = np.concatenate([Tc[:i], Tb])
                        lnp = np.concatenate([lnp[:i], self.lnp_sample(Tb, x, True)])
                        form_b = True
                    f = np.exp(lnp) * np.sqrt(self.Ms / (2 * np.pi * R * Tc))
                total += C * (f[0] * wait if kind == "wait" else np.trapezoid(f, t)) * 1e6
            out[j] = total
        return out

    def residuals(self, x):
        sig = np.sqrt(SIGMA_MG ** 2 + (SIGMA_REL * self.obs) ** 2)
        return ((self.predict(x) - self.obs) / sig)[self.use]


def fit(m, fixed=None):
    """Least-squares fit of the model; fixed: {parameter index: value} held fixed."""
    nf = len(m.foils)
    x0 = np.array([90.0, np.log(1e-2)] + [np.log(5.0)] * nf + [0.0] * m.n_spare + [0.0] * m.lead_free)
    # conductances first, from the reference stays alone (the sample parameters do not enter them)
    ref_stays = [j for j, s in enumerate(m.stays)
                 if m.use[j] and all(m.by[i]["compound"] == "reference" for i in s["runs"])]
    for k, f in enumerate(m.foils):
        js = [j for j in ref_stays if m.by[m.stays[j]["runs"][0]]["foil"] == f]
        x = x0.copy()
        x[2 + k] = 0.0
        x0[2 + k] = float(np.log(np.mean(m.obs[js] / m.predict(x)[js]))) if js else np.nan
    x0[2:2 + nf] = np.where(np.isnan(x0[2:2 + nf]), np.nanmean(x0[2:2 + nf]), x0[2:2 + nf])
    fixed = fixed or {}
    for i, v in fixed.items():
        x0[i] = v
    best = None
    m.coarse = True
    for dH in np.arange(50.0, 160.0, 2.5):
        for lnp in np.arange(-14.0, 4.0, 0.5):
            x = x0.copy()
            x[0], x[1] = dH, lnp
            with np.errstate(over="ignore", invalid="ignore"):
                pr = m.predict(x)
            if not np.all(np.isfinite(pr)) or np.any(pr <= 0):
                continue
            c = float(np.sum(np.log(pr / np.abs(m.obs))[m.use] ** 2))
            if best is None or c < best[0]:
                best = (c, dH, lnp)
    x0[0], x0[1] = best[1], best[2]
    m.coarse = False
    free = [i for i in range(len(x0)) if i not in fixed]

    def res(z):
        x = x0.copy()
        x[free] = z
        return m.residuals(x)
    f = least_squares(res, x0[free], x_scale="jac", xtol=1e-12, ftol=1e-12)
    x = x0.copy()
    x[free] = f.x
    return x, float(2 * f.cost / int(m.use.sum()))


def solve(folder):
    m = Model(folder)
    x, chi2 = fit(m)
    nf = len(m.foils)
    return {"dH_sub_298_kJ_mol": float(x[0]), "p_298_Pa": float(np.exp(x[1])), "chi2_per_stay": chi2,
            "conductance_m2": {str(m.foils[k]): float(np.exp(x[2 + k]) * 1e-7) for k in range(nf)},
            "spare_reads_high_K": [float(v) for v in x[2 + nf:2 + nf + m.n_spare]],
            "lead_resistance_ohm": float(x[-1]) if m.lead_free else m.lead}


def main():
    Path(sys.argv[2]).write_text(json.dumps(solve(sys.argv[1]), indent=1) + "\n")


if __name__ == "__main__":
    main()
