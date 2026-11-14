"""Reduce one titration-calorimetry run to its dissociation constant and binding enthalpy.

    python3 itc_reduce.py <run folder> <answer.json>

Steps:
1. Every record (each syringe fill, and the blank) is fitted as a whole: baseline (piecewise linear, knots every
   600 s) plus one peak per injection, each peak the response of a first-order instrument lag to a release spread
   over the injection and mixed in with a short lag of its own, and one peak per calibration pulse with the
   instrument lag only. Areas and baseline enter linearly; the two lags of each record are found by least squares
   over that record, so each record keeps the lag of its own feedback mode.
2. The power scale belongs to the instrument on the day: it is the mean of pulse area over stated energy over every
   pulse fired in any record of the day.
3. Titrant concentrations: as measured, or from the weighing record with the molar mass of the solid actually
   weighed (parent formula plus counter-ions and water of the named salt or hydrate), times its purity.
4. Delivered volumes are the scheduled ones unless an operator remark gives the delivered volumes of injections.
5. The blank gives the titrant's heat of dilution per mole once its own dimer dissociation (if the compound
   self-associates) is accounted for; the titration model then carries the same dilution heat, its dimer
   dissociation computed with the binding, and the binding itself.
6. Binding isotherm: overfilled cell of fixed volume, each injection expelling its own volume of the pre-injection
   solution; one site per active macromolecule; only the ligand monomer binds; a factor on the cell species'
   concentration is fitted with Kd and dH. The first injection of every fill is left out; its volume is kept.
"""
import json
import os
import re
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import brentq, least_squares  # noqa: E402

UCAL_PER_KJ_X_PMOL = 1e-3 / 4.184
ATOM = {"H": 1.008, "Cl": 35.45, "Br": 79.904, "Na": 22.990, "K": 39.098, "O": 15.999}
H2O = 2 * ATOM["H"] + ATOM["O"]
COUNT = {"": 1, "mono": 1, "di": 2, "tri": 3, "tetra": 4, "hemi": 0.5, "sesqui": 1.5}


def load_record(path):
    return np.genfromtxt(path, delimiter=",", names=True)


# ---------------------------------------------------------------- records
def ramp_integral(x, ta, tb):
    x = np.maximum(x, 0.0)
    if tb <= 0.0:
        return x - ta + ta * np.exp(-x / ta)
    if abs(ta - tb) < 1e-6 * ta:
        tb = ta * (1 + 1e-6)
    return x - ta - tb + (ta * ta * np.exp(-x / ta) - tb * tb * np.exp(-x / tb)) / (ta - tb)


def response_column(t, period, start, dur, ta, tb):
    S = ramp_integral
    a = t - start
    return (S(a, ta, tb) - S(a - period, ta, tb) - S(a - dur, ta, tb) + S(a - dur - period, ta, tb)) / (period * dur)


class RecordModel:
    def __init__(self, t, y, period, peaks, pulses, knot_s=600.0):
        self.t, self.y, self.period = t, y, period
        self.peaks, self.pulses = peaks, pulses
        n_k = max(2, int(np.ceil((t[-1] - t[0]) / knot_s)) + 1)
        knots = np.linspace(t[0], t[-1], n_k)
        eye = np.eye(n_k)
        self.base = np.array([np.interp(t, knots, eye[k]) for k in range(n_k)]).T

    def design(self, ta, tm):
        cols = [response_column(self.t, self.period, s, d, ta, tm) for s, d in self.peaks]
        cols += [response_column(self.t, self.period, s, d, ta, 0.0) for s, d in self.pulses]
        return np.column_stack(cols + [self.base])

    def solve(self, ta, tm):
        A = self.design(ta, tm)
        coef, *_ = np.linalg.lstsq(A, self.y, rcond=None)
        return coef, self.y - A @ coef

    def fit(self):
        best = None
        for ta in (6.0, 10.0, 16.0, 25.0, 40.0, 60.0, 90.0):
            r = self.solve(ta, 4.0)[1]
            if best is None or r @ r < best[0]:
                best = (r @ r, ta)
        sol = least_squares(lambda p: self.solve(*np.exp(p))[1], np.log([best[1], 4.0]), x_scale=[0.2, 0.3],
                            bounds=([np.log(2.0), np.log(0.3)], [np.log(400.0), np.log(60.0)]))
        ta, tm = np.exp(sol.x)
        coef, _ = self.solve(ta, tm)
        n = len(self.peaks)
        return dict(tau_inst=ta, tau_mix=tm, areas=coef[:n], pulse_areas=coef[n:n + len(self.pulses)])


def decompose_record(t, y, period, injections, pulses):
    return RecordModel(t, y, period, [(j["start_s"], j["duration_s"]) for j in injections],
                    [(p["start_s"], p["duration_s"]) for p in pulses]).fit()


# ---------------------------------------------------------------- solutions
def weighed_form_mw(parent_mw, solid_form):
    """Molar mass of the weighed solid: parent plus the counter-ions and water its name gives."""
    name = solid_form.lower()
    mw = parent_mw
    m = re.search(r"\b(mono|di|tri|tetra|)hydro(chloride|bromide)\b", name)
    if m:
        mw += COUNT[m.group(1)] * (ATOM["H"] + (ATOM["Cl"] if m.group(2) == "chloride" else ATOM["Br"]))
    m = re.search(r"\b(mono|di|tri|tetra|)(sodium|potassium) salt\b", name)
    if m:
        mw += COUNT[m.group(1)] * ((ATOM["Na"] if m.group(2) == "sodium" else ATOM["K"]) - ATOM["H"])
    m = re.search(r"\b(hemi|mono|sesqui|di|tri|tetra)hydrate\b", name)
    if m:
        mw += COUNT[m.group(1)] * H2O
    return mw


def remark_checks(sheet):
    """Later checks the operator reported that replace sheet values: a stock re-measured, the heater re-measured."""
    stocks, pulse = {}, 1.0
    for text in sheet.get("remarks", []):
        for sid, c in re.findall(r"\b(S\d+) re-measured[^:]*:\s*([\d.]+)\s*uM", text):
            stocks[sid] = float(c)
        m = re.search(r"pulse put ([\d.]+) % (more|less) energy", text)
        if m:
            pulse = 1.0 + (1 if m.group(2) == "more" else -1) * float(m.group(1)) / 100.0
    return stocks, pulse


def stock_um(sheet, sid):
    sol = sheet["solutions"][sid]
    remeasured = remark_checks(sheet)[0]
    if sid in remeasured:
        return remeasured[sid]
    if "conc_um" in sol:
        return float(sol["conc_um"])
    p = sol["prep"]
    mw = weighed_form_mw(float(sheet["compounds"][sol["compound"]]["mw_g_per_mol"]), p["solid_form"])
    return p["mass_mg"] * p["purity_pct"] / 100.0 * 1e6 / (mw * p["flask_ml"])


def volumes_from_notes(sheet):
    """{(fill, injection number): delivered volume} from remarks such as
    'injections 2, 3 and 4 of fill 2 delivered 5.9, 6.6 and 8.1 uL'."""
    out = {}
    pat = re.compile(r"injections?\s+([\d,\sand]+?)\s+of fill\s+(\d+)\s+delivered\s+([\d.,\sand]+?)\s*(?:uL|µL|ul)", re.I)
    for text in sheet.get("remarks", []):
        for m in pat.finditer(text):
            nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
            vols = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", m.group(3))]
            if len(nums) == len(vols):
                for n, v in zip(nums, vols):
                    out[(int(m.group(2)), n)] = v
    return out


# ---------------------------------------------------------------- isotherm
class Cell:
    """Overfilled cell. Titrant enters fill by fill from its own stock; forward: macromolecule in the cell."""

    def __init__(self, V_ul, cell_conc, forward, fills, dimer):
        self.V, self.X0, self.forward = V_ul, cell_conc, forward
        self.fills = fills                    # list of (titrant_conc, volumes)
        self.kdim, self.dhdis = dimer if dimer else (None, 0.0)

    def _eq(self, sites, lig, kd):
        if self.kdim is None or not self.forward:
            s = sites + lig + kd
            return 0.5 * (s - np.sqrt(max(s * s - 4.0 * sites * lig, 0.0))), 0.0
        if lig <= 0.0:
            return 0.0, 0.0
        k2 = self.kdim
        m = brentq(lambda m: m + 2 * m * m / k2 + sites * m / (kd + m) - lig, 0.0, lig, xtol=1e-14)
        return sites * m / (kd + m), m * m / k2

    def heats(self, g, kd, dh, h_mol):
        V = self.V
        X, Y, C, D = g * self.X0, 0.0, 0.0, 0.0
        out = []
        for conc, vols in self.fills:
            Ds = self._eq(0.0, conc, kd)[1] if self.forward else 0.0
            for v in vols:
                keep = 1.0 - v / V
                Xn, Yn = X * keep, Y * keep + conc * v / V
                sites, lig = (Xn, Yn) if self.forward else (Yn, Xn)
                Cn, Dn = self._eq(sites, lig, kd)
                q = (dh * V * (Cn - C * keep) - self.dhdis * (V * (Dn - D * keep) - v * Ds)) * UCAL_PER_KJ_X_PMOL
                out.append(q + h_mol * conc * v * 1e-3)
                X, Y, C, D = Xn, Yn, Cn, Dn
        return np.array(out)


# ionisation of buffering species at 25 C: dH (kJ/mol), dCp (J/(K mol)), as in run_format.md
BUFFER_ION = {"phosphate": (3.6, -187.0), "hepes": (20.4, 47.0), "tris": (47.45, -59.0), "mops": (21.1, 25.0),
              "pipes": (11.2, 22.0), "imidazole": (36.6, -9.0)}


def ionisation_dh(buffer, T_c):
    """Ionisation enthalpy of the run's buffer at T_c, carried from 25 C with its heat capacity change."""
    names = [b for b in BUFFER_ION if b in buffer.lower()]
    h25, cp = BUFFER_ION[names[0]]
    return h25 + cp * 1e-3 * (T_c - 25.0)


def analyse_run(folder):
    sheet = json.load(open(os.path.join(folder, "run_sheet.json")))
    period = float(sheet["instrument"]["filter_period_s"])
    V = 1000.0 * float(sheet["instrument"]["cell_volume_ml"])
    trace = load_record(os.path.join(folder, "thermogram.csv"))
    blank = sheet["blank"]
    btrace = load_record(os.path.join(folder, blank["file"]))

    ratios, fill_areas, taus = [], [], []
    for fill in sheet["fills"]:
        m = trace["fill"] == fill["fill"]
        fr = decompose_record(trace["time_s"][m], trace["dp_ucal_per_s"][m], period, fill["injections"], fill["pulses"])
        fill_areas.append(fr["areas"])
        taus.append(fr["tau_inst"])
        ratios += [-a / p["energy_ucal"] for a, p in zip(fr["pulse_areas"], fill["pulses"])]
    br = decompose_record(btrace["time_s"], btrace["dp_ucal_per_s"], period, blank["injections"], blank["pulses"])
    ratios += [-a / p["energy_ucal"] for a, p in zip(br["pulse_areas"], blank["pulses"])]
    # a pulse whose heater was found to deliver more than stated calibrates against what it delivered
    scale = float(np.mean(ratios)) / remark_checks(sheet)[1] if ratios else 1.0

    cell_sol = sheet["solutions"][sheet["cell"]["solution"]]
    forward = sheet["compounds"][cell_sol["compound"]]["kind"] == "macromolecule"
    lig_name = cell_sol["compound"] if not forward else sheet["solutions"][sheet["fills"][0]["syringe_solution"]]["compound"]
    dim = sheet["compounds"][lig_name].get("dimer")
    T_cell = float(np.mean(trace["cell_temp_c"]))
    dimer = None
    if dim:
        # the dimer constants are stated at 25 C: van 't Hoff to the cell temperature
        dhd = float(dim["dh_dissociation_kj_per_mol"])
        dimer = (float(dim["kd_um"]) * float(np.exp(-dhd * 1e3 / 8.314462618 * (1 / (T_cell + 273.15) - 1 / 298.15))),
                 dhd)

    def titrant(sid):
        c = stock_um(sheet, sid)
        sol = sheet["solutions"][sid]
        return c * float(sol.get("active_fraction", 1.0)) if not forward else c

    cell_conc = stock_um(sheet, sheet["cell"]["solution"])
    if forward:
        cell_conc *= float(cell_sol.get("active_fraction", 1.0))
    delivered = volumes_from_notes(sheet)
    fills, first = [], []
    for fill in sheet["fills"]:
        vols = [delivered.get((fill["fill"], j["number"]), j["volume_ul"]) for j in fill["injections"]]
        fills.append((titrant(fill["syringe_solution"]), vols))
        first += [True] + [False] * (len(vols) - 1)
    q = np.concatenate(fill_areas) / scale
    use = ~np.array(first)
    vflat = np.concatenate([np.asarray(v, float) for _, v in fills])

    # blank: buffer in the cell, so only dimer dissociation (if any) and the dilution heat per mole
    bconc = titrant(blank["syringe_solution"])
    bvols = [j["volume_ul"] for j in blank["injections"]]
    bq = br["areas"] / scale
    bcell = Cell(V, 0.0, forward, [(bconc, bvols)], dimer)

    dim_part = bcell.heats(1.0, 1.0, 0.0, 0.0)[1:]          # no macromolecule: Kd plays no part
    amount = bconc * np.asarray(bvols[1:]) * 1e-3
    h_mol = float(np.sum((bq[1:] - dim_part) * amount) / np.sum(amount * amount))

    cell = Cell(V, cell_conc, forward, fills, dimer)

    def resid(p):
        kd = 10.0 ** p[1]
        mq = cell.heats(p[0], kd, p[2], h_mol)
        sd = np.sqrt(0.03 ** 2 + (0.004 * mq) ** 2)
        return ((mq - q) / sd)[use]

    t0 = fills[0][0]
    dh0 = float(np.mean(q[use][:3] / (vflat[use][:3] * t0)) / UCAL_PER_KJ_X_PMOL)
    grid = np.linspace(np.log10(max(cell_conc, 1.0)) - 4.0, np.log10(max(cell_conc, 1.0)) + 1.0, 26 if dimer is None else 11)
    best = None
    for lkd in grid:
        try:
            s = least_squares(resid, [1.0, lkd, dh0], x_scale=[0.05, 0.2, max(abs(dh0), 1.0) * 0.1])
        except (ValueError, ZeroDivisionError):
            continue
        if best is None or s.cost < best.cost:
            best = s
    p = best.x
    # protons taken up from the buffer add its ionisation heat to every mole of complex formed
    n_h = float(sheet.get("binding", {}).get("proton_uptake", 0.0))
    dh_bind = float(p[2]) - n_h * ionisation_dh(sheet["buffer"], T_cell) if n_h else float(p[2])
    return {"kd_um": float(10.0 ** p[1]), "dh_kj_per_mol": dh_bind, "dh_observed_kj_per_mol": float(p[2]),
            "cell_factor": float(p[0]),
            "power_scale": scale, "tau_inst_s": [float(x) for x in taus] + [float(br["tau_inst"])],
            "dilution_ucal_per_nmol": h_mol, "chi2": float(2 * best.cost), "points": int(use.sum())}


def main():
    folder, out = sys.argv[1], sys.argv[2]
    res = analyse_run(folder)
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
