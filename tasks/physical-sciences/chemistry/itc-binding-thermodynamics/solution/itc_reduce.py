"""Reduce one titration-calorimetry run to its dissociation constant and binding enthalpy.

    python3 itc_reduce.py <run folder> <answer.json>

Steps:
1. Each record (every syringe fill, and the blank) is fitted as a whole: baseline (piecewise linear, knots every
   600 s) plus one peak per injection, each peak the response of a first-order instrument lag to a release spread
   over the injection and mixed in with a short lag of its own; the calibration pulse, where the record has one,
   is a further peak with the instrument lag only. Peak areas and baseline enter linearly; the two lags are found
   by least squares over the whole record (variable projection), so they come from the peaks themselves.
2. The power scale is the pulse area over its stated energy, or the last check when the run has no pulse.
3. The blank's heats give the titrant's heat of dilution per microlitre; without a blank it is fitted with the
   binding constants.
4. Binding isotherm: overfilled cell of fixed volume, each injection expelling its own volume of the
   pre-injection solution; one site per active macromolecule; the titrant concentration is taken as stated (the
   active part, for a macromolecule) and a factor on the cell species' concentration is fitted. The first injection
   of every fill is left out of the fit; its volume is kept.
"""
import json
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

UCAL_PER_KJ_PMOL = 1e-3 / 4.184


def read_trace(path):
    raw = np.genfromtxt(path, delimiter=",", names=True)
    return raw


def integrated_step(x, ta, tb):
    """Second integral of the unit impulse response of lags ta, tb (tb = 0: single lag), zero for x <= 0."""
    x = np.maximum(x, 0.0)
    if tb <= 0.0:
        return x - ta + ta * np.exp(-x / ta)
    if abs(ta - tb) < 1e-6 * ta:
        tb = ta * (1 + 1e-6)
    return x - ta - tb + (ta * ta * np.exp(-x / ta) - tb * tb * np.exp(-x / tb)) / (ta - tb)


def peak_column(t, period, start, dur, ta, tb):
    """Unit-area peak averaged over each filter period that ends at t."""
    S = integrated_step
    a = t - start
    return (S(a, ta, tb) - S(a - period, ta, tb) - S(a - dur, ta, tb) + S(a - dur - period, ta, tb)) / (period * dur)


class TraceFit:
    def __init__(self, t, y, period, peaks, pulse=None, knot_s=600.0):
        self.t, self.y, self.period = t, y, period
        self.peaks = peaks                 # list of (start, duration)
        self.pulse = pulse                 # (start, duration) or None
        n_k = max(2, int(np.ceil((t[-1] - t[0]) / knot_s)) + 1)
        knots = np.linspace(t[0], t[-1], n_k)
        eye = np.eye(n_k)
        self.base = np.array([np.interp(t, knots, eye[k]) for k in range(n_k)]).T

    def design(self, ta, tm):
        cols = [peak_column(self.t, self.period, s, d, ta, tm) for s, d in self.peaks]
        if self.pulse:
            cols.append(peak_column(self.t, self.period, self.pulse[0], self.pulse[1], ta, 0.0))
        return np.column_stack(cols + [self.base])

    def solve(self, ta, tm):
        A = self.design(ta, tm)
        coef, *_ = np.linalg.lstsq(A, self.y, rcond=None)
        return coef, self.y - A @ coef

    def fit(self):
        best = None
        for ta in (6.0, 10.0, 16.0, 25.0, 40.0, 60.0, 90.0):
            r = self.solve(ta, 4.0)[1]
            ss = float(r @ r)
            if best is None or ss < best[0]:
                best = (ss, ta)
        p0 = np.log([best[1], 4.0])
        sol = least_squares(lambda p: self.solve(*np.exp(p))[1], p0, x_scale=[0.2, 0.3],
                            bounds=([np.log(2.0), np.log(0.3)], [np.log(400.0), np.log(60.0)]))
        ta, tm = np.exp(sol.x)
        coef, res = self.solve(ta, tm)
        n = len(self.peaks)
        areas = coef[:n]
        pulse_area = coef[n] if self.pulse else None
        return dict(tau_inst=ta, tau_mix=tm, areas=areas, pulse_area=pulse_area, rms=float(res.std()))


def fit_record(t, y, period, injections, pulse):
    peaks = [(j["start_s"], j["duration_s"]) for j in injections]
    pul = (pulse["start_s"], pulse["duration_s"]) if pulse else None
    return TraceFit(t, y, period, peaks, pul).fit()


class Isotherm:
    """Overfilled-cell titration: cell species total X, titrant total Y; complex from the exact quadratic."""

    def __init__(self, V_ul, cell_um, titrant_um, titrant_is_ligand, fills):
        self.V, self.X0, self.Ys = V_ul, cell_um, titrant_um
        self.lig_titrant = titrant_is_ligand
        self.vols = [np.asarray(f, float) for f in fills]

    def heats(self, g, kd, dh):
        V, Ys = self.V, self.Ys
        X, Y = g * self.X0, 0.0
        C = 0.0
        out = []
        for vols in self.vols:
            for v in vols:
                keep = 1.0 - v / V
                Xn, Yn = X * keep, Y * keep + Ys * v / V
                s = Xn + Yn + kd
                Cn = 0.5 * (s - np.sqrt(max(s * s - 4.0 * Xn * Yn, 0.0)))
                out.append(dh * V * (Cn - C * keep) * UCAL_PER_KJ_PMOL)
                X, Y, C = Xn, Yn, Cn
        return np.array(out)


def reduce_run(folder):
    sheet = json.load(open(os.path.join(folder, "run_sheet.json")))
    inst = sheet["instrument"]
    period = float(inst["filter_period_s"])
    trace = read_trace(os.path.join(folder, "thermogram.csv"))
    calib = sheet["calibration"]
    pulse = calib.get("pulse")

    areas, vols, first_flags, taus = [], [], [], []
    scale = None
    for fill in sheet["fills"]:
        m = trace["fill"] == fill["fill"]
        own_pulse = pulse if (pulse and fill["fill"] == 1) else None
        fr = fit_record(trace["time_s"][m], trace["dp_ucal_per_s"][m], period, fill["injections"], own_pulse)
        taus.append(fr["tau_inst"])
        areas += list(fr["areas"])
        vols.append([j["volume_ul"] for j in fill["injections"]])
        first_flags += [True] + [False] * (len(fill["injections"]) - 1)
        if own_pulse:
            scale = -fr["pulse_area"] / own_pulse["energy_ucal"]
    if scale is None:
        scale = float(calib["last_check"]["power_scale"])
    q = np.array(areas) / scale

    blank = sheet.get("blank")
    dil = None
    if blank:
        b = read_trace(os.path.join(folder, blank["file"]))
        br = fit_record(b["time_s"], b["dp_ucal_per_s"], period, blank["injections"], None)
        bq = br["areas"][1:] / scale
        bv = np.array([j["volume_ul"] for j in blank["injections"]][1:])
        dil = float(np.sum(bq) / np.sum(bv))

    V = 1000.0 * float(inst["cell_volume_ml"])
    cell, syr = sheet["cell"], sheet["syringe"]
    if syr["species"] == "ligand":
        titrant = float(syr["conc_um"])
        cell_conc = float(cell["conc_um"]) * float(cell["active_fraction"])
    else:
        titrant = float(syr["conc_um"]) * float(syr["active_fraction"])
        cell_conc = float(cell["conc_um"])
    iso = Isotherm(V, cell_conc, titrant, syr["species"] == "ligand", vols)
    vflat = np.concatenate(iso.vols)
    use = ~np.array(first_flags)

    def model(p):
        h = dil if dil is not None else p[3]
        return iso.heats(p[0], 10.0 ** p[1], p[2]) + h * vflat

    def resid(p):
        mq = model(p)
        sd = np.sqrt(0.03 ** 2 + (0.004 * mq) ** 2)
        return ((mq - q) / sd)[use]

    # start values from the data: dH from the early heats, Kd scanned on a grid
    early = q[use][:3] - (dil or 0.0) * vflat[use][:3]
    dh0 = float(np.mean(early / (vflat[use][:3] * titrant)) / UCAL_PER_KJ_PMOL * 1e0)
    best = None
    for lkd in np.linspace(np.log10(cell_conc) - 4.0, np.log10(cell_conc) + 1.0, 26):
        p0 = [1.0, lkd, dh0] + ([] if dil is not None else [float(np.mean(q[-3:] / vflat[-3:]))])
        try:
            s = least_squares(resid, p0, x_scale=[0.05, 0.2, max(abs(dh0), 1.0) * 0.1] + ([] if dil is not None else [0.01]))
        except ValueError:
            continue
        if best is None or s.cost < best.cost:
            best = s
    p = best.x
    return {"kd_um": float(10.0 ** p[1]), "dh_kj_per_mol": float(p[2]), "cell_factor": float(p[0]),
            "power_scale": float(scale), "tau_inst_s": [float(x) for x in taus],
            "dilution_ucal_per_ul": float(dil if dil is not None else p[3]),
            "chi2": float(2 * best.cost), "points": int(use.sum())}


def main():
    folder, out = sys.argv[1], sys.argv[2]
    res = reduce_run(folder)
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
