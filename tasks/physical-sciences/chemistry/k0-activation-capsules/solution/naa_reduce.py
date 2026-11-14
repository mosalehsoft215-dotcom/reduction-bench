"""Reduces one k0-NAA capsule to the sodium and manganese mass fractions of its sample.

    python3 naa_reduce.py <capsule directory> <result.json>

1. Every peak area is divided by its count's pulser recovery, pulser_area / (pulser_hz t_live): pile-up removes
   gamma and pulser peaks alike, beyond what the live-time clock restores.  Every count of every item then gives
   the specific saturation count rate of its line,
   A = N / (t_live S D C m), D = exp(-lambda t_d), C = (1 - exp(-lambda t_real)) / (lambda t_real), m the item's
   mass, and S from the reactor's power log: S = (lambda / P_n) sum over the logged minutes of P_i (exp(-lambda
   (t_irr - t_i - 60 s)) - exp(-lambda (t_irr - t_i))) / lambda, the power held over each minute (for a steady
   reactor this is 1 - exp(-lambda t_irr)).
   In a sample count that also holds the 27Mg 1014.52 keV line, the 27Mg 843.76 keV peak share inside the Mn-56
   region of interest (846.76 +- 4 keV, Gaussian peaks of the stated FWHM) is computed from it and removed from the
   Mn-56 area first.  Repeated counts of one line are averaged.
2. The gold rate at a position is interpolated linearly along the capsule axis between the two gold wires; that
   gives the Cd wire's.  The sample, the foil and the wire share one inner holder whose recorded position is not
   used: the pure iron wire's 59Fe, read as iron through the k0 equation, must be 1 g/g, which fixes the gold rate
   where the holder was (A_Au,holder = A_Fe-59,wire / k0,Au(Fe) x (f + Q0,Au) / (f + Q0,Fe) x eps(E_Au) / eps(E_Fe)).
3. The channel's f during the irradiation from the cadmium ratio: R_Cd = A_Au(bare, at the Cd wire's position) /
   A_Au(Cd-covered), f = (R_Cd - 1) Q0,Au (alpha = 0, cadmium transmission 1); capsule.json's f is not used.
   A line's element mass fraction (g/g) is A / A_Au(z) / k0,Au(a) x (f + Q0,Au) / (f + Q0,a) x eps(E_Au) / eps(E_a).
4. The aluminium foil and the iron wire beside the sample are pure Al and Fe: the 24Na of the foil and the 56Mn of
   the wire, read as Na and Mn in the same way, are the apparent Na per gram of Al and Mn per gram of Fe from fast
   neutrons. The sample's Al (from 28Al) and Fe (from 59Fe) times these are subtracted from its apparent Na and Mn.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

NUC = {"Au-198": (2.6947 * 86400, 1.0, 15.71), "Na-24": (14.997 * 3600, 4.68e-2, 0.59),
       "Mn-56": (2.5789 * 3600, 0.496, 1.053), "Al-28": (2.2414 * 60, 1.75e-2, 0.71),
       "Fe-59": (44.495 * 86400, 7.77e-5, 0.975), "Mg-27": (9.458 * 60, 5.33e-5, 0.64)}
MG_OTHER = (843.76, 1.37e-4)                     # second 27Mg line: energy keV, k0 against Au
ROI_HALF = 4.0                                   # keV


def fwhm(e_kev):
    return 1.0 + 0.0011 * e_kev


def share_inside(line, centre):
    sig = fwhm(line) / math.sqrt(8 * math.log(2))
    return 0.5 * (math.erf((centre + ROI_HALF - line) / (sig * math.sqrt(2)))
                  - math.erf((centre - ROI_HALF - line) / (sig * math.sqrt(2))))


def reduce_capsule(folder, gradient=True, interference=True, pulser=True, magnesium=True, power=True, cadmium=True,
                   holder="wire"):
    """Options other than the defaults are the shortcuts tabulated by authoring/evidence/bench.py (power=False: the
    steady-reactor saturation from the clock's irradiation time; cadmium=False: capsule.json's f; holder="sheet": the
    gold interpolated at the holder's recorded position; pulser="worked": the live time with the worked capsules'
    pile-up share 0.5 d^2 instead of the pulser)."""
    folder = Path(folder)
    sheet = json.loads((folder / "capsule.json").read_text())
    ef = sheet["efficiency"]
    eff = lambda e: np.exp(ef["a0"] + ef["a1"] * np.log(e / 1000.0) + ef["a2"] * np.log(e / 1000.0) ** 2)
    rates, lines = {}, {}
    with open(folder / "counts.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    with open(folder / "reactor_power.csv", newline="") as fh:
        log = [(float(r["minute"]), float(r["power_kw"])) for r in csv.DictReader(fh)]
    t_irr = float(sheet["irradiation_s"])

    def saturation(lam):
        if not power:
            return 1 - np.exp(-lam * t_irr)
        out = 0.0
        for minute, kw in log:
            t0 = minute * 60.0
            if t0 >= t_irr:
                continue
            t1 = min(t0 + 60.0, t_irr)
            out += kw / sheet["nominal_power_kw"] * (np.exp(-lam * (t_irr - t1)) - np.exp(-lam * (t_irr - t0)))
        return out

    def factor(nuc, r):
        lam = np.log(2) / NUC[nuc][0]
        t_real, t_live = float(r["real_s"]), float(r["live_s"])
        s = saturation(lam)
        d = np.exp(-lam * float(r["decay_s"]))
        c = (1 - np.exp(-lam * t_real)) / (lam * t_real)
        return t_live * s * d * c * sheet["masses_g"][r["item"]]

    def area(r):
        n = float(r["net_area"])
        if pulser == "worked":
            n /= 1 - 0.5 * (1 - float(r["live_s"]) / float(r["real_s"])) ** 2
        elif pulser:
            n *= sheet["pulser_hz"] * float(r["live_s"]) / float(r["pulser_area"])
        return n

    by_count = {}
    for r in rows:
        by_count.setdefault(r["count"], {})[r["nuclide"]] = r
    for r in rows:
        n = area(r)
        if magnesium and r["nuclide"] == "Mn-56" and "Mg-27" in by_count[r["count"]]:
            m = by_count[r["count"]]["Mg-27"]
            a_mg = area(m) / factor("Mg-27", m)                 # 1014.52 keV, per gram, saturation
            a_other = (a_mg * MG_OTHER[1] / NUC["Mg-27"][1] * eff(MG_OTHER[0]) / eff(float(m["line_kev"])))
            n -= share_inside(MG_OTHER[0], float(r["line_kev"])) * a_other * factor("Mg-27", r)
        rates.setdefault((r["item"], r["nuclide"]), []).append(n / factor(r["nuclide"], r))
        lines[r["nuclide"]] = float(r["line_kev"])
    rate = {k: float(np.mean(v)) for k, v in rates.items()}
    z = sheet["axis_mm"]
    lo, hi = rate[("au_bottom", "Au-198")], rate[("au_top", "Au-198")]

    def axis_gold(item):
        if not gradient:
            return 0.5 * (lo + hi)
        u = (z[item] - z["au_bottom"]) / (z["au_top"] - z["au_bottom"])
        return lo + u * (hi - lo)

    f = sheet["f_thermal_epithermal"]
    if cadmium:
        r_cd = axis_gold("au_cd") / rate[("au_cd", "Au-198")]
        f = (r_cd - 1.0) * NUC["Au-198"][2]
    k0_fe, q0_fe = NUC["Fe-59"][1], NUC["Fe-59"][2]
    holder_gold = (rate[("fe_wire", "Fe-59")] / k0_fe * (f + NUC["Au-198"][2]) / (f + q0_fe)
                   * eff(lines["Au-198"]) / eff(lines["Fe-59"]))

    def gold(item):
        if holder == "wire" and item in ("sample", "al_foil", "fe_wire"):
            return holder_gold
        return axis_gold(item)

    def fraction(item, nuc):
        k0, q0 = NUC[nuc][1], NUC[nuc][2]
        return (rate[(item, nuc)] / gold(item) / k0 * (f + NUC["Au-198"][2]) / (f + q0)
                * eff(lines["Au-198"]) / eff(lines[nuc]))

    na, mn = fraction("sample", "Na-24"), fraction("sample", "Mn-56")
    if interference:
        per_al = fraction("al_foil", "Na-24")            # apparent Na per gram of Al
        per_fe = fraction("fe_wire", "Mn-56")            # apparent Mn per gram of Fe
        na -= fraction("sample", "Al-28") * per_al
        mn -= fraction("sample", "Fe-59") * per_fe
    return {"sodium_ug_per_g": float(na * 1e6), "manganese_ug_per_g": float(mn * 1e6)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_capsule(sys.argv[1])))
