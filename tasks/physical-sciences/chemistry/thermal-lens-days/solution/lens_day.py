"""Reference reduction of one thermal-lens day: absorption coefficient and thermal diffusivity of the sample.

    python3 /app/solution/lens_day.py <day folder> <answer.json>

1. Probe geometry from the setup: zc = pi w0^2 / lambda_p, z1 = sample - waist, m = (w(z1)/w_e)^2, V = z1/zc.
2. Every transient: the excitation onset t0 from the pump photodiode (first crossing of half its plateau), the probe
   normalised by its mean before t = 0, and the mode-mismatched model fitted for theta and tc.
3. The fraction eta of the meter reading that reaches the liquid from the reference liquids. A reference's alpha is
   its certificate unless the same-day absorbance disagrees by more than three standard uncertainties.
4. Temperature: with references in two solvents, a common offset to the cell thermometer is chosen so that both
   give the same eta (the guide: one eta for every run with the same optics); otherwise the cell reading.
5. Sample: P_abs = -theta k lambda_p / (dn/dT); alpha = -ln(1 - P_abs / (eta P_meter)) / l; D = w_e^2 / (4 tc).
LENS_ROUTE (comma-separated) switches single steps for authoring/proofs/missteps.py.
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import least_squares, minimize_scalar

ROUTE = set(filter(None, os.environ.get("LENS_ROUTE", "").split(",")))
LP = 632.8e-9


def csv_rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def shen(t, theta, tc, m, V):
    arg = 2 * m * V / (((1 + 2 * m) ** 2 + V ** 2) * tc / (2 * t) + 1 + 2 * m + V ** 2)
    return (1 - 0.5 * theta * np.arctan(arg)) ** 2


class Day:
    def __init__(self, folder):
        self.d = Path(folder)
        self.setup = json.loads((self.d / "setup.json").read_text())
        self.runs = csv_rows(self.d / "runs.csv")
        self.props = csv_rows(self.d / "properties.csv")
        self.absb = {r["liquid"]: (float(r["absorbance_per_cm"]), float(r["uncertainty"])) for r in csv_rows(self.d / "absorbance.csv")}
        pr = self.setup["probe"]
        w0, zc_lam = pr["waist_radius_um"] * 1e-6, LP
        zc = math.pi * w0 ** 2 / zc_lam
        z1 = (pr["sample_position_mm"] - pr["waist_position_mm"]) * 1e-3
        self.we = self.setup["excitation"]["waist_radius_um_at_sample"] * 1e-6
        self.m = (w0 ** 2 * (1 + (z1 / zc) ** 2)) / self.we ** 2
        self.V = z1 / zc

    def prop(self, solvent, t):
        rows = [r for r in self.props if r["solvent"] == solvent]
        T = [float(r["t_C"]) for r in rows]
        return (float(np.interp(t, T, [float(r["k_W_per_m_K"]) for r in rows])),
                float(np.interp(t, T, [float(r["dn_dT_per_K"]) for r in rows])))

    def fit(self, run):
        rows = csv_rows(self.d / "transients" / f"{run}.csv")
        t = np.array([float(r["t_ms"]) for r in rows]) * 1e-3
        pv = np.array([float(r["probe_V"]) for r in rows])
        pump = np.array([float(r["pump_V"]) for r in rows])
        I0 = pv[t < 0].mean()
        if "trigger_as_zero" in ROUTE:
            lo, hi = 0.0, 1e-7
        else:
            plateau = np.median(pump[-200:])
            first = t[np.argmax(pump > 0.5 * plateau)]
            lo, hi = first - 1.0e-4, first          # onset lies within one sample before the first lit one
        sel = t > hi
        x, y = t[sel], pv[sel] / I0
        res = least_squares(lambda p: shen(x - p[2], p[0], p[1], self.m, self.V) - y, x0=[0.3, 0.008, 0.5 * (lo + hi)],
                            x_scale=[0.1, 0.005, 1e-5], bounds=([0, 1e-5, lo], [5, 1, hi]))
        return res.x[:2]

    def reduce(self):
        liq = self.setup["liquids"]
        fits = {r["run"]: self.fit(r["run"]) for r in self.runs}
        refs = [r for r in self.runs if liq[r["liquid"]]["role"] == "reference"]

        def ref_alpha(name):
            cert = liq[name]["certified_alpha_per_cm"]
            if name in self.absb and "trust_certificate" not in ROUTE:
                a, u = self.absb[name]
                if abs(a * math.log(10) - cert) > 3 * u * math.log(10):
                    return a * math.log(10)
            return cert

        def etas(dT):
            out = {}
            for r in refs:
                solvent = liq[r["liquid"]]["solvent"]
                k, dndt = self.prop(solvent, float(r["cell_C"]) + dT)
                A, l = ref_alpha(r["liquid"]), float(r["path_cm"])
                P_abs = -fits[r["run"]][0] * k * LP / dndt
                out.setdefault(solvent, []).append(P_abs / (float(r["pump_meter_mW"]) * 1e-3 * (1 - math.exp(-A * l))))
            return out

        dT = 0.0
        if len(etas(0.0)) > 1 and "cell_thermometer" not in ROUTE:
            spread = lambda d: float(np.std([np.mean(v) for v in etas(d).values()]))
            dT = minimize_scalar(spread, bounds=(-15, 15), method="bounded", options={"xatol": 1e-4}).x
        eta = float(np.mean([e for v in etas(dT).values() for e in v]))
        alphas, Ds = [], []
        for r in self.runs:
            if liq[r["liquid"]]["role"] != "sample":
                continue
            k, dndt = self.prop(liq[r["liquid"]]["solvent"], float(r["cell_C"]) + dT)
            theta, tc = fits[r["run"]]
            frac = (-theta * k * LP / dndt) / (eta * float(r["pump_meter_mW"]) * 1e-3)
            l = 1.0 if "unit_path" in ROUTE else float(r["path_cm"])
            alphas.append(frac / l if "linear_absorption" in ROUTE else -math.log(1 - frac) / l)
            Ds.append(self.we ** 2 / (4 * tc))
        return {"alpha_per_cm": float(np.mean(alphas)), "diffusivity_mm2_per_s": float(np.mean(Ds)) * 1e6}


def main():
    Path(sys.argv[2]).write_text(json.dumps(Day(sys.argv[1]).reduce(), indent=2) + "\n")


if __name__ == "__main__":
    main()
