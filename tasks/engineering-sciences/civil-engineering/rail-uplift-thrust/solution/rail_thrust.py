"""Destressing temperature and axial thrust of continuous welded rail from one jack uplift test.

    python3 rail_thrust.py <test folder> <result file>

Model (displacement w measured from the bed, which follows the vertical alignment of radius R, crest positive):
EI w'''' + P w'' = -(q - P/R), q the catalogue weight per metre, P the thrust (compression positive); a rail bent to
the alignment carries the radial share P/R of its thrust as well as its weight.
- Free lift-off (half-length a inside half the undone length L): w'(0) = 0, w(a) = w'(a) = w''(a) = 0, mid shear F/2;
  a grid of a gives the (rise, F) curve.
- Held span (the free lift-off point would lie beyond L): the first fastened sleeper holds the rail, w(L) = w'(L) = 0,
  with w(0) = rise; F is linear in the rise, so two solves per thrust value cover all held steps.
UpliftFit.thrust() minimises the squared force misfit over all steps (50 kN scan, then bounded refinement).
T_destress = mean probe temperature - P / (E A alpha).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

E, ALPHA, G = 210e9, 1.15e-5, 9.80665
SECTIONS = {"60E1": (76.7, 3038.3, 60.21), "54E1": (69.77, 2337.9, 54.77), "49E1": (62.92, 1816.0, 49.39)}


def shape_rows(x, P, EI, q):
    """(4x4 homogeneous rows w..w''', particular w..w''') at x for EI w'''' + P w'' = -q"""
    if abs(P) < 1e-6:
        hom = np.array([[1, x, x * x, x ** 3], [0, 1, 2 * x, 3 * x * x], [0, 0, 2, 6 * x], [0, 0, 0, 6]], float)
        part = np.array([-q * x ** 4 / (24 * EI), -q * x ** 3 / (6 * EI), -q * x * x / (2 * EI), -q * x / EI])
        return hom, part
    k = np.sqrt(abs(P) / EI)
    if P > 0:
        c, s = np.cos(k * x), np.sin(k * x)
        hom = np.array([[1, x, c, s], [0, 1, -k * s, k * c], [0, 0, -k * k * c, -k * k * s], [0, 0, k ** 3 * s, -k ** 3 * c]])
    else:
        c, s = np.cosh(k * x), np.sinh(k * x)
        hom = np.array([[1, x, c, s], [0, 1, k * s, k * c], [0, 0, k * k * c, k * k * s], [0, 0, k ** 3 * s, k ** 3 * c]])
    return hom, np.array([-q * x * x / (2 * P), -q * x / P, -q / P, 0.0])


class UpliftFit:
    SPANS = np.linspace(0.6, 26.0, 420)

    def __init__(self, folder):
        sheet = json.loads((folder / "test_sheet.json").read_text())
        with open(folder / "uplift_log.csv") as f:
            log = [(float(r["jack_rise_mm"]), float(r["cell_reading_mv"])) for r in csv.DictReader(f)]
        self.rise = np.array([a for a, _ in log]) / 1000.0
        self.force = sheet["load_cell"]["factor_N_per_mV"] * np.array([b for _, b in log])
        area, inertia, mass = SECTIONS[sheet["section"]]
        self.area, self.EI, self.weight = area * 1e-4, E * inertia * 1e-8, mass * G
        self.radius = sheet["alignment"]["vertical_radius_m"]
        self.half_undone = sheet["fastenings_undone_m"] / 2.0
        self.t_rail = float(np.mean(sheet["probe_temperatures_C"]))

    def load(self, P):
        return self.weight - (P / self.radius if self.radius else 0.0)

    def _solve(self, P, q, end, end_rows, rise=None):
        H0, h0 = shape_rows(0.0, P, self.EI, q)
        He, he = shape_rows(end, P, self.EI, q)
        lhs, rhs = np.zeros((5, 5)), np.zeros(5)
        lhs[0, :4], rhs[0] = H0[1], -h0[1]
        lhs[1, :4], lhs[1, 4], rhs[1] = self.EI * H0[3] + P * H0[1], -0.5, -(self.EI * h0[3] + P * h0[1])
        for j, r in enumerate(end_rows):
            lhs[2 + j, :4], rhs[2 + j] = He[r], -he[r]
        if rise is not None:
            lhs[4, :4], rhs[4] = H0[0], rise - h0[0]
        sol = np.linalg.solve(lhs, rhs)
        return sol[4], H0[0] @ sol[:4] + h0[0]

    def predict(self, P):
        q = self.load(P)
        free = np.array([self._solve(P, q, a, (0, 1, 2)) for a in self.SPANS])
        order = np.argsort(free[:, 1])
        out = np.interp(self.rise, free[order, 1], free[order, 0])
        held = np.interp(self.rise, free[order, 1], self.SPANS[order]) > self.half_undone
        if held.any():
            f0 = self._solve(P, q, self.half_undone, (0, 1), rise=0.0)[0]
            f1 = self._solve(P, q, self.half_undone, (0, 1), rise=1.0)[0]
            out[held] = f0 + self.rise[held] * (f1 - f0)
        return out

    def misfit(self, p_kn):
        return float(np.sum((self.predict(p_kn * 1e3) - self.force) ** 2))

    def thrust(self):
        scan = np.arange(-900.0, 900.1, 50.0)
        start = scan[int(np.argmin([self.misfit(p) for p in scan]))]
        return minimize_scalar(self.misfit, bounds=(start - 50, start + 50), method="bounded", options={"xatol": 0.02}).x


def reduce(folder):
    fit = UpliftFit(Path(folder))
    p_kn = fit.thrust()
    return fit.t_rail - p_kn * 1e3 / (E * fit.area * ALPHA), p_kn


def main():
    t_d, p_kn = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps({"destress_temp_c": round(t_d, 3), "rail_thrust_kn": round(p_kn, 2)}) + "\n")


if __name__ == "__main__":
    main()
