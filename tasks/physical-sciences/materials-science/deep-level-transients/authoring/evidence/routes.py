"""Wrong routes, one step at a time: the reference with a single step replaced by the shortcut a hurried analyst
takes. Prints, per route and diode, the worst of the three errors in allowances (level 0.010 eV, density 6 %,
sigma 0.20 in log10) and writes routes.txt beside this file.

    python routes.py [<cases root with D-?? folders and truth.json>]
"""
import contextlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[1]
sys.path.insert(0, str(BUNDLE / "solution"))
import trap_reduce as R  # noqa: E402


def cases(root):
    if root:
        root = Path(root)
        truth = json.loads((root / "truth.json").read_text())
        return {k: root / k for k in truth}, truth
    truth = json.loads((BUNDLE / "environment/data/diodes/accepted.json").read_text())
    truth.update(json.loads((BUNDLE / "tests/diode_truth.json").read_text()))
    where = {k: (BUNDLE / "environment/data/diodes" / k) if k in ("D-21", "D-22") else (BUNDLE / "tests/diodes" / k)
             for k in truth}
    return where, truth


@contextlib.contextmanager
def patched(obj, name, value):
    old = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, old)


def misses(ans, tru):
    return (abs(ans["level_ev"] - tru["level_ev"]) / 0.010,
            abs(ans["density_cm3"] / tru["density_cm3"] - 1.0) / 0.06,
            abs(math.log10(ans["sigma_300k_cm2"] / tru["sigma_300k_cm2"])) / 0.20)


# ---- the shortcuts ---------------------------------------------------------------------------------------------------
def meter_as_junction(self, c_pf, g_us):
    return np.asarray(c_pf) * 1e-12


def one_kind_only(folder, talk=False):
    """Bias series never consulted: every centre treated as neutral after emission."""
    orig = R.Layer.pf_gamma
    with patched(R.Layer, "pf_gamma", lambda self, T, x, x2, w: np.zeros_like(x)):
        return R.reduce(folder)


def no_capture_barrier(bench, layer):
    return {}


def uniform_doping(bench, layer):
    lay = R.first_cv(bench)
    return lay.with_trap(layer.level, layer.nt, layer.donor)


def small_signal_density(layer, tr, fa=1.0):
    """N_T = 2 N_d (dC / C) W^2 / (x2^2 - x1^2) with dC from the 1/C^2 amplitude, N_d at the edge."""
    w_inf = tr["w_inf"]
    c_inf = R.EPS * tr["scale"] / w_inf
    c0 = 1.0 / math.sqrt(tr["a"] + tr["b"])
    dc = c_inf - c0
    return 2.0 * layer.nd(w_inf) * dc / c_inf * w_inf ** 2 / (tr["x2"] ** 2 - tr["x1"] ** 2)


def no_lambda_density(layer, tr, fa=1.0):
    c_inf = R.EPS * tr["scale"] / tr["w_inf"]
    c0 = 1.0 / math.sqrt(tr["a"] + tr["b"])
    return 2.0 * layer.nd(tr["w_inf"]) * (c_inf - c0) / c_inf


def arrhenius_without_t2(points):
    T = np.array([p["T"] for p in points])
    e = np.array([p["e0"] for p in points])
    X = 1.0 / (R.KB * T)
    coef = np.polyfit(X, np.log(e), 1)
    ea = -coef[0]
    gamma = R.vth(300.0) * R.nc(300.0) / 300.0 ** 2
    return float(ea), float(math.exp(coef[1]) / gamma / 300.0 ** 2)


def stage_is_sample(self):
    return lambda t: t


ROUTES = [
    ("reference", []),
    ("stage reading taken as the sample temperature", [(R.Bench, "sample_temperatures", stage_is_sample)]),
    ("meter C taken as the junction capacitance", [(R.Bench, "junction_c", meter_as_junction)]),
    ("every centre neutral after emission (bias series unused)", [(R, "reduce", None)]),
    ("capture series unused (no capture barrier)", [(R, "capture_rates", no_capture_barrier)]),
    ("uniform doping from the 1/C^2 line", [(R, "fit_cv", uniform_doping)]),
    ("small-signal density, lambda kept", [(R, "density_from", small_signal_density)]),
    ("density 2 N_d dC/C, no lambda", [(R, "density_from", no_lambda_density)]),
    ("Arrhenius of e_n without T^2", [(R, "arrhenius", arrhenius_without_t2)]),
]


def run_route(patches, folder):
    with contextlib.ExitStack() as st:
        for obj, name, val in patches:
            if val is None:
                return one_kind_only(folder)
            st.enter_context(patched(obj, name, val))
        return R.reduce(folder)


def main():
    where, truth = cases(sys.argv[1] if len(sys.argv) > 1 else None)
    names = sorted(truth)
    lines = ["route | " + " | ".join(names)]
    for label, patches in ROUTES:
        cells = []
        for n in names:
            m = misses(run_route(patches, where[n]), truth[n])
            cells.append(" ".join(f"{v:.1f}" for v in m) + ("!" if max(m) > 1.0 else ""))
        lines.append(label + " | " + " | ".join(cells))
        print(lines[-1], flush=True)
    if len(sys.argv) == 1:
        (HERE / "routes.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
