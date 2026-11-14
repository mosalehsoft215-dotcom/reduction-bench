"""Uplift bench.  Every number is an error in units of the hidden limits (1.5 K on destress_temp_c, 25 kN on
rail_thrust_kn): first the reference and five simplified routes on the seven shipped tests, then the reference and
second_uplift.py on eight redraws.

Plain routes (all reuse the reference's UpliftFit, so only the modelling choice differs):
  free lift-off, weight only     - lift-off model at every step, catalogue weight, no alignment term
  free lift-off + alignment      - adds the radial share of the thrust, still free lift-off everywhere
  held span, weight only         - switches to the fastened sleeper when needed, no alignment term
  first six steps                - free lift-off, weight only, fitted to the six lowest rises
  worked-test offset             - free lift-off, weight only, shifted by its mean error on the two worked tests

    python bench.py <bundle root>    (result: authoring/evidence/bench.json)
"""
from pathlib import Path
import sys, json, tempfile, subprocess

import numpy as np

LIMIT = {"destress_temp_c": 1.5, "rail_thrust_kn": 25.0}


def plain(S, folder, alignment, held, first=None):
    fit = S.UpliftFit(folder)
    if not alignment:
        fit.radius = None
    if not held:
        fit.half_undone = 1e9
    fit.rise, fit.force = fit.rise[:first], fit.force[:first]
    p_kn = fit.thrust()
    return {"destress_temp_c": fit.t_rail - p_kn * 1e3 / (S.E * fit.area * S.ALPHA), "rail_thrust_kn": p_kn}


def excess(got, want):
    return max(abs(got[key] - want[key]) / LIMIT[key] for key in LIMIT)


def run(root):
    b = Path(root).resolve()
    sys.path.insert(0, str(b / "solution"))
    import rail_thrust as S
    data = b / "environment" / "data" / "uplift_tests"
    worked = json.loads((data / "reference_results.json").read_text())
    truth = {**json.loads((b / "tests" / "uplift_truth.json").read_text()), **worked}
    where = {u: data / u if (data / u).exists() else b / "tests" / "sealed_uplifts" / u for u in truth}
    ref = {u: dict(zip(("destress_temp_c", "rail_thrust_kn"), S.reduce(where[u]))) for u in truth}
    base = {u: plain(S, where[u], False, False) for u in truth}
    shift = {k: np.mean([worked[u][k] - base[u][k] for u in worked]) for k in LIMIT}
    routes = {"reference": ref, "free lift-off, weight only": base,
              "free lift-off + alignment": {u: plain(S, where[u], True, False) for u in truth},
              "held span, weight only": {u: plain(S, where[u], False, True) for u in truth},
              "first six steps": {u: plain(S, where[u], False, False, 6) for u in truth},
              "worked-test offset": {u: {k: base[u][k] + shift[k] for k in LIMIT} for u in truth}}
    res = {"routes": {n: {u: round(float(excess(r[u], truth[u])), 2) for u in sorted(truth)} for n, r in routes.items()},
           "redraws": {}}
    for n, row in res["routes"].items():
        print(n, row, flush=True)
    worst = [0.0, 0.0]
    for k in range(1, 9):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            subprocess.run([sys.executable, str(b / "authoring/provenance/make_uplifts.py"), str(tmp), "--redraw", str(k)],
                           check=True, capture_output=True)
            tt = json.loads((tmp / "truth.json").read_text())
            tt.update(json.loads((tmp / "data/uplift_tests/reference_results.json").read_text()))
            row = {}
            for u in sorted(tt):
                folder = tmp / "data/uplift_tests" / u
                folder = folder if folder.exists() else tmp / "sealed" / u
                m1 = excess(dict(zip(("destress_temp_c", "rail_thrust_kn"), S.reduce(folder))), tt[u])
                subprocess.run([sys.executable, str(b / "authoring/evidence/second_uplift.py"), str(folder), str(tmp / "o.json")],
                               check=True)
                m2 = excess(json.loads((tmp / "o.json").read_text()), tt[u])
                row[u] = [round(float(m1), 3), round(float(m2), 3)]
                worst = [max(worst[0], m1), max(worst[1], m2)]
            res["redraws"][f"draw {k}"] = row
            print(k, row, flush=True)
    res["worst_redraw"] = {"reference": round(float(worst[0]), 3), "second": round(float(worst[1]), 3)}
    (Path(__file__).parent / "bench.json").write_text(json.dumps(res, indent=1) + "\n", newline="\n")
    print(res["worst_redraw"])


if __name__ == "__main__":
    run(sys.argv[1])
