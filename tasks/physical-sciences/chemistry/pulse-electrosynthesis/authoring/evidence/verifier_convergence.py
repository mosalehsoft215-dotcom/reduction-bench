"""Grid refinement study for the verifier's physics (tests/physics.py).

The submitted verifier used a uniform 80-cell grid; the run audit showed its
charge bins up to 1.2e-4 from the converged solution, i.e. at the 1e-4
tolerance itself. This script compares, for the reference design and two
adversarial pulse trains (potential steps 0.5 ms before and after bin edges),
in the fastest and the slowest scenario:
  - the submitted uniform 80-cell grid (physics_uniform80_submitted.py)
  - the stretched grids n = 80 (r = 1.06, used by the verifier), 160 (r = 1.03)
  - the truth: stretched n = 320 (r = 1.015), checked against uniform n = 1280.
Run from the bundle root: python authoring/evidence/verifier_convergence.py
Writes verifier_convergence.txt next to itself.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests"))
import physics as new  # noqa: E402

spec = importlib.util.spec_from_file_location("old", os.path.join(HERE, "physics_uniform80_submitted.py"))
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)

cell = json.load(open(os.path.join(ROOT, "tests", "case", "cell.json")))
scen = json.load(open(os.path.join(ROOT, "tests", "case", "scenarios.json")))
ref = json.load(open(os.path.join(HERE, "reference_output.json")))


def program(rows):
    return np.array([[r["duration_s"], r["potential_V"]] for r in rows])


def adversarial(offset):
    rows = [[0.05 + offset, -1.25]]
    for k in range(46):
        rows.append([0.1, -1.30 if k % 2 == 0 else -0.90])
    rows.append([6.0 - sum(r[0] for r in rows), -1.20])
    return np.array(rows)


designs = {"reference A": program(ref["programs"]["A"]), "reference B": program(ref["programs"]["B"]),
           "steps 0.5 ms before bin edges": adversarial(-0.0005),
           "steps 0.5 ms after bin edges": adversarial(+0.0005)}
out = []


def log(s):
    print(s, flush=True)
    out.append(s)


log("max |difference| from the truth (stretched n=320, r=1.015): charge bins; outcomes")
for sname in ("fast_transfer", "slow_transport"):
    pars = [scen[sname][k] for k in cell["scenario_parameter_order"]]
    for dname, prog in designs.items():
        t0 = time.time()
        mt, bt = new.trajectory(prog, pars, n=320, ratio=1.015)
        t_truth = time.time() - t0
        mu, bu = old.trajectory(prog, pars, n=1280)
        rows = [("uniform 1280 (check of the truth)", mu, bu)]
        t0 = time.time()
        m80, b80 = old.trajectory(prog, pars, n=80)
        t_old = time.time() - t0
        rows.append(("submitted uniform 80", m80, b80))
        t0 = time.time()
        ms, bs = new.trajectory(prog, pars)
        t_new = time.time() - t0
        rows.append(("verifier stretched 80, r=1.06", ms, bs))
        m16, b16 = new.trajectory(prog, pars, n=160, ratio=1.03)
        rows.append(("stretched 160, r=1.03", m16, b16))
        log("%s / %s (seconds: uniform80 %.1f, stretched80 %.1f, truth %.1f)"
            % (sname, dname, t_old, t_new, t_truth))
        for label, m, b in rows:
            log("   %-36s bins %.2e   outcomes %.2e" % (label, np.max(abs(b - bt)), np.max(abs(m - mt))))
with open(os.path.join(HERE, "verifier_convergence.txt"), "w", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
