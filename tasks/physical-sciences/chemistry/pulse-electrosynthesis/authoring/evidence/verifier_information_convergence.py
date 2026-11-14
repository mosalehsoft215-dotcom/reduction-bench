"""Refinement of the verifier's profiled information matrix and outcomes for
the reference design: stretched grids n = 80, 160 against n = 320.
Run from the bundle root; writes verifier_information_convergence.txt."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests"))
from information import evaluate_pair  # noqa: E402
cell = json.load(open(os.path.join(ROOT, "tests", "case", "cell.json")))
scen = json.load(open(os.path.join(ROOT, "tests", "case", "scenarios.json")))
ref = json.load(open(os.path.join(HERE, "reference_output.json")))
progs = [np.array([[r["duration_s"], r["potential_V"]] for r in ref["programs"][k]]) for k in ("A", "B")]
out = []
for sname in sys.argv[1:]:
    pars = [scen[sname][k] for k in cell["scenario_parameter_order"]]
    res = {}
    for n in (80, 160, 320):
        t0 = time.time()
        res[n] = evaluate_pair(progs, pars, cell, n=n) + (time.time() - t0,)
    mt, bt, ht, _ = res[320]
    val, vec = np.linalg.eigh(ht); wh = (vec / np.sqrt(val)) @ vec.T
    for n in (80, 160, 320):
        m, b, h, sec = res[n]
        err = np.linalg.norm(wh @ (h - ht) @ wh, ord=2)
        out.append("%s n=%d: lambda_min %.4f, scaled H error vs n=320 %.2e, bins %.1e, outcomes %.1e, %.0f s"
                   % (sname, n, np.linalg.eigvalsh(h)[0], err, np.max(abs(b - bt)), np.max(abs(m - mt)), sec))
        print(out[-1], flush=True)
with open(os.path.join(HERE, "verifier_information_convergence_%s.txt" % "_".join(sys.argv[1:])), "w", newline="\n") as fh:
    fh.write("\n".join(out) + "\n")
