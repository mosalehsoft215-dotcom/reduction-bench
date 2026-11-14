"""Fresh noise draws of all seven samples (generator draws 1-30), reduced by the reference and by other_route.py.

    python noise_redraws.py      -> noise_redraws.csv

Each draw rebuilds every record of every sample with new noise (authoring/forge/samples.py, same truths and
instrument settings) in a temporary folder, then both programs reduce it.
"""
import csv
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[1]
sys.path.insert(0, str(BUNDLE / "solution"))
sys.path.insert(0, str(BUNDLE / "authoring/forge"))
sys.path.insert(0, str(HERE))
import other_route  # noqa: E402
import samples  # noqa: E402
import vsm_reduce  # noqa: E402

DRAWS = range(1, 31)
OUT = HERE / "noise_redraws.csv"
if len(sys.argv) == 3:                       # a slice of the draws, for running slices in parallel
    DRAWS = range(int(sys.argv[1]), int(sys.argv[2]) + 1)
    OUT = HERE / f"noise_redraws_part{sys.argv[1]}.csv"
truth = json.loads((BUNDLE / "authoring/forge/truth.json").read_text())
limit = json.loads((BUNDLE / "tests/answer_key.json").read_text())["tolerance_percent"]
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for draw in DRAWS:
        for s in samples.SAMPLES:
            where = Path(tmp) / f"d{draw}" / s["folder"]
            samples.build(s, draw, where)
            t = truth[s["folder"]]
            for solver, fn in (("reference", vsm_reduce.reduce_sample), ("other_route", other_route.reduce)):
                got = fn(where)
                err = {k: 100.0 * (got[k] / t[k] - 1.0) for k in limit}
                rows.append([draw, s["folder"], solver, f"{err['ms_kA_per_m']:+.3f}", f"{err['hk_kA_per_m']:+.3f}",
                             f"{max(abs(err[k]) / limit[k] for k in limit):.3f}"])
        print("draw", draw, "done", flush=True)
with OUT.open("w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["draw", "sample", "solver", "ms_error_percent", "hk_error_percent", "allowances"])
    w.writerows(rows)
for solver in ("reference", "other_route"):
    sub = [r for r in rows if r[2] == solver]
    worst = max(sub, key=lambda r: float(r[5]))
    print(solver, "worst", worst, "runs", len(sub))
    for name in truth:
        print("   ", name, max(float(r[5]) for r in sub if r[1] == name))
