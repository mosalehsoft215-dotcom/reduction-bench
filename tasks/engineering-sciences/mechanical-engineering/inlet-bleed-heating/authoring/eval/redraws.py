"""Fresh-noise re-runs of all seven units, scored against the generator's truth (authoring evidence).

    python3 redraws.py OUT.json FIRST LAST [ROUTE ...]

For draws FIRST..LAST-1: the reference, the second reduction, and the reference along each named ROUTE
(comma-joined switches of bleed_check.py). OUT.json holds per program and unit the % errors [bleed_b, rise_a].
"""
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, "..", "plant", "inlet_sim.py")
REF = os.path.join(HERE, "..", "..", "solution", "bleed_check.py")
TWO = os.path.join(HERE, "second_route.py")
KEYS = ("bleed_b_kg_s", "rise_a_k")


def draw(job):
    k, routes = job
    res = {}
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([sys.executable, SIM, tmp, str(k)], check=True, capture_output=True)
        truth = json.load(open(os.path.join(tmp, "truth.json")))
        progs = [("second", TWO, "")] if os.environ.get("ONLY_SECOND") else [("reference", REF, ""), ("second", TWO, "")]
        for name, prog, route in progs + [(r, REF, r) for r in routes]:
            for tag, tru in truth.items():
                out = os.path.join(tmp, "a.json")
                p = subprocess.run([sys.executable, prog, os.path.join(tmp, tag), out], capture_output=True,
                                   env=dict(os.environ, BLEED_ROUTE=route))
                if p.returncode:
                    err = [None, None]
                else:
                    got = json.load(open(out))
                    err = [round(100 * (got[x] / tru[x] - 1), 3) for x in KEYS]
                res.setdefault(name, {}).setdefault(tag, []).append(err)
    return res


if __name__ == "__main__":
    dest, a, b, routes = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4:]
    total = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("WORKERS", "4"))) as pool:
        for part in pool.map(draw, [(k, routes) for k in range(a, b)]):
            for name, tags in part.items():
                for tag, errs in tags.items():
                    total.setdefault(name, {}).setdefault(tag, []).extend(errs)
    json.dump(total, open(dest, "w"), indent=0)
    for name, tags in total.items():
        cells = []
        for tag, errs in tags.items():
            ok = [e for e in errs if e[0] is not None]
            hi = [max(abs(e[i]) for e in ok) if ok else None for i in (0, 1)]
            lo = [min(abs(e[i]) for e in ok) if ok else None for i in (0, 1)]
            cells.append(f"{tag}:{hi[0]}/{hi[1]} (min {lo[0]}/{lo[1]})")
        print(f"{name[:44]:44s}", " ".join(cells), "crashes", sum(e[0] is None for v in tags.values() for e in v))
