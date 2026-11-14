"""Repeatability: every well filed again with fresh measurement noise, reduced by both programs.

    python redraws.py 20       -> redraws.txt, redraws.json   (SciPy needed for second_solver.py)

Formations, drilling plan, probe and calibration stay as modelled; the resistance noise, the encoder noise and the
core-sample scatter are drawn anew with seed 91000 + 37 i + well index. For each program, well and figure the
largest error over the redraws is kept, in the figure's own unit.
"""
import json
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

EVIDENCE = Path(__file__).resolve().parent
BUNDLE = EVIDENCE.parents[1]
KEYS = ("P1", "P2", "H1", "H2", "H3", "H4", "H5")
FIGURES = ("basal_heat_flow_mw_m2", "temperature_at_td_c", "report_gradient_k_km", "report_conductivity_w_mk")


def one(job):
    key, i = job
    for sub in ("authoring/provenance", "authoring/evidence", "solution"):
        if str(BUNDLE / sub) not in sys.path:
            sys.path.append(str(BUNDLE / sub))
    import file_well
    import second_solver
    import well_reduce
    with tempfile.TemporaryDirectory() as tmp:
        here = Path(tmp)
        truth = file_well.file_well(key, here / "w", noise_seed=91000 + 37 * i + KEYS.index(key))
        ref = well_reduce.reduce(here / "w")
        second_solver.main(str(here / "w"), str(here / "s.json"))
        alt = json.loads((here / "s.json").read_text())
    return key, {f: abs(ref[f] - truth[f]) for f in FIGURES}, {f: abs(alt[f] - truth[f]) for f in FIGURES}


if __name__ == "__main__":
    n = int(sys.argv[1])
    top = {}
    with Pool(6) as pool:
        for key, ref, alt in pool.imap_unordered(one, [(k, i) for k in KEYS for i in range(n)]):
            for prog, errs in (("reference", ref), ("second", alt)):
                for f, e in errs.items():
                    top[f"{prog}|{key}|{f}"] = max(top.get(f"{prog}|{key}|{f}", 0.0), e)
    lines = [f"largest error in {n} noise redraws per well (mW/m2, C, K/km, W/m/K)"]
    for prog in ("reference", "second"):
        for key in KEYS:
            lines.append(f"{prog} {key}: " + "  ".join(f"{f.split('_')[0]}_{f.split('_')[1]} {top[f'{prog}|{key}|{f}']:.3f}"
                                                         for f in FIGURES))
    (EVIDENCE / "redraws.txt").write_bytes(("\n".join(lines) + "\n").encode())
    (EVIDENCE / "redraws.json").write_bytes((json.dumps(top, indent=1, sort_keys=True) + "\n").encode())
    print("\n".join(lines))
