"""Noise redraws of every session, calibrated by both programs.

    python redraws.py 20         -> redraws.txt, redraws.json      (second_solver.py needs SciPy)

The physical parameters of each session stay fixed; thermal noise, the drive phase, the bath readings and the power
meter readings are drawn again with seed 50000 + 31 i + session index. For each program, session and figure the
largest relative error over the redraws is kept, in percent.
"""
import json
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HOME = Path(__file__).resolve().parent
TOP = HOME.parents[1]
ORDER = ("P1", "P2", "H1", "H2", "H3", "H4", "H5")


def redraw(job):
    key, i = job
    for sub in ("authoring/provenance", "authoring/evidence", "solution"):
        if str(TOP / sub) not in sys.path:
            sys.path.append(str(TOP / sub))
    import file_session
    import second_solver
    import trap_reduce
    with tempfile.TemporaryDirectory() as tmp:
        where = Path(tmp) / "session"
        truth = file_session.file_session(key, where, noise_seed=50000 + 31 * i + ORDER.index(key))
        a = trap_reduce.reduce(where)
        second_solver.main(str(where), str(Path(tmp) / "b.json"))
        b = json.loads((Path(tmp) / "b.json").read_text())
    return key, {f: 100 * abs(a[f] / truth[f] - 1) for f in truth}, {f: 100 * abs(b[f] / truth[f] - 1) for f in truth}


if __name__ == "__main__":
    count = int(sys.argv[1])
    largest = {}
    jobs = [(k, i) for k in ORDER for i in range(count)]
    with ProcessPoolExecutor(max_workers=6) as pool:
        for key, a, b in pool.map(redraw, jobs):
            for program, errs in (("reference", a), ("second", b)):
                for f, v in errs.items():
                    largest[f"{program}|{key}|{f}"] = max(largest.get(f"{program}|{key}|{f}", 0.0), v)
    text = [f"largest relative error (%) in {count} noise redraws per session: per watt, reference stiffness, "
            "sensitivity, radius"]
    for program in ("reference", "second"):
        for key in ORDER:
            vals = [largest[f"{program}|{key}|{f}"] for f in ("stiffness_per_watt_pn_nm_w", "reference_stiffness_pn_nm",
                                                              "reference_sensitivity_nm_v", "bead_radius_nm")]
            text.append(f"{program} {key}: " + " ".join(f"{v:.2f}" for v in vals))
    (HOME / "redraws.txt").write_bytes(("\n".join(text) + "\n").encode())
    (HOME / "redraws.json").write_bytes((json.dumps(largest, indent=1, sort_keys=True) + "\n").encode())
    print("\n".join(text))
