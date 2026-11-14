"""Fresh noise on every survey, reduced by the reference and by the independent program (SciPy needed).

    python redraws.py 30       -> redraws.txt
"""
import concurrent.futures as cf
import json
import pathlib
import sys
import tempfile

BUNDLE = pathlib.Path(__file__).resolve().parents[2]
LIMITS = {"gl_t": 0.0015, "b6_units": 0.15, "axis_x_mm": 0.015, "axis_y_mm": 0.015, "roll_mrad": 0.08}
ORDER = ("P1", "P2", "H1", "H2", "H3", "H4", "H5")


def redraw(task):
    key, k = task
    for sub in ("authoring/sim", "authoring/proof", "solution"):
        if str(BUNDLE / sub) not in sys.path:
            sys.path.insert(0, str(BUNDLE / sub))
    import coil_reduce
    import file_survey
    import independent
    with tempfile.TemporaryDirectory() as tmp:
        folder = pathlib.Path(tmp) / "s"
        truth = file_survey.file_survey(key, folder, noise_seed=52000 + 41 * k + ORDER.index(key))
        first = coil_reduce.reduce(folder)
        independent.main(str(folder), str(pathlib.Path(tmp) / "i.json"))
        second = json.loads((pathlib.Path(tmp) / "i.json").read_text())
    return key, {name: max(abs(ans[q] - truth[q]) / lim for q, lim in LIMITS.items())
                 for name, ans in (("reference", first), ("independent", second))}


if __name__ == "__main__":
    n = int(sys.argv[1])
    worst = {}
    with cf.ProcessPoolExecutor(max_workers=6) as pool:
        for key, errs in pool.map(redraw, [(key, k) for key in ORDER for k in range(n)]):
            for name, e in errs.items():
                worst[(name, key)] = max(worst.get((name, key), 0.0), e)
    text = [f"{n} noise redraws per survey; largest error over all five figures, in tolerance units"]
    for name in ("reference", "independent"):
        text.append(name + ": " + ", ".join(f"{key} {worst[(name, key)]:.2f}" for key in ORDER))
    (pathlib.Path(__file__).parent / "redraws.txt").write_bytes(("\n".join(text) + "\n").encode())
    print("\n".join(text))
