"""Fresh-noise redraws: every test regenerated with new instrument noise and reduced by the reference and by
second_solver.py. Noise 0 is the bundle itself, and only there is the quick replica scored too.

    python redraws.py <work folder> <csv out> <noise> [noise ...]

One CSV row per program, test and control head, rewritten after each noise level so partial runs are kept.
"""
import csv, json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "authoring" / "plant" / "make_surge_tests.py"
SOLVERS = (("reference", ROOT / "solution" / "surge_line.py"),
           ("second_solver", ROOT / "authoring" / "checks" / "second_solver.py"),
           ("replica_low", ROOT / "authoring" / "checks" / "replica_low" / "surge_line.py"))


def answer(program, test_dir, scratch):
    done = subprocess.run([sys.executable, str(program), str(test_dir), str(scratch)], capture_output=True)
    return json.loads(scratch.read_text()) if done.returncode == 0 else {}


work, target = Path(sys.argv[1]), Path(sys.argv[2])
table = [["noise", "program", "test", "key", "found", "true", "error_pct"]]
for level in ("0", *sys.argv[3:]):
    tests = work / f"noise_{level}"
    subprocess.run([sys.executable, str(GENERATOR), "--noise", level, "--out", str(tests)], check=True,
                   stdout=subprocess.DEVNULL)
    truths = json.loads((tests / "truth.json").read_text())
    for label, program in SOLVERS:
        if label == "replica_low" and level != "0":
            continue
        for test in truths:
            got = answer(program, tests / test, tests / "_scratch.json")
            for key, true in truths[test].items():
                found = got.get(key, float("nan"))
                table.append([level, label, test, key, found, true, f"{(found - true) / true * 100:+.3f}"])
    with target.open("w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows(table)
    print("noise", level, "written", flush=True)
