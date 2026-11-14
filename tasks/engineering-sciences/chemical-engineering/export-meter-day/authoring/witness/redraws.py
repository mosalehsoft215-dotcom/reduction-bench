"""Instrument-noise redraws: same seven gas days, fresh meter, transmitter, GC and laboratory noise.

    python redraws.py <scratch folder> <out.csv> <noise number> [noise number ...]

For each noise number the days are rewritten by authoring/skid/make_days.py --noise N into the scratch
folder, and the reference program and the second solver are run on all seven. One row per noise number,
program, day and key with the relative error in %.
"""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[1]
PROGRAMS = {"reference": BUNDLE / "solution" / "meter_day.py", "second": HERE / "second_solver.py"}


def main():
    scratch, out = Path(sys.argv[1]), Path(sys.argv[2])
    with open(out, "w", newline="\n") as fh:
        wr = csv.writer(fh, lineterminator="\n")
        wr.writerow(["noise", "program", "day", "key", "error_pct"])
        for noise in map(int, sys.argv[3:]):
            days = scratch / f"noise_{noise}"
            subprocess.run([sys.executable, str(BUNDLE / "authoring" / "skid" / "make_days.py"), "--noise",
                            str(noise), "--out", str(days)], check=True, cwd=BUNDLE / "authoring" / "skid",
                           stdout=subprocess.DEVNULL)
            truth = json.loads((days / "truth.json").read_text())
            for label, prog in PROGRAMS.items():
                for day, tr in truth.items():
                    with tempfile.TemporaryDirectory() as tmp:
                        ans = Path(tmp) / "a.json"
                        subprocess.run([sys.executable, str(prog), str(days / day), str(ans)], check=True,
                                       capture_output=True)
                        got = json.loads(ans.read_text())
                    for k, v in tr.items():
                        wr.writerow([noise, label, day, k, f"{(got[k] / v - 1) * 100:+.4f}"])
                    fh.flush()
            print("noise", noise, "done", flush=True)


if __name__ == "__main__":
    main()
