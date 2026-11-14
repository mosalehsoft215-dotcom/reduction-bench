"""Redraw all seven lines with fresh instrument noise and class errors, then score both solvers.

    python redraws.py [N_SEEDS] [SCRATCH_DIR]    -> writes redraw_spread.txt next to this file

Each seed regenerates every line from its fixed design (authoring/provenance), runs
solution/reduce_run.py and second_solver.py, and records the error of each figure against
that draw's truth (CE in percentage points, energy in percent).
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parent.parent
sys.path.insert(0, str(BUNDLE / "authoring" / "provenance"))
sys.path.insert(0, str(BUNDLE / "solution"))
sys.path.insert(0, str(HERE))
import write_lines  # noqa: E402
import reduce_run  # noqa: E402
import second_solver  # noqa: E402


def errors(res, tr):
    return (res["current_efficiency_pct"] - tr["current_efficiency_pct"],
            100.0 * (res["dc_energy_kwh_per_t"] / tr["dc_energy_kwh_per_t"] - 1.0))


def main():
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    scratch = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.mkdtemp(prefix="redraw_"))
    worst = {}
    lines_out = []
    for seed in range(1, n_seeds + 1):
        out = scratch / f"seed{seed:02d}"
        truths = write_lines.make_all(out, out, out / "truth.json", 7000 + seed)
        cells = []
        for tag, tr in truths.items():
            a = errors(reduce_run.Line(out / tag).reduce(), tr)
            b = errors(second_solver.solve(out / tag), tr)
            for name, e in (("reference", a), ("second", b)):
                key = (name, tag)
                w = worst.get(key, (0.0, 0.0))
                worst[key] = (max(w[0], abs(e[0])), max(w[1], abs(e[1])))
            cells.append(f"{tag} {a[0]:+.3f}/{a[1]:+.3f} {b[0]:+.3f}/{b[1]:+.3f}")
        lines_out.append(f"seed {7000 + seed}: " + " | ".join(cells))
        print(lines_out[-1], flush=True)
        shutil.rmtree(out, ignore_errors=True)
    report = ["Redraws of all seven lines (CE error in points / energy error in %), reference then second solver",
              *lines_out, "", "Worst absolute error per line over the redraws:"]
    for name in ("reference", "second"):
        for tag in sorted({k[1] for k in worst}):
            ce, en = worst[(name, tag)]
            report.append(f"{name:9s} {tag}: CE {ce:.3f} pt, energy {en:.3f} %")
        ce = max(v[0] for k, v in worst.items() if k[0] == name)
        en = max(v[1] for k, v in worst.items() if k[0] == name)
        report.append(f"{name:9s} all lines: CE {ce:.3f} pt, energy {en:.3f} %")
    (HERE / "redraw_spread.txt").write_text("\n".join(report) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
