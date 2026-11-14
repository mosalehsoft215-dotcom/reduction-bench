"""Shared table writer for the replica and misreading scripts: runs one pipeline configuration on all seven runs and
reports, per run, the error of each graded value in per cent and in allowances (error over the limit)."""
import csv
import json
from pathlib import Path

from variants import BUNDLE, patient

TRUTH = json.loads((BUNDLE / "authoring" / "simulator" / "truth.json").read_text())
LIMITS = json.loads((BUNDLE / "tests" / "constants.json").read_text())["limits_percent"]


def score_rows(label, options, chi2_ref=None):
    rows = []
    for run, truth in TRUTH.items():
        out = patient(run, **options)
        row = {"reading": label, "run": run}
        worst = 0.0
        for key, lim in LIMITS.items():
            err = (out[key] / truth[key] - 1.0) * 100.0
            row[f"{key}_err_pct"] = round(err, 2)
            worst = max(worst, abs(err) / lim)
        row["allowances"] = round(worst, 2)
        row["passes"] = worst <= 1.0
        row["chi2"] = round(out["chi2"], 1)
        if chi2_ref:
            row["chi2_correct"] = round(chi2_ref[run], 1)
        row["repaired"] = out["repaired"]
        rows.append(row)
        print(row, flush=True)
    return rows


def write(path, rows):
    with open(path, "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
