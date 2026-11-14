"""Run other_analysis.py on every programme of a root holding truth.json; write other_analysis.csv.
usage: python3 other_analysis_all.py <programme root> <limits json> <out csv>"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from other_analysis import analyse  # noqa: E402

root, lim_path, out_path = sys.argv[1:4]
limits = json.load(open(lim_path))["limits_percent"]
truth = json.load(open(os.path.join(root, "truth.json")))
lines = []
for key in sorted(truth):
    res = analyse(os.path.join(root, key))
    devs = [100 * (res[k] / truth[key][k] - 1) for k in limits]
    over = max(abs(d) / limits[k] for d, k in zip(devs, limits))
    lines.append([key] + [f"{d:+.3f}" for d in devs] + [f"{over:.2f}", res["points"]])
    print(*lines[-1], flush=True)
with open(out_path, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["programme", "n_dev_percent", "Q_dev_percent", "rate_dev_percent", "worst_over_limit", "points"])
    w.writerows(lines)
