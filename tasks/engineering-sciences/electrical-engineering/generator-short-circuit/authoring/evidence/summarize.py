"""Print the numbers the README quotes, straight from the tables.

    python authoring/evidence/summarize.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

T = Path(__file__).resolve().parent / "tables"
KEYS = ["Rfd", "Lfd", "R1d", "L1d", "R2d", "L2d", "Lf1d", "Xd"]


sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabulate import allowances  # noqa: E402

ALLOW = allowances()


def rows(name):
    """Rows of a table, with the allowance column recomputed from err_pct and the current allowances (the table
    is rewritten when it differs, so the files always carry the allowances the tests use)."""
    p = T / name
    if not p.exists():
        return []
    with open(p) as fh:
        data = list(csv.DictReader(fh))
    if data and "err_pct" in data[0]:
        changed = False
        for r in data:
            new = f"{abs(float(r['err_pct'])) / 100.0 / ALLOW[r['key']]:.2f}" if r["err_pct"] else ""
            changed |= new != r["allowances"]
            r["allowances"] = new
        if changed:
            with open(p, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(data[0]), lineterminator="\n")
                w.writeheader()
                w.writerows(data)
    return data


def correct(name, who):
    data = rows(name)
    if not data:
        return
    worst = defaultdict(lambda: (0.0, ""))
    secs = []
    for r in data:
        if not r["err_pct"]:
            print(f"  {who}: {r['label']} {r['record']} has no answer")
            continue
        e = abs(float(r["err_pct"]))
        if e > worst[r["key"]][0]:
            worst[r["key"]] = (e, f"{r['label']} {r['record']}")
        secs.append(float(r["seconds"]))
    print(f"{who}: worst |error| per key over {len({r['label'] for r in data})} draws x 7 records")
    for k in KEYS:
        print(f"  {k:5s} {worst[k][0]:7.3f} %   ({worst[k][1]})")
    print(f"  wall time per record {min(secs):.0f}-{max(secs):.0f} s (authoring PC, shared)")


def routes():
    data = rows("wrong_routes.csv")
    if not data:
        return
    by = defaultdict(dict)
    for r in data:
        if r["allowances"]:
            cur = by[(r["label"], r["record"])].get("w", (0.0, ""))
            a = float(r["allowances"])
            if a > cur[0]:
                by[(r["label"], r["record"])]["w"] = (a, f"{r['key']} {float(r['err_pct']):+.2f} %")
        else:
            by[(r["label"], r["record"])]["w"] = (float("nan"), "no answer")
    labels = sorted({k[0] for k in by}, key=lambda s: [r["label"] for r in data].index(s))
    recs = ["W-1", "W-2", "H-1", "H-2", "H-3", "H-4", "H-5"]
    print("wrong routes: worst key in allowances per record")
    print("  " + " " * 22 + "".join(f"{r:>24s}" for r in recs))
    for lab in labels:
        cells = []
        for r in recs:
            a, txt = by[(lab, r)].get("w", (float("nan"), "-"))
            cells.append(f"{a:6.2f} ({txt})")
        print(f"  {lab:22s}" + "".join(f"{c:>24s}" for c in cells))


def constants():
    data = rows("constant_answers.csv")
    if data:
        low = min(data, key=lambda r: float(r["worst_allowances"]))
        print(f"constant answers: smallest miss {low['worst_allowances']} allowances ({low['constant']} on "
              f"{low['record']}, key {low['worst_key']})")


def probes():
    data = rows("judge_probes.csv")
    if data:
        ok = sum(1 for r in data if r["expected"] == r["reward"])
        print(f"judge probes: {ok}/{len(data)} as expected")
        for r in data:
            print(f"  {r['probe']:45s} expected {r['expected']} reward {r['reward']}  {r['seconds']} s")


if __name__ == "__main__":
    correct("reference.csv", "reference")
    correct("second_solver.csv", "second solver")
    routes()
    constants()
    probes()
