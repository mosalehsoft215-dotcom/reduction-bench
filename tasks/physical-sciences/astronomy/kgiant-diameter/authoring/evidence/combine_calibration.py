"""Merges calibration_results_night_*.json into the graded quantities of the task: diameter and flux ratio averaged
over the three nights, separation and position angle per night.   -> calibration_results.json, .txt"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRUTH = json.loads((HERE.parent.parent / "tests" / "ground_truth.json").read_text())
REL = {"diameter_mas": 0.3, "separation_mas": 0.05, "flux_ratio_percent": 0.6}
ABS = 0.06
NIGHTS = ("night_2024", "night_2025", "night_2026")
per = {n: {r["route"]: r for r in json.loads((HERE / f"calibration_results_{n}.json").read_text())} for n in NIGHTS}
rows = []
for route in per["night_2026"]:
    if not all(route in per[n] for n in NIGHTS):
        continue
    v = {}
    for n in NIGHTS:
        y = n.split("_")[1]
        v[f"separation_mas_{y}"] = per[n][route]["values"]["separation_mas"]
        v[f"position_angle_deg_{y}"] = per[n][route]["values"]["position_angle_deg"]
    v["diameter_mas"] = sum(per[n][route]["values"]["diameter_mas"] for n in NIGHTS) / 3
    v["flux_ratio_percent"] = sum(per[n][route]["values"]["flux_ratio_percent"] for n in NIGHTS) / 3
    e = {}
    for k in TRUTH:
        if k.startswith("position_angle_deg"):
            e[k] = (v[k] - TRUTH[k] + 180) % 360 - 180
        else:
            e[k] = 100 * (v[k] / TRUTH[k] - 1)
    ratio = {k: abs(e[k]) / (ABS if k.startswith("position") else REL[k.rsplit("_", 1)[0] if k[-4:].isdigit() else k]) for k in TRUTH}
    worst = max(ratio, key=ratio.get)
    rows.append(dict(label=per["night_2026"][route]["label"], route=route, values=v, errors=e, worst=worst,
                     worst_over_band=ratio[worst], inside=all(x <= 1 for x in ratio.values())))
(HERE / "calibration_results.json").write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8", newline="\n")
lines = ["the three nights together; errors in percent (position angles in degrees), keys in the order of tests/ground_truth.json"]
for r in rows:
    lines.append(f"{r['label']:4s} {r['route']:110s} " + " ".join(f"{r['errors'][k]:+7.3f}" for k in TRUTH)
                 + f"  worst {r['worst']} {r['worst_over_band']:.1f}x band  {'INSIDE' if r['inside'] else 'outside'}")
(HERE / "calibration_results.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
print("\n".join(lines))
