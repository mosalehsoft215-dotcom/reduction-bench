"""Measure precision and discriminate shortcuts on records and independent count redraws."""
import argparse
import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / part) for part in ("solution", "authoring/provenance", "tests")]
from decay_fit import fit, read
from make_records import RECORDS, fabricate, write_capture, write_json
from split_decay import reduce_case
from test_afterglow import WINDOWS


def scaled(answer, truth):
    return {key: abs(answer[key] - truth[key]) / float(window)
            for key, window in WINDOWS.items()}


def measure(redraws):
    rows = {}
    for name, config in RECORDS.items():
        row = rows[name] = {"redraws": []}
        for redraw in range(redraws + 1):
            with tempfile.TemporaryDirectory() as location:
                folder = Path(location)
                info, counts, truth = fabricate(config, redraw)
                write_json(folder / "cell.json", info)
                write_capture(folder / "decays.npz", counts)
                args = read(folder)
                results = {}
                for label, calculate in (("reference", lambda: fit(*args)),
                                         ("splitting", lambda: reduce_case(folder))):
                    begin = time.perf_counter()
                    values = calculate()
                    results[label] = {"scaled_errors": scaled(values, truth),
                                      "seconds": round(time.perf_counter() - begin, 4)}
                row["redraws"].append(results)
                if redraw == 0:
                    raw = read(folder, normalize=False)
                    row["logged_values"] = {
                        "gauge_readings_as_pressures": scaled(fit(dict(args[0], evacuated_reading_torr=0.0), *args[1:]), truth),
                        "no_outside_chord_normalisation": scaled(fit(*raw), truth),
                        "both_as_logged": scaled(fit(dict(raw[0], evacuated_reading_torr=0.0), *raw[1:]), truth),
                    }
                if redraw == 0 and name not in {"AG-04", "AG-09"}:
                    info, observed, sigma = args
                    simplified = dict(info, ray_offsets_mm=[0.0], ray_weights=[1.0])
                    centers = np.mean(info["gates_ms"], axis=1)
                    midpoint = dict(info, gates_ms=np.column_stack((centers, centers)).tolist())
                    row["shortcuts"] = {
                        "single_ray": scaled(fit(simplified, observed, sigma), truth),
                        "gate_midpoint": scaled(fit(midpoint, observed, sigma), truth),
                        "no_exchange": scaled(fit(*args, free=[0, 1, 2, 4, 5], fixed={3: 0}), truth),
                        "no_pair_loss": scaled(fit(*args, free=[0, 1, 2, 3, 5], fixed={4: 0}), truth),
                        "linear_bulk": scaled(fit(*args, free=[0, 1, 2, 3, 4], fixed={5: 0}), truth),
                    }
        print(name, [round(max(max(d[s]["scaled_errors"].values()) for d in row["redraws"]), 3)
                     for s in ("reference", "splitting")], flush=True)
    return rows


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("output", type=Path)
    cli.add_argument("--redraws", type=int, default=6)
    options = cli.parse_args()
    write_json(options.output, measure(options.redraws))
