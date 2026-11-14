"""Switches on the reference, both programs on the shipped stations, and blind answers.

    python3 switch_runs.py          -> switches.txt

Each row leaves one step of the reference out. A cell is direction error (degrees, modulo 180) / delay error (s), with
! when it falls outside the windows of 11 degrees and 0.3 s.
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

EV = Path(__file__).resolve().parent
TOP = EV.parents[1]
sys.path[:0] = [str(TOP / "solution"), str(TOP / "authoring/provenance"), str(EV)]
import eigen_split  # noqa: E402
import make_stations  # noqa: E402
import split_measure  # noqa: E402

ROWS = {"reference": {}, "calibration trusted, no phase check": {"use_phase": False},
        "sheet trusted, no phase check": {"use_calibration": False, "use_phase": False},
        "one wiring for the whole deployment": {"one_wiring": True}, "rotation only, no reversed channel": {"flip_allowed": False},
        "no orientation step": {"use_orientation": False}, "mean of single-event results": {"stack": False}}
WIDE = {"dir": 11.0, "dt": 0.3}


def where(name):
    worked = make_stations.STATIONS[name]["role"] == "worked"
    return TOP / ("environment/data/vaults" if worked else "tests/sealed") / name


def work(job):
    name, row = job
    return eigen_split.main(where(name)) if row == "second" else split_measure.measure(where(name), **ROWS[row])


def error(answer, truth):
    return ((answer["fast_direction_deg"] - truth[0] + 90.0) % 180.0 - 90.0, answer["split_delay_s"] - truth[1])


if __name__ == "__main__":
    names = list(make_stations.STATIONS)
    truth = {n: (make_stations.STATIONS[n]["phi"], make_stations.STATIONS[n]["dt"]) for n in names}
    rows = list(ROWS) + ["second"]
    with ProcessPoolExecutor(6) as pool:
        got = dict(zip([(n, r) for r in rows for n in names], pool.map(work, [(n, r) for r in rows for n in names])))
    lines = ["switch | " + " | ".join(names)]
    for r in rows:
        cells = []
        for n in names:
            d, e = error(got[(n, r)], truth[n])
            cells.append(f"{d:+.0f}/{e:+.2f}" + ("!" if abs(d) > WIDE["dir"] or abs(e) > WIDE["dt"] else ""))
        lines.append(f"{r} | " + " | ".join(cells))
    sealed = [n for n in names if make_stations.STATIONS[n]["role"] != "worked"]
    inside = lambda a, n: abs((a[0] - truth[n][0] + 90) % 180 - 90) <= WIDE["dir"] and abs(a[1] - truth[n][1]) <= WIDE["dt"]
    worked = [truth[n] for n in names if n not in sealed]
    guesses = {"copy of SV-36": worked[0], "copy of SV-38": worked[1],
               "mean of the worked stations": (100.6, sum(w[1] for w in worked) / 2),
               "other mean of the worked stations": (10.6, sum(w[1] for w in worked) / 2)}
    for name, g in guesses.items():
        lines.append(f"{name}: inside on {sum(inside(g, n) for n in sealed)} of 5")
    best = max(sum(inside((p, q / 100.0), n) for n in sealed) for p in range(0, 180, 2) for q in range(40, 240, 5))
    lines.append(f"one shared answer, best over a grid: inside on {best} of 5")
    (EV / "switches.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
