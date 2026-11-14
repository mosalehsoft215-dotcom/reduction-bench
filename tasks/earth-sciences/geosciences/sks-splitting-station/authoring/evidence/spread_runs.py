"""Spread over new seeds: every station regenerated twenty times and measured by both programs.

    python3 spread_runs.py <scratch folder>

Seed j of a station is (its seed + 7907 j). spread.json keeps every answer and truth. spread.txt gives per station
the largest direction error (degrees, taken modulo 180) and delay error (s) of each program.
"""
import json
import shutil
import sys
from multiprocessing import Pool
from pathlib import Path

EV = Path(__file__).resolve().parent
TOP = EV.parents[1]
sys.path[:0] = [str(TOP / "solution"), str(TOP / "authoring" / "provenance"), str(EV)]
import eigen_split  # noqa: E402
import make_stations  # noqa: E402
import split_measure  # noqa: E402

N_SEEDS = 20


def one(job):
    name, j, scratch = job
    where = Path(scratch) / f"{name}-{j}"
    truth = make_stations.write(name, where, seed=make_stations.STATIONS[name]["seed"] + 7907 * j)
    row = {"station": name, "seed_index": j, "truth": truth, "reference": split_measure.measure(where),
           "second": eigen_split.main(where)}
    shutil.rmtree(where, ignore_errors=True)
    return row


def misses(ans, truth):
    d = (ans["fast_direction_deg"] - truth["fast_direction_deg"] + 90.0) % 180.0 - 90.0
    return abs(d), abs(ans["split_delay_s"] - truth["split_delay_s"])


if __name__ == "__main__":
    scratch = Path(sys.argv[1])
    scratch.mkdir(parents=True, exist_ok=True)
    with Pool(8) as pool:
        rows = pool.map(one, [(n, j, str(scratch)) for n in make_stations.STATIONS for j in range(1, N_SEEDS + 1)])
    (EV / "spread.json").write_text(json.dumps(rows, indent=1))
    text = [f"{N_SEEDS} seeds per station; largest direction error (deg) and delay error (s)"]
    for name in make_stations.STATIONS:
        for prog in ("reference", "second"):
            m = [misses(r[prog], r["truth"]) for r in rows if r["station"] == name]
            text.append(f"{name} {prog:9s} {max(x[0] for x in m):6.2f} deg  {max(x[1] for x in m):5.3f} s")
    (EV / "spread.txt").write_text("\n".join(text) + "\n")
    print("\n".join(text))
