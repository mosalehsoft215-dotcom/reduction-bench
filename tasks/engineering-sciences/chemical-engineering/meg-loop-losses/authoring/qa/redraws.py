"""Fresh-noise copies of all seven loops (process noise and instrument noise redrawn), reduced by the reference and
by the second solver; deviations in units of the limits.

    python redraws.py <copies per loop> <workers>
"""
import csv
import json
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HERE / "authoring" / "model"))
sys.path.insert(0, str(HERE / "solution"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
LIM = {"formation_water_t_d": 0.03, "meg_loss_kg_d": 0.05}
KEYS = ["P1", "P2", "U1", "U2", "U3", "U4", "U5"]


def one(job):
    import write_records as W
    import meg_balance as R
    import second_solver as S2
    key, i = job
    seed = 50000 + 100 * i + KEYS.index(key)
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp) / key
        _, truth = W.write(key, folder, noise_seed=seed)
        R.main(str(folder), str(Path(tmp) / "ref.json"))
        S2.main(str(folder), str(Path(tmp) / "s2.json"))
        ref = json.loads((Path(tmp) / "ref.json").read_text())
        s2 = json.loads((Path(tmp) / "s2.json").read_text())
    row = {"loop": key, "copy": i}
    for k, lim in LIM.items():
        row[f"ref_{k}"] = (ref[k] / truth[k] - 1) / lim
        row[f"s2_{k}"] = (s2[k] / truth[k] - 1) / lim
    return row


def main(n, workers):
    jobs = [(k, i) for i in range(n) for k in KEYS]
    with Pool(workers) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    out = Path(__file__).parent / "redraws.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    lines = []
    for who in ("ref", "s2"):
        for k in LIM:
            vals = [abs(r[f"{who}_{k}"]) for r in rows]
            worst = max(rows, key=lambda r: abs(r[f"{who}_{k}"]))
            lines.append(f"{who} {k}: worst {max(vals):.2f} of the limit ({worst['loop']} copy {worst['copy']}), "
                         f"mean {sum(vals) / len(vals):.2f}")
    (Path(__file__).parent / "redraws.txt").write_text(
        f"{len(rows)} fresh-noise loops ({n} of each of the seven)\n" + "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
