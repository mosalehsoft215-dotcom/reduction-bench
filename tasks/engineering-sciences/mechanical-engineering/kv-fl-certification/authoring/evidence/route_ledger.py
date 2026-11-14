"""Every route on every specimen, and both correct reductions over noise redraws.

    python route_ledger.py <bundle> <scratch> [redraws]

Writes authoring/evidence/route_ledger.json; misses in allowances (|relative error| / limit, larger key).
  reference        solution/kv_fl_reduce.py (every rated-travel sample, converter damping fitted)
  second           authoring/evidence/second_capacity.py (damping from the meter's run-downs, 1 s blocks, Nelder-Mead)
  lab_convention   settled holds and their last 45 s as the lab reports them, dp zeros from the spells with the meter
                   at 0, the shown flow as it is, then the full fit (n, height, meter zero): the reference of the
                   earlier version
  undamped_samples the reference's sample fit with the shown flow taken as the flow (no damping)
  damped_holds     hold windows between throttle moves, each window mean lifted by tau x the shown flow's change over
                   it, dp zeros from the spells with the meter at 0
  as_logged        the sample fit with damping, every other reading taken as it stands: an abrupt change between the
                   limits, no transmitter height, no meter zero
  abrupt           everything right except that the change between the limits is abrupt (n fixed large)
  no_height        everything right except that the p1 transmitter is taken to read the tap pressure
  no_meter_zero    everything right except that the meter's zero is taken as 0
  one_dp_zero      everything right except that the first spell's dp zero serves the whole record
  one_column       everything right except that the upstream column is taken as the same all through the record
  no_boyle         everything right except that the air taken into the lines at a stop changes the readings with the
                   pump running by what it changed them with the pump off (no compression by the tap pressure)
  one_meter_zero   everything right except that one meter zero serves every pump run of the record
  shut_as_off      everything right except that every closed-valve spell is taken as pump-off (dp zero and column)
  runs_at_trips    everything right except that the meter's zero changes only where the pump stops, not where the
                   valve is shut with the pump running
  installed        the reference's fit, but Kv and FL taken as installed (no reducer correction)
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIMIT = {"kv_rated_m3h": 0.015, "fl_rated": 0.010}


def miss(got, truth):
    return round(max(abs(got[k] / truth[k] - 1.0) / LIMIT[k] for k in LIMIT), 3)


def routes(vr, folder):
    out = {"as_logged": vr.reduce(folder, (), boyle=False, stretch_zero=False), "abrupt": vr.reduce(folder, ("e", "q0")),
           "no_height": vr.reduce(folder, ("n", "q0")), "no_meter_zero": vr.reduce(folder, ("n", "e"))}
    out["lab_convention"] = vr.reduce(folder, damping=False, samples=False)
    out["undamped_samples"] = vr.reduce(folder, damping=False)
    out["damped_holds"] = vr.reduce(folder, samples=False)
    out["no_boyle"] = vr.reduce(folder, boyle=False)
    out["one_meter_zero"] = vr.reduce(folder, stretch_zero=False)
    keep = vr.spells

    def flat(k):
        def f(rows, rate):
            got = list(keep(rows, rate))
            got[k] = np.full(len(got[k]), got[k][0])
            return tuple(got)
        return f
    vr.spells = flat(1)
    out["one_dp_zero"] = vr.reduce(folder)
    vr.spells = flat(2)
    out["one_column"] = vr.reduce(folder)

    def as_off(rows, rate):
        index, dz, pr, run = keep(rows, rate)
        shut = rows[:, 1] < 1.0
        found, start = [], None
        for i, ok in enumerate(np.append(shut, False)):
            if ok and start is None:
                start = i
            elif not ok and start is not None:
                if i - start >= 20 * rate:
                    found.append((start, float(np.median(rows[start + int(5 * rate):i, 4])),
                                  float(np.median(rows[start + int(5 * rate):i, 3]))))
                start = None
        return run, np.array([f[1] for f in found]), np.array([f[2] for f in found]), run
    vr.spells = as_off
    out["shut_as_off"] = vr.reduce(folder)
    vr.spells = lambda rows, rate: (lambda g: (g[0], g[1], g[2], g[0]))(keep(rows, rate))
    out["runs_at_trips"] = vr.reduce(folder)
    vr.spells = keep
    sheet, rows = vr.load(folder)
    s2 = json.loads(json.dumps(sheet))
    s2["valve_end_bore_mm"] = s2["loop_pipe"]["inside_diameter_mm"]
    kv, fl, _ = vr.sample_fit(s2, rows)
    out["installed"] = {"kv_rated_m3h": kv, "fl_rated": fl}
    return out


def main():
    root, scratch = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    redraws = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    sys.path[:0] = [str(root / "solution"), str(HERE)]
    import kv_fl_reduce as vr
    import second_capacity as sc
    truth = json.loads((root / "environment/data/loop_logs/worked_capacities.json").read_text())
    truth.update(json.loads((root / "tests/specimen_truth.json").read_text()))
    report = {"limits": LIMIT, "specimens": {}, "redraws": {}}
    for name, t in truth.items():
        f = root / "environment/data/loop_logs" / name
        f = f if f.is_dir() else root / "tests/sealed_specimens" / name
        row = {"reference": miss(vr.reduce(f), t), "second": miss(sc.reduce(f), t)}
        row.update({k: miss(v, t) for k, v in routes(vr, f).items()})
        report["specimens"][name] = row
        print(name, row, flush=True)
    for k in range(1, redraws + 1):
        out = scratch / f"valve_redraw_{k}"
        subprocess.run([sys.executable, str(root / "authoring/provenance/make_specimens.py"), str(out), "--redraw", str(k)], check=True)
        rt = json.loads((out / "truth.json").read_text())
        report["redraws"][str(k)] = {n: {"reference": miss(vr.reduce(out / n), rt[n]), "second": miss(sc.reduce(out / n), rt[n]),
                                         "lab_convention": miss(vr.reduce(out / n, damping=False, samples=False), rt[n]),
                                         "undamped_samples": miss(vr.reduce(out / n, damping=False), rt[n]),
                                         "no_boyle": miss(vr.reduce(out / n, boyle=False), rt[n]),
                                         "one_meter_zero": miss(vr.reduce(out / n, stretch_zero=False), rt[n])} for n in rt}
        print("redraw", k, report["redraws"][str(k)], flush=True)
    report["worst_over_redraws"] = {w: max(v[w] for r in report["redraws"].values() for v in r.values())
                                    for w in ("reference", "second")} if redraws else {}
    print("worst over redraws", report["worst_over_redraws"])
    (HERE / "route_ledger.json").write_text(json.dumps(report, indent=1) + "\n", newline="\n")


if __name__ == "__main__":
    main()
