"""Every wrong reading of the README table on the runs it concerns, eight refits at a time ("+search": followed by
the one-knob repair of wrong_routes_one.py --search).

    python3 wrong_routes.py > wrong_routes.csv

Output in per cent of the true efficiency and in allowances (deviation / tolerance); a run fails when either
efficiency is more than one allowance off.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parents[1]
TOL = json.loads((BUNDLE / "tests" / "expected.json").read_text())["tolerance_percent"]
RUNS = {p.name: p for d in (BUNDLE / "environment" / "data" / "runs", BUNDLE / "tests" / "columns") for p in d.iterdir()}
PUBLIC = ("run_C201_0312", "run_C305_0605")
ROWS = [("none", r) for r in sorted(RUNS)] + [("raoult", r) for r in PUBLIC] + [("satreflux", r) for r in PUBLIC] + [
    ("isobaric", "run_C201_0312")] + [(t, r) for t in ("tau15", "tau30", "nolag", "steam_noloss") for r in PUBLIC] + [
    ("nocal", "h1"), ("check_inverted", "h1"), ("ua_zero", "h1"), ("ua_nosteam", "h1"), ("nocal+search", "h1"),
    ("nocal", "h2"), ("notc", "h2"), ("check_rho20", "h2"), ("check_inverted", "h2"), ("rho20", "h2"),
    ("satreflux", "h2"), ("nocal+search", "h2"),
    ("tau15", "h3"), ("tau30", "h3"), ("nolag", "h3"), ("molbasis", "h3"), ("noside", "h3"),
    ("tf_feedin", "h4"), ("tf_free", "h4"), ("feedfirst", "h4"), ("toppress", "h4"), ("tf_feedin+search", "h4"),
    ("steam_nosub", "h5"), ("steam_nosub+search", "h5"), ("steam_noloss", "h5"), ("r_nosteam", "h5"),
    ("noflash", "h5"), ("flash_std", "h5"), ("antoine_public", "h5")]
env = dict(os.environ, OPENBLAS_NUM_THREADS="1")


def one(row):
    trap, run = row
    args = [trap.split("+")[0]] + (["--search"] if trap.endswith("+search") else [])
    p = subprocess.run([sys.executable, str(HERE / "wrong_routes_one.py"), str(RUNS[run])] + args, env=env,
                       capture_output=True, text=True)
    return row, (p.stdout.strip().splitlines() or ["failed"])[-1]


print("reading,run,emv_rectifying_dev_pct,emv_stripping_dev_pct,allowances_worst,chi2,run_fails,knob_kept")
with ThreadPoolExecutor(8) as pool:
    for (trap, run), line in pool.map(one, ROWS):
        if "rect" not in line:
            print(f"{trap},{run},,,,,yes (program stops),")
            continue
        w = line.split()
        dr, ds, chi2 = float(w[3]), float(w[6]), float(w[9])
        kept = line.split(" kept ", 1)[1].split(" (")[0] if " kept " in line else ""
        a = max(abs(dr) / TOL["emv_rectifying"], abs(ds) / TOL["emv_stripping"])
        print(f"{trap},{run},{dr:+.2f},{ds:+.2f},{a:.2f},{chi2:.0f},{'yes' if a > 1 else 'no'},{kept}")
