"""Natural wrong readings, each run on the seven units at the bundle noise.

    python wrong_routes.py <units folder from `make_captures.py --out`> <csv out> [route ...]

Routes:
  raw_amplitudes     each direction's recorded amplitude read directly (no runout, no filter), X and Y taken
                     as the two modes;
  probe_reading      each direction's Bode curve compensated (filter undone, slow-roll vector removed) and
                     read directly, X and Y taken as the two modes (probe_reading.py);
  hv_split           every bearing's two probes resolved into a horizontal and a vertical response with the
                     stated probe angles and rotation, as if the supports' principal directions were
                     horizontal and vertical, and each read directly;
  one_mode           the reference's fit held to a single mode, reported for both;
  filter_ignored     the reference with the tracking filter left out.
"""
import csv, importlib.util, json, math, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOP = HERE.parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


REF = load("reference", TOP / "solution" / "bode_reduce.py")
PR = load("probe_reading", HERE / "probe_reading.py")


def write(slot, pairs):
    out = {}
    for name, (peak, af) in zip(("low", "high"), sorted(pairs)):
        out[f"peak_{name}_rpm"] = round(float(peak), 1)
        out[f"af684_{name}"] = round(float(af), 3)
    Path(slot).write_text(json.dumps(out) + "\n")


def raw_amplitudes(folder, slot):
    c = PR.table(folder)
    found = []
    for axis in ("x", "y"):
        both = [PR.api_684(c["speed_rpm"], c[f"{side}_{axis}_um"]) for side in ("de", "nde")]
        found.append(tuple(np.mean(both, axis=0)))
    write(slot, found)


def hv_split(folder, slot):
    sheet = json.loads((Path(folder) / "unit_sheet.json").read_text())
    c = PR.table(folder)
    tau = 1.0 / (2.0 * math.pi * float(sheet["tracking_filter"]["bandwidth_hz"]))
    rho = 1.0 if sheet["rotation"] == "clockwise" else -1.0
    curves = {"h": [], "v": []}
    for side in ("de", "nde"):
        nx, zx = PR.compensated(c, sheet, f"{side}_x", tau)
        ny, zy = PR.compensated(c, sheet, f"{side}_y", tau)
        ox, oy = np.argsort(nx), np.argsort(ny)
        nx, zx, ny, zy = nx[ox], zx[ox], ny[oy], zy[oy]
        common = np.intersect1d(nx, ny)
        zx = np.interp(common, nx, zx.real) + 1j * np.interp(common, nx, zx.imag)
        zy = np.interp(common, ny, zy.real) + 1j * np.interp(common, ny, zy.imag)
        tx = math.radians(sheet["probes"][f"{side}_x"]["angle_deg"])
        ty = math.radians(sheet["probes"][f"{side}_y"]["angle_deg"])
        m = np.array([[math.cos(tx), 1j * rho * math.sin(tx)], [math.cos(ty), 1j * rho * math.sin(ty)]])
        hv = np.linalg.solve(m, np.vstack([zx, zy]))
        curves["h"].append(PR.api_684(common, np.abs(hv[0])))
        curves["v"].append(PR.api_684(common, np.abs(hv[1])))
    write(slot, [tuple(np.mean(v, axis=0)) for v in curves.values()])


def one_mode(folder, slot):
    rec = REF.Record(folder)
    n1, z1 = REF.coarse_start(rec)
    modes, _ = REF.fit(rec, [[n1, z1]])
    grid = np.arange(300.0, 3200.0, 0.05)
    pk = REF.peak_and_af(grid, np.abs(REF.mode(grid, *modes[0])))
    write(slot, [pk, pk])


def filter_ignored(folder, slot):
    saved = REF.filter_tau
    REF.filter_tau = lambda sheet: 0.0
    try:
        REF.main(folder, slot)
    finally:
        REF.filter_tau = saved


ROUTES = {"raw_amplitudes": raw_amplitudes, "probe_reading": lambda f, s: PR.main(f, s), "hv_split": hv_split,
          "one_mode": one_mode, "filter_ignored": filter_ignored}


if __name__ == "__main__":
    base, out = Path(sys.argv[1]), Path(sys.argv[2])
    chosen = sys.argv[3:] or list(ROUTES)
    truth = json.loads((base / "truth.json").read_text())
    table = [("route", "unit", "worst_crit_pct", "worst_af_pct")]
    for route in chosen:
        for unit, true in sorted(truth.items()):
            ROUTES[route](base / unit, base / "_route.json")
            got = json.loads((base / "_route.json").read_text())
            crit = max(abs(got[k] / true[k] - 1) * 100 for k in true if k.startswith("peak_"))
            af = max(abs(got[k] / true[k] - 1) * 100 for k in true if k.startswith("af684"))
            table.append((route, unit, f"{crit:.2f}", f"{af:.2f}"))
            print(*table[-1], sep=",", flush=True)
    with open(out, "w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows(table)
