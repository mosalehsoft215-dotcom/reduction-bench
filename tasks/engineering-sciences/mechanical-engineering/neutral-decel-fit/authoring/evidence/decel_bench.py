"""Route bench of the deceleration task: reference, replay reduction and every plain or partial route.

    python decel_bench.py <bundle> <scratch> [redraws] [--replay] [--drift-replay]

Misses are max over the keys of |got / truth - 1| / limit (A 5 %, C 3 %). Routes on the shipped days:
  pairs       textbook: logged wheel speed as truth, one balance A rT + C q V^2 + m g h s, headings averaged
  v1_flat     the earlier winning route: wheel speed as truth, along-track wind w0 + w1 t fitted, uniform slope
              (sign with heading), temperature and density corrections, all passes in one least-squares fit
  v1_flat_warm     the same plus warm-up A (1 + E e^(-t/tau))
  v1_spline        slope as a linear B-spline in GPS position (5 nodes, sign with heading), no warm-up
  v1_spline_warm   B-spline slope plus warm-up
  mast_along  mast wind resolved along the track, no yaw growth, wheel speed as truth
  yaw_only    mast wind with yaw growth K fitted, wheel speed as truth
  scale_only  wheel speed put on the GPS scale, mast wind along the track, no yaw growth
  scaled      pairs times the mean truth / pairs ratio of the two shown days
  per_kg      shown days' coefficients per kg of test mass (no data)
Writes authoring/evidence/decel_bench.json.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "solution"), str(HERE)]
import decel_fit as ref  # noqa: E402
import gps_rk4_replay as rep  # noqa: E402

BAND = {"a_roll_N": 0.05, "c_air_N_per_kmh2": 0.03}
SHOWN = ("DR-23", "DR-08")


def miss(got, want):
    return round(max(abs(got[k] / want[k] - 1.0) / BAND[k] for k in BAND), 3)


def pair(a, c):
    return {"a_roll_N": float(a), "c_air_N_per_kmh2": float(c)}


class Samples:
    """per-sample arrays of one day, with wheel speed either raw or on the GPS scale"""

    def __init__(self, folder, scaled):
        car = json.load(open(folder / "vehicle_sheet.json"))
        self.mass = float(car["mass_kg"])
        me = self.mass + float(car["rotating_mass_kg"])
        _, mast, wheel, gps = rep.load(str(folder))
        w_np = {k: (np.array([r[0] for r in v]), np.array([r[1] for r in v])) for k, v in wheel.items()}
        g_np = {k: (np.array([r[0] for r in v]), np.array([r[1] for r in v])) for k, v in gps.items()}
        a, b = ref.speed_scale(w_np, g_np) if scaled else (1.0, 0.0)
        bearing = math.radians(float(car["track"]["bearing_deg"]))
        t0 = min(rep.clock(r["start_time"]) for r in mast)
        cols = {k: [] for k in ("force", "V", "h", "rt", "q", "wa", "wc", "clock", "x")}
        for r in mast:
            t, vw = w_np[r["run"]]
            v = vw * (a + b * vw ** 2)
            dvdt = ref.savgol_filter(v / 3.6, 51, 2, deriv=1, delta=0.1)[25:-25]
            vs = ref.savgol_filter(v, 51, 2)[25:-25]
            tg, xg = g_np[r["run"]]
            n = len(vs)
            temp, kpa = float(r["temp_C"]), float(r["press_kPa"])
            to = math.radians(float(r["wind_dir_deg"])) + math.pi - bearing
            w = float(r["wind_speed_ms"]) * 3.6
            cols["force"].append(-me * dvdt)
            cols["V"].append(vs)
            cols["h"].append(np.full(n, 1.0 if r["direction"] == "eastbound" else -1.0))
            cols["rt"].append(np.full(n, 1.0 / (1.0 + ref.K_T * (temp - 20.0))))
            cols["q"].append(np.full(n, kpa * 1000.0 / (ref.R_AIR * (temp + 273.15)) / ref.RHO_REF))
            cols["wa"].append(np.full(n, w * math.cos(to)))
            cols["wc"].append(np.full(n, abs(w * math.sin(to))))
            cols["clock"].append(t[25:-25] + rep.clock(r["start_time"]) - t0)
            cols["x"].append(np.interp(t[25:-25], tg, xg))
        for k, v in cols.items():
            setattr(self, k, np.concatenate(v))
        self.span = float(car["track"]["span_m"])

    def lsq(self, cols):
        coef, *_ = np.linalg.lstsq(np.column_stack(cols), self.force, rcond=None)
        return coef, float(np.sum((np.column_stack(cols) @ coef - self.force) ** 2))

    def hill(self):
        return self.mass * ref.GRAV * self.h

    def pairs(self):
        c, _ = self.lsq([self.rt, self.q * self.V ** 2, self.hill()])
        return pair(c[0], c[1])

    def mast(self, yaw):
        ux = self.V - self.h * self.wa
        drag = self.q * ux * np.abs(ux)
        cols = [self.rt, drag, self.hill()] + ([drag * np.arctan2(self.wc, ux) ** 2] if yaw else [])
        c, _ = self.lsq(cols)
        return pair(c[0], c[1])

    def v1_route(self, nodes, warm):
        if nodes > 1:
            knots = np.linspace(0.0, self.span, nodes)
            hat = np.clip(1.0 - np.abs(self.x[:, None] - knots[None, :]) / (knots[1] - knots[0]), 0.0, None)
            slope_cols = [self.mass * ref.GRAV * self.h * hat[:, j] for j in range(nodes)]
        else:
            slope_cols = [self.hill()]

        def build(p):
            ux = self.V - self.h * (p[0] + p[1] * self.clock / 3600.0)
            cols = [self.rt, self.q * ux * np.abs(ux)] + slope_cols
            if warm:
                cols.append(self.rt * np.exp(-self.clock / math.exp(p[2])))
            return cols

        cost = lambda p: self.lsq(build(p))[1]
        grid = [(w0, w1, lt) for w0 in np.arange(-30, 31, 3.0) for w1 in (-10.0, 0.0, 10.0)
                for lt in (np.log([400, 1200, 3000]) if warm else [math.log(1000.0)])]
        bounds = [(-60.0, 60.0), (-60.0, 60.0), (math.log(200.0), math.log(6000.0))]
        got = minimize(cost, min(grid, key=cost), method="Nelder-Mead", bounds=bounds,
                       options={"maxiter": 800, "xatol": 1e-4, "fatol": 1e-6})
        c, _ = self.lsq(build(got.x))
        return pair(c[0], c[1])


def shipped(root, truth, place, drift_replay):
    out = {}
    for n in truth:
        raw, cal = Samples(place(n), False), Samples(place(n), True)
        r = {"reference": ref.solve(str(place(n))), "pairs": raw.pairs(), "v1_flat": raw.v1_route(1, False), "v1_flat_warm": raw.v1_route(1, True),
             "v1_spline": raw.v1_route(5, False), "v1_spline_warm": raw.v1_route(5, True),
             "mast_along": raw.mast(False), "yaw_only": raw.mast(True), "scale_only": cal.mast(False)}
        if drift_replay:
            r["v1_route_replay"] = rep.reduce(str(place(n)), "drift")
        out[n] = r
    ratio = {k: np.mean([truth[n][k] / out[n]["pairs"][k] for n in SHOWN]) for k in BAND}
    mass = {n: json.load(open(place(n) / "vehicle_sheet.json"))["mass_kg"] for n in truth}
    per = {k: np.mean([truth[n][k] / mass[n] for n in SHOWN]) for k in BAND}
    for n in out:
        out[n]["scaled"] = {k: out[n]["pairs"][k] * ratio[k] for k in BAND}
        out[n]["per_kg"] = {k: per[k] * mass[n] for k in BAND}
    return {n: {route: miss(got, truth[n]) for route, got in r.items()} for n, r in out.items()}


def main():
    root = ROOT
    scratch = Path(sys.argv[2])
    redraws = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else 8
    data = root / "environment" / "data" / "decel_runs"
    truth = {**json.load(open(data / "reference_pairs.json")), **json.load(open(root / "tests" / "held_pairs.json"))}
    place = lambda n: (data if n in SHOWN else root / "tests" / "held_runs") / n
    report = {"limits": BAND, "shipped": shipped(root, truth, place, "--drift-replay" in sys.argv), "redraws": {}}
    for n, r in report["shipped"].items():
        print(n, r, flush=True)
    if "--replay" in sys.argv:
        report["replay_shipped"] = {n: miss(rep.reduce(str(place(n))), truth[n]) for n in truth}
        print("replay shipped", report["replay_shipped"], flush=True)
    for k in range(1, redraws + 1):
        folder = scratch / f"dr5_redraw_{k}"
        subprocess.run([sys.executable, str(root / "authoring" / "provenance" / "make_decel.py"), str(folder),
                        "--redraw", str(k)], check=True)
        t = json.load(open(folder / "truth.json"))
        line = {"reference": {n: miss(ref.solve(str(folder / n)), t[n]) for n in t}}
        if "--replay" in sys.argv:
            line["replay"] = {n: miss(rep.reduce(str(folder / n)), t[n]) for n in t}
        report["redraws"][k] = line
        print("redraw", k, line, flush=True)
    for key in ("reference", "replay"):
        vals = [v for line in report["redraws"].values() for v in line.get(key, {}).values()]
        if vals:
            report[key + "_worst_redraw"] = max(vals)
    with open(HERE / "decel_bench.json", "w", newline=chr(10)) as fh:
        fh.write(json.dumps(report, indent=1) + chr(10))
    print({k: v for k, v in report.items() if k.endswith("worst_redraw")})


if __name__ == "__main__":
    main()
