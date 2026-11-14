"""Second, separately written reduction of a basalt (authoring evidence only).

    python3 cross_reduction.py <basalt folder> <answer file>

Different choices from the reference at every step:
  - tensor: half the difference of each opposite pair of positions, not a least-squares fit;
  - Arai points: components along the specimen's high-temperature direction and along the oven axis, not vector norms;
  - slope: orthogonal regression, not ordinary least squares;
  - alteration: the running sum of the check differences, as a capacity curve by temperature reached, its median
    over the specimens serving the whole flow;
  - a reheating break: a fixed comparison of the whole slope and the slope above 515 C, not a search for the break;
  - field law: a shape-preserving interpolation of the field test, not a tanh;
  - flow value: the mean after dropping the highest and lowest, not the median;
  - tilt: Rodrigues rotation about the strike.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def load(folder, name):
    with open(Path(folder) / name, newline="") as fh:
        return list(csv.DictReader(fh))


def xyz(r):
    return np.array([float(r[k]) for k in ("x_Am2", "y_Am2", "z_Am2")])


def pair_tensor(rows):
    got = {r["position"]: xyz(r) for r in rows}
    cols = [(got["+" + a] - got["-" + a]) / 2.0 for a in "xyz"]
    k = np.column_stack(cols)
    k = (k + k.T) / 2.0
    return k * 3.0 / np.trace(k)


def orthogonal_slope(x, y):
    a, b = x - x.mean(), y - y.mean()
    sxx, syy, sxy = a @ a, b @ b, a @ b
    return (syy - sxx + math.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)


def arai(rows, inverse, use_checks=True, shared=None, low_cut=250.0):
    z_at, i_at, order = {}, {}, []
    for r in rows:
        t = float(r["temperature_C"])
        if r["step"] == "Z":
            z_at[t] = inverse @ xyz(r)
            order.append(("Z", t))
        elif r["step"] == "I":
            i_at[t] = inverse @ xyz(r)
        elif r["step"] == "P":
            order.append(("P", t, inverse @ xyz(r)))
    temps = sorted(z_at)
    gained = {t: i_at[t] - z_at[t] for t in temps}
    axis = gained[temps[-1]] / np.linalg.norm(gained[temps[-1]])
    top = [t for t in temps if 300 <= t <= 530]
    line = np.linalg.svd(np.array([z_at[t] for t in top]) - np.mean([z_at[t] for t in top], axis=0))[2][0]
    line = line if line @ (z_at[top[0]] - z_at[top[-1]]) > 0 else -line
    # running capacity: every check compares a pTRM made earlier with the same step repeated at a higher temperature
    reached, cum = [temps[0]], [0.0]
    top_sum = bottom_sum = 0.0
    total = float(gained[temps[-1]] @ axis)
    hot = None
    for item in order:
        if item[0] == "Z":
            hot = item[1]
        else:
            t = item[1]
            first = float(gained[t] @ axis)
            again = float((item[2] - z_at[hot]) @ axis)
            if first > 0.05 * total:
                reached.append(hot)
                cum.append(cum[-1] + 0.5 * (again / first - 1.0))     # each interval is crossed by two checks
                top_sum += first * first * (again / first - 1.0)
                bottom_sum += first * first
    excess = np.array(cum[1:]) if len(cum) > 1 else np.zeros(1)
    if shared is None:
        use = use_checks and len(cum) > 1 and float(np.mean(np.diff(cum))) > 0.012
        capacity = (lambda t: 1.0 + float(np.interp(t, reached, cum))) if use else (lambda t: 1.0)
    else:
        capacity = lambda t: 1.0 + float(np.interp(t, shared[0], shared[1]))
    pick = [t for t in temps if t >= low_cut]
    nrm = np.array([z_at[t] @ line for t in pick])
    ptrm = np.array([gained[t] @ axis / capacity(t) for t in pick])
    return -orthogonal_slope(ptrm, nrm), float(excess.max()), line, reached, cum, top_sum / max(bottom_sum, 1e-40)


def run(folder):
    sheet = json.loads((Path(folder) / "basalt_sheet.json").read_text())
    cores, steps, ani = load(folder, "specimens.csv"), load(folder, "heating_steps.csv"), load(folder, "atrm_positions.csv")
    of = lambda rows, n: [r for r in rows if r["specimen"] == n]
    lab = sheet["oven"]["field_ut"]
    cool = load(folder, "cooling_rate_test.csv")
    xs, ys = [], []
    for n in sorted({r["specimen"] for r in cool}):
        mine = of(cool, n)
        t = [float(r["cooling_min"]) for r in mine]
        m = [float(r["moment_Am2"]) for r in mine]
        xs += [math.log10(a / t[0]) for a in t[1:]]
        ys += [b / m[0] - 1.0 for b in m[1:]]
    rate = float(np.dot(xs, ys) / np.dot(xs, xs))
    flow = sheet["outcrop"]
    decades = math.log10(flow["natural_cooling_days"] * 24.0 * 60.0 / sheet["oven"]["cooling_min"])
    nature = 1.0 + (rate * decades if rate > 0.005 else 0.0)

    true_lab = lab
    for c in cores:
        if c["solenoid_field_ut"]:
            own = np.linalg.inv(pair_tensor(of(ani, c["specimen"])))
            seen = arai(of(steps, c["specimen"]), own, use_checks=False)[0] * lab
            own_time = math.log10(float(c["solenoid_cooling_min"]) / sheet["oven"]["cooling_min"])
            ratio = float(c["solenoid_field_ut"]) * (1.0 + (rate * own_time if rate > 0.005 else 0.0)) / seen
            true_lab = lab * ratio if abs(ratio - 1.0) > 0.015 else lab

    test = load(folder, "trm_field_test.csv")
    grid = sorted({float(r["field_ut"]) for r in test})
    curve = []
    for n in sorted({r["specimen"] for r in test}):
        m = np.array([float(r["moment_Am2"]) for r in of(test, n)])
        curve.append(m / m[1])
    curve = np.median(curve, axis=0)
    law = PchipInterpolator([0.0] + [b * true_lab / lab for b in grid], [0.0] + list(curve))
    span = np.linspace(0.0, grid[-1] * true_lab / lab * 1.3, 4000)

    rock = []                                    # specimens whose own checks stay moderate, and their tensors
    for c in cores:
        if not c["solenoid_field_ut"]:
            inv = np.linalg.inv(pair_tensor(of(ani, c["specimen"])))
            first = arai(of(steps, c["specimen"]), inv)
            rock.append((c, inv, first))
    typical = float(np.median([r[2][5] for r in rock]))
    rock = [r for r in rock if r[2][5] <= max(0.02, 3.0 * typical)]     # a specimen changing faster than the flow goes
    # one running capacity for the flow: the median of the specimens' running sums on a common temperature scale
    grid = np.array([350.0, 400.0, 430.0, 460.0, 490.0, 510.0, 530.0, 550.0, 565.0, 580.0])
    curves = np.array([np.interp(grid, r[2][3], r[2][4]) for r in rock])
    curve = np.median(curves, axis=0)
    changed = typical > 0.008                                      # checks weighted by their size, over the flow
    shared = (grid, curve) if changed else (grid, np.zeros_like(grid))
    # a break shared by the flow: the slope from 250 to 400 C against the one above 430 C
    lows = [arai(of(steps, r[0]["specimen"]), r[1], shared=shared)[0] for r in rock]
    highs = [arai(of(steps, r[0]["specimen"]), r[1], shared=shared, low_cut=515.0)[0] for r in rock]
    cut = 515.0 if abs(float(np.median(lows)) / float(np.median(highs)) - 1.0) > 0.07 else 250.0
    values, lines = [], []
    for c, inv, _ in rock:
        b, line = (lambda got: (got[0], got[2]))(arai(of(steps, c["specimen"]), inv, shared=shared, low_cut=cut))
        target = b / nature * float(law(true_lab))
        vals = law(span)
        values.append(float(np.interp(target, vals, span)) if rate >= -1 else 0.0)
        az, pl = math.radians(float(c["azimuth_deg"])), math.radians(float(c["plunge_deg"]))
        zc = np.array([math.cos(pl) * math.cos(az), math.cos(pl) * math.sin(az), math.sin(pl)])
        xc = np.array([math.sin(pl) * math.cos(az), math.sin(pl) * math.sin(az), -math.cos(pl)])
        yc = np.cross(zc, xc)
        lines.append(line[0] * xc + line[1] * yc + line[2] * zc)
    values = sorted(values)
    field = float(np.mean(values[1:-1]))
    mean = np.mean(lines, axis=0)
    strike, dip = math.radians(flow["bedding_strike_deg"]), math.radians(flow["bedding_dip_deg"])
    u = np.array([math.cos(strike), math.sin(strike), 0.0])
    a = -dip
    mean = mean * math.cos(a) + np.cross(u, mean) * math.sin(a) + u * (u @ mean) * (1 - math.cos(a))
    inc = math.atan2(mean[2], math.hypot(mean[0], mean[1]))
    colat = math.pi / 2 - math.atan(math.tan(inc) / 2)
    return {"paleofield_ut": field,
            "vdm_zam2": 4 * math.pi * 6.371e6 ** 3 / (4e-7 * math.pi) * field * 1e-6 / math.sqrt(1 + 3 * math.cos(colat) ** 2) / 1e21}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(run(sys.argv[1])))
