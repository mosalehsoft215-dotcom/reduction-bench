"""Second, separately written reduction of one occupation (shares no code with solution/gravity_reduce.py).

    python3 second_gravity.py <occupation directory> <result file>

Differences: the profile is the Lagrange parabola through the survey written out explicitly; the drop model is the
exact solution of z'' = g0 + a z (cosh / sinh, linear in z0, v0, g0) with the curvature part added by numerical double
integration of b z(t)^2 over the drop; the wavelength follows each drop's own logged beat; the nominal pressure from the
station elevation; g and the tide factor from one least-squares fit of all drops on [1, logged tide].
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

LIGHT = 299792458.0


def lagrange(hs, rs):
    def value(h, order=0):
        total = 0.0
        for i in range(3):
            j, k = [m for m in range(3) if m != i]
            den = (hs[i] - hs[j]) * (hs[i] - hs[k])
            if order == 0:
                total += rs[i] * (h - hs[j]) * (h - hs[k]) / den
            elif order == 1:
                total += rs[i] * (2 * h - hs[j] - hs[k]) / den
            else:
                total += rs[i] * 2.0 / den
        return total
    return value


def fit_drop(t, z, a, b):
    s = math.sqrt(a)
    ch, sh = np.cosh(s * t), np.sinh(s * t)
    cols = np.stack([ch, sh / s, (ch - 1) / a], 1)
    corr = np.zeros_like(t)
    for _ in range(3):
        p = np.linalg.lstsq(cols, z - corr, rcond=None)[0]
        base = cols @ p
        acc = b * base ** 2                                  # extra acceleration from the curvature
        vel = np.concatenate([[0.0], np.cumsum((acc[1:] + acc[:-1]) / 2 * np.diff(t))])
        corr = np.concatenate([[0.0], np.cumsum((vel[1:] + vel[:-1]) / 2 * np.diff(t))])
    return p[2]


def reduce(folder):
    folder = Path(folder)
    sheet = json.loads((folder / "occupation.json").read_text())
    prof = lagrange(sheet["survey"]["height_m"], sheet["survey"]["reading_ugal"])
    h1 = sheet["first_tag_height_m"]
    a = -prof(h1, 1) * 1e-8
    b = prof(h1, 2) / 2 * 1e-8
    log = {}
    for line in (folder / "drop_log.csv").read_text().splitlines()[1:]:
        s, d, _, td, p, beat = line.split(",")
        log[(int(s), int(d))] = (float(td), float(p), float(beat))
    pn = 1013.25 * (1 - 0.0065 * sheet["station_elevation_m"] / 288.15) ** 5.2559
    obs, tides = [], []
    for path in sorted((folder / "drops").glob("set_*.csv")):
        s = int(path.stem.split("_")[1])
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        for d in np.unique(data[:, 0]).astype(int):
            rows = data[data[:, 0] == d]
            rows = rows[np.argsort(rows[:, 1])]
            td, p, beat = log[(s, d)]
            lam = LIGHT / ((sheet["reference_laser_mhz"] + beat) * 1e6)
            z = rows[:, 1] * sheet["fringes_per_tag"] * lam / 2
            g = fit_drop(rows[:, 2] - z / LIGHT, z, a, b) * 1e8
            obs.append(g - sheet["pressure_admittance_ugal_per_hpa"] * (p - pn))
            tides.append(td)
    # tide factor and g together: least squares of the drop values on [1, logged tide]
    A = np.stack([np.ones(len(tides)), np.array(tides)], 1)
    g_top = float(np.linalg.lstsq(A, np.array(obs), rcond=None)[0][0])
    return {"g_top_ugal": g_top, "g_mark_ugal": g_top + prof(0.0) - prof(h1)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + "\n")
