"""Second, separately written reduction of one hot-wire cell (no code shared with solution/hotwire_reduce.py).

    python3 second_hotwire.py <cell directory> <result file>

Differences: the spool factor comes from the short wire's pre-heat reading; the section temperature by Newton
iteration on the certificate polynomial; the kernel is the exact line-source exponential integral E1 with the
wire-capacity (Healy) term added; one global least-squares fit over all three runs at once (Levenberg-Marquardt on
log lambda and log rho c_p).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import exp1

GAMMA_E = 0.5772156649015329


def load(folder):
    cell = json.loads((folder / "cell.json").read_text())
    rows = [dict(zip(*[l.split(",") for l in (lines[0], line)])) for lines in [(folder / "run_log.csv").read_text().split()]
            for line in lines[1:]]
    return cell, rows


def section_rise(per_len, factor, cert, bath):
    a, b, r0 = cert["A_per_k"], cert["B_per_k2"], cert["r0_ohm_per_m_at_0c"]
    target = per_len / (factor * r0) - 1
    d = np.full_like(per_len, bath - 273.15)
    for _ in range(30):
        d = d - (a * d + b * d * d - target) / (a + 2 * b * d)
    return d - (bath - 273.15)


def reduce(folder):
    folder = Path(folder)
    cell, rows = load(folder)
    cert = cell["wire_certificate"]
    span = cell["long_wire_length_m"] - cell["short_wire_length_m"]
    radius = cell["wire_radius_um"] * 1e-6
    runs = []
    for row in rows:
        bath = float(row["bath_k"])
        db = bath - 273.15
        factor = float(row["r0_short_ohm"]) / (cell["short_wire_length_m"] * cert["r0_ohm_per_m_at_0c"]
                                                * (1 + cert["A_per_k"] * db + cert["B_per_k2"] * db * db))
        amps = float(row["v_std_mv"]) / 1000.0 / cell["standard_resistor_ohm"]
        data = np.loadtxt(folder / "runs" / f"run_{row['run']}.csv", delimiter=",", skiprows=1)
        t, per_len = data[:, 0], (data[:, 1] - data[:, 2]) / span
        q = amps * amps * per_len
        steps = np.diff(np.concatenate([[0.0], q]))
        lag = t[:, None] - np.concatenate([[0.0], t[:-1]])[None, :]
        keep = t >= 0.05
        runs.append((lag[keep], steps, section_rise(per_len, factor, cert, bath)[keep]))

    def residual(p):
        lam, rc = math.exp(p[0]), math.exp(p[1])
        diff = lam / rc
        out = []
        for lag, steps, rise in runs:
            pos = lag > 0
            u = np.where(pos, radius * radius / (4 * diff * np.where(pos, lag, 1.0)), 1.0)
            ln_term = np.log(1 / u) - GAMMA_E
            kern = np.where(pos, exp1(u) + 2 * u * ((1 - 2.85e6 / rc) * ln_term + 1), 0.0)
            out.append(kern @ steps / (4 * math.pi * lam) - rise)
        return np.concatenate(out)

    fit = least_squares(residual, [math.log(0.1), math.log(1.5e6)], method="lm")
    return {"conductivity_mw_per_m_k": math.exp(fit.x[0]) * 1e3,
            "volumetric_heat_capacity_mj_per_m3_k": math.exp(fit.x[1]) / 1e6}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + "\n")
