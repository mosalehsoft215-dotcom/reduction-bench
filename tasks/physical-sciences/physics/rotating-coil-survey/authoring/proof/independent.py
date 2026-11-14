"""A second reduction of a survey, written apart from the reference (needs SciPy).

    python independent.py <survey folder> <answer.json>

Different choices: harmonics by least squares of the increments on differenced cos/sin bases instead of a DFT;
the shaft twist removed by regressing each segment's quadrupole phase on the rotation direction turn by turn; the
eddy decay fitted with scipy curve_fit and kept only when an F-test says it helps; the stray dipole taken as the
intercept of the dipole against the quadrupole over the low plateaus (up to 600 A); the rising branch found from
the highest current reached before a block.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

NMAX = 15


def design(steps):
    th = np.linspace(0.0, 2 * np.pi, steps + 1)
    cols = []
    for n in range(1, NMAX + 1):
        cols.append(np.diff(np.cos(n * th)))
        cols.append(np.diff(np.sin(n * th)))
    return np.column_stack(cols + [np.ones(steps)])


def harmonics(inc, coil):
    seg = coil["segments"][0]
    R = coil["reference_radius_m"]
    D = design(coil["steps_per_turn"])
    sol = np.linalg.lstsq(D, inc.T, rcond=None)[0]           # Phi = sum a cos + b sin
    out = np.zeros((inc.shape[0], NMAX + 1), complex)
    for n in range(1, NMAX + 1):
        kappa = seg["turns"] * (seg["outer_radius_m"] ** n - seg["inner_radius_m"] ** n) / (n * R ** (n - 1))
        a, b = sol[2 * n - 2], sol[2 * n - 1]
        out[:, n] = (a - 1j * b) / kappa                     # Re[C e^{in theta}] = a cos + b sin
    return out


def decay_level(t, y):
    def flat(t, g):
        return g + 0 * t

    def decay(t, g, a, tau):
        return g * (1 + a * np.exp(-t / tau))
    try:
        p, _ = curve_fit(decay, t, y, p0=(y[-1], 1e-3, 3.0), bounds=([0, -0.05, 0.2], [np.inf, 0.05, 20.0]))
    except RuntimeError:
        return float(np.mean(y))
    s1 = np.sum((y - decay(t, *p)) ** 2)
    s0 = np.sum((y - y.mean()) ** 2)
    F = ((s0 - s1) / 2) / (s1 / max(len(y) - 3, 1)) if s1 > 0 else np.inf
    return float(p[0]) if stats.f.sf(F, 2, max(len(y) - 3, 1)) < 0.01 else float(np.mean(y))


def main(folder, out):
    f = Path(folder)
    coil = json.loads((f / "coil.json").read_text())
    sheet = json.loads((f / "magnet_sheet.json").read_text())
    T = np.genfromtxt(f / "turns.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    raw = np.loadtxt(f / "increments.csv", delimiter=",", skiprows=1)
    H = harmonics(raw[:, 3:], coil)
    R = coil["reference_radius_m"]
    index = {(int(a), int(b), int(c)): i for i, (a, b, c) in enumerate(raw[:, :3])}
    blocks = list(dict.fromkeys(int(b) for b in T["block"]))
    level = {b: float(T["current_set_a"][T["block"] == b][0]) for b in blocks}

    def integrated(b):
        rows = T[T["block"] == b]
        dirs = np.where(rows["direction"].astype(str) == "cw", 1.0, -1.0)
        norm = level[b] / rows["current_readback_a"].astype(float) if level[b] > 0 else np.ones(len(rows))
        total = np.zeros(NMAX + 1, complex)
        series = np.zeros(len(rows))
        for s in (1, 2, 3):
            c = np.array([H[index[(b, int(t), s)]] for t in rows["turn"]]) * norm[:, None]
            series += np.abs(c[:, 2])
            ph = np.unwrap(np.angle(c[:, 2]))
            if len(set(dirs)) == 2:
                A = np.column_stack([np.ones_like(dirs), dirs])
                base, tw = np.linalg.lstsq(A, ph, rcond=None)[0]
            else:
                base, tw = float(np.mean(ph)), 0.0
            n = np.arange(NMAX + 1)
            untwisted = c * np.exp(-1j * n[None, :] * (tw / 2) * dirs[:, None])
            total += untwisted.mean(axis=0)
        if level[b] > 0:
            t = rows["start_after_ramp_s"].astype(float)
            g = decay_level(t, series) if t.min() < 5.0 else float(series.mean())
            total *= g / float(series.mean())
        zero = float(np.mean(rows["encoder_zero_mrad"])) * 1e-3
        return total * np.exp(-1j * np.arange(NMAX + 1) * zero)

    rising = []
    for k, b in enumerate(blocks):
        before = [level[x] for x in blocks[:k]]
        rising.append(b if not before or max(before) < level[b] else None)
    i_op = sheet["operating_current_a"]
    op = [b for b in rising if b is not None and abs(level[b] - i_op) < 0.5][0]
    low = [b for b in rising if b is not None and 0 < level[b] <= 600.0]
    c2 = np.array([integrated(b)[2] for b in low])
    c1 = np.array([integrated(b)[1] for b in low])
    A = np.column_stack([np.ones(len(low)), c2]).astype(complex)
    bg = np.linalg.lstsq(A, c1, rcond=None)[0][0]
    g = integrated(op)
    roll = -0.5 * float(np.angle(g[2]))
    centre = -R * (g[1] - bg) / g[2]
    b6 = (g[6] / g[2] * np.exp(4j * roll)).real * 1e4
    Path(out).write_text(json.dumps({"gl_t": float(abs(g[2]) / R), "b6_units": float(b6),
                                     "axis_x_mm": float(centre.real * 1e3), "axis_y_mm": float(centre.imag * 1e3),
                                     "roll_mrad": roll * 1e3}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
