#!/usr/bin/env python3
"""A model-free reading of each probe direction (the misreading `probe_reading` of wrong_routes.py): API 684 read
directly off each direction's compensated Bode curve, the X and Y results then taken as the two modes.

    python3 probe_reading.py <unit folder> <answer.json>

  - filter: a single-pole tracking filter of time constant tau = 1/(2 pi B) gives y' = (x - y)/tau, so the
    vector before the filter is x = y + tau dy/dt; dy/dt comes from a local cubic in time (7 rows);
  - runout: the slow-roll vector the sheet gives for each probe (barring run, 300 rpm), subtracted
    vectorially;
  - curve: the compensated amplitude is smoothed by a local quadratic in speed (+-40 rpm), sampled every
    0.5 rpm, and read for the peak and the two half-power speeds; DE and NDE averaged per direction.
With unequal supports each probe sees a blend of both support modes, so these are not the modes' own values.
"""
import csv, json, math, sys
from pathlib import Path

import numpy as np


def table(folder):
    with open(Path(folder) / "trip.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def local_quadratic(x, y, at, half):
    out = np.empty(at.size, dtype=y.dtype)
    for i, x0 in enumerate(at):
        m = np.abs(x - x0) <= half
        if m.sum() < 5:
            m = np.argsort(np.abs(x - x0))[:5]
        dx = x[m] - x0
        A = np.column_stack([np.ones(dx.size), dx, dx * dx])
        out[i] = np.linalg.lstsq(A, y[m], rcond=None)[0][0]
    return out


def slope_in_time(t, z):
    d = np.empty_like(z)
    for i in range(t.size):
        lo, hi = max(0, i - 3), min(t.size, i + 4)
        dt = t[lo:hi] - t[i]
        A = np.column_stack([np.ones(dt.size), dt, dt * dt, dt ** 3])
        d[i] = np.linalg.lstsq(A, z[lo:hi], rcond=None)[0][1]
    return d


def api_684(n, amp):
    grid = np.arange(n.min() + 40, n.max() - 40, 0.5)
    a = local_quadratic(n, amp, grid, 40.0)
    k = int(np.argmax(a))
    half = a[k] / math.sqrt(2)
    i = k
    while a[i] > half:
        i -= 1
    j = k
    while a[j] > half:
        j += 1
    n1 = grid[i] + (half - a[i]) / (a[i + 1] - a[i]) * 0.5
    n2 = grid[j] - (half - a[j]) / (a[j - 1] - a[j]) * 0.5
    return grid[k], grid[k] / (n2 - n1)


def compensated(c, sheet, probe, tau):
    """Speeds and compensated amplitudes of one probe: filter undone, slow-roll vector removed, bad gaps out."""
    gap_lo, gap_hi = sheet["gap_healthy_v"]
    t, n = c["time_s"], c["speed_rpm"]
    falling = np.r_[True, np.diff(n) < -1.0]           # rows of a steady hold add nothing to the curve
    z = c[f"{probe}_um"] * np.exp(-1j * np.radians(c[f"{probe}_lag_deg"]))
    if tau > 0:
        z = z + tau * slope_in_time(t, z)
    gap = c[f"{probe}_gap_v"]
    ok = falling & (gap >= gap_lo) & (gap <= gap_hi)
    sr = sheet["slow_roll"]["vectors"][probe]
    return n[ok], z[ok] - sr["um"] * np.exp(-1j * np.radians(sr["lag_deg"]))


def main(folder, answer):
    sheet = json.loads((Path(folder) / "unit_sheet.json").read_text())
    c = table(folder)
    bw = float(sheet["tracking_filter"]["bandwidth_hz"])
    tau = 1.0 / (2.0 * math.pi * bw) if bw > 0 else 0.0
    found = []
    for axis in ("x", "y"):
        both = []
        for probe in (f"de_{axis}", f"nde_{axis}"):
            n_ok, z_ok = compensated(c, sheet, probe, tau)
            both.append(api_684(n_ok, np.abs(z_ok)))
        found.append(tuple(np.mean(both, axis=0)))
    out = {}
    for name, (peak, af) in zip(("low", "high"), sorted(found)):
        out[f"peak_{name}_rpm"] = round(float(peak), 1)
        out[f"af684_{name}"] = round(float(af), 3)
    Path(answer).write_text(json.dumps(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
