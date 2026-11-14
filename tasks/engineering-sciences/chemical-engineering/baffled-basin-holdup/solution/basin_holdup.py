"""Holdup volume and spread ratio of one baffled basin from a salt-dose trial.

usage: python3 basin_holdup.py <trial folder> <report file>

Axis: inflow volume u passed since the dose (trapezoid of the logged inflow); the basin's internal exchange follows
the inflow, so its residence distribution E is fixed on that axis.  E belongs to a line plug Vp, n equal stirred
compartments holding Vm and stagnant pockets Vs trading with the last compartment at x times the inflow.  The salt
term is pinned to the dosed mass times the cell factor, because every gram leaves eventually even when the record
stops early.  The inflow's own conductivity is carried through the same basin, so each change of it appears as a
scaled copy of F (running integral of E): one change before the dose (w1 m3 ahead of it) and one after (w2 m3).
    reading = k0 + d1 F(u + w1) + d2 F(u - w2) + c M E(u)
The dose itself is a short box at the inflow of the moment.  Stage 1 fits the early record (n = 1..4, no change
after the dose); stage 2 fits the whole record for the best two or three n, the late change started at four places.
k0, d1, d2 are solved linearly at every step of the bounded least squares.  Reported: the first moment of E and variance over mean squared,
from the cumulants of the transfer function.
"""
import json
import sys

import numpy as np
from scipy.optimize import least_squares

POINTS = 1 << 14


def read_trial(folder):
    with open(f"{folder}/trial_card.json") as fh:
        card = json.load(fh)
    table = np.loadtxt(f"{folder}/outlet_trace.csv", delimiter=",", skiprows=1)
    return card, table[:, 0], table[:, 1], table[:, 2] * 3.6


def g_of_s(s, vp, vm, n, vs, x):
    vc = vm / n
    return np.exp(-vp * s) / (1.0 + vc * s) ** (n - 1) / (vc * s + 1.0 + x * vs * s / (x + vs * s))


def moments(vp, vm, n, vs, x):
    h = 1e-4 / (vp + vm + vs)
    k = [np.log(g_of_s(z, vp, vm, n, vs, x).real) for z in (-h, 0.0, h)]
    mean = -(k[2] - k[0]) / (2 * h)
    return mean, (k[2] - 2 * k[1] + k[0]) / h ** 2


def build_residual(u, k, salt, scale, box, keep, n, late):
    """misfit function of the log-volumes, log-exchange and the two change positions (in units of scale)"""
    du = 50.0 * scale / POINTS
    freq = np.fft.rfftfreq(POINTS, du)
    spin, smear, axis = 2j * np.pi * freq, np.sinc(freq * box), np.arange(POINTS) * du
    uu, kk = u[keep], k[keep]

    def misfit(z):
        vp, vm, vs, x = np.exp(z[:4])
        g = g_of_s(spin, vp, vm, n, vs, x)
        pulse = salt * np.interp(uu, axis, np.fft.irfft(g * smear, POINTS) / du, left=0.0, right=0.0)
        bare = np.fft.irfft(g, POINTS) / du
        run = np.r_[0.0, np.cumsum(0.5 * (bare[1:] + bare[:-1]) * du)]
        basis = [np.ones_like(uu), np.interp(uu + z[4] * scale, axis, run, left=0.0, right=1.0)]
        if late:
            basis.append(np.interp(uu - z[5] * scale, axis, run, left=0.0, right=1.0))
        basis = np.column_stack(basis)
        rest = kk - pulse
        return basis @ np.linalg.lstsq(basis, rest, rcond=None)[0] - rest

    return misfit


def bounded_fit(misfit, z0, scale, late):
    lower = np.r_[np.log([1e-3 * scale, 0.05 * scale, 1e-3 * scale, 0.01]), 0.0, 0.05]
    upper = np.r_[np.log([0.5 * scale, 2.0 * scale, 2.0 * scale, 5.0]), 2.5, 8.0]
    width = 6 if late else 5
    lower, upper = lower[:width], upper[:width]
    sol = least_squares(misfit, np.clip(np.asarray(z0)[:width], lower + 1e-9, upper - 1e-9), bounds=(lower, upper),
                        x_scale="jac", max_nfev=300)
    return sol.cost, sol.x


def guess(scale, vp, vm, vs, x, w1, w2):
    return np.r_[np.log([vp * scale, vm * scale, vs * scale, x]), w1, w2]


def evaluate(folder):
    card, t, k, q = read_trial(folder)
    t_dose = card["salt_dose"]["mid_time_s"]
    u = np.r_[0.0, np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))] / 3600.0
    u -= np.interp(t_dose, t, u)
    box = np.interp(t_dose, t, q) / 3600.0 * card["salt_dose"]["duration_s"]
    salt = card["conductivity_cell"]["uS_cm_per_g_m3"] * 1000.0 * card["salt_dose"]["kg"]
    scale = max(salt / max(float(np.max(k - np.median(k[u < 0]))), 1e-9) * 0.6, 0.5)
    # stage 1: record up to 1.6 scale volumes, only the change of the inflow conductivity ahead of the dose
    head = u < 1.6 * scale
    stage_one = []
    for n in (1, 2, 3, 4):
        misfit = build_residual(u, k, salt, scale, box, head, n, False)
        for z0 in (guess(scale, 0.05, 0.6, 0.4, 0.2, 0.4, 1.5), guess(scale, 0.05, 0.5, 0.8, 0.1, 0.8, 1.5)):
            cost, z = bounded_fit(misfit, z0, scale, False)
            stage_one.append((cost, n, z))
    stage_one.sort(key=lambda item: item[0])
    # stage 2: whole record with a second change after the dose, started across the late record
    whole = np.ones(u.size, bool)
    chosen = (np.inf, None, None)
    tried = []
    for _, n, z in stage_one[:3]:
        if n in tried:
            continue
        tried.append(n)
        misfit = build_residual(u, k, salt, scale, box, whole, n, True)
        for w2 in np.linspace(1.4, max(u[-1] / scale - 0.3, 1.5), 4):
            cost, zz = bounded_fit(misfit, np.r_[z, w2], scale, True)
            if cost < chosen[0]:
                chosen = (cost, n, zz)
    _, n, z = chosen
    vp, vm, vs, x = np.exp(z[:4])
    mean, var = moments(vp, vm, n, vs, x)
    return float(mean), float(var / mean ** 2)


def main(folder, target):
    hold, spread = evaluate(folder)
    with open(target, "w") as fh:
        json.dump({"holdup_volume_m3": round(hold, 4), "spread_ratio": round(spread, 5)}, fh)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
