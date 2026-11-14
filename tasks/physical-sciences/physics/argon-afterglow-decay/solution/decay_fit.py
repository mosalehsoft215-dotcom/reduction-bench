"""Fit both afterglow populations in detector space, then integrate wall loss."""
import json
import os
import sys
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


NAMES = ("diffusion_mm2_ms", "wall_speed_m_s", "bulk_loss_ms_inv",
         "transfer_ms_inv", "cross_loss_m3_s", "wall_fraction")


def read(folder, normalize=True):
    info = json.loads((folder / "cell.json").read_text())
    with np.load(folder / "decays.npz", allow_pickle=False) as data:
        probe, reference, dark = (data[k].astype(float) for k in
                                  ("probe", "reference", "dark"))
    transmission = (probe - dark) / (reference - dark)
    sigma = np.sqrt(probe + transmission**2 * reference
                    + (1 - transmission)**2 * dark) / (reference - dark)
    # chords whose every ray passes outside the wall see no metastables: they give the probe power at acquisition
    reach = np.asarray(info["impact_mm"]) - np.max(np.abs(info["ray_offsets_mm"]))
    outside = reach > info["radial_edges_mm"][-1]
    scale = transmission[..., outside].mean(axis=-1, keepdims=True)
    if not normalize:
        return info, transmission, sigma
    return info, transmission / scale, sigma / scale


def pressure_ratio(info):
    """Manometer readings less the evacuated reading, over the reference pressure."""
    readings = np.asarray(info["gauge_readings_torr"])
    return (readings - info["evacuated_reading_torr"]) / info["reference_pressure_torr"]


def geometry(info):
    edges = np.asarray(info["radial_edges_mm"])
    offsets = np.asarray(info["ray_offsets_mm"])
    chords = np.asarray(info["impact_mm"])[:, None] + offsets
    half = np.sqrt(np.maximum(edges**2 - chords[:, :, None]**2, 0))
    lengths = 2 * np.diff(half, axis=-1)
    return lengths[None, :, :, :] * np.asarray(info["opacity_per_mm"])[:, None, None, None]


def light(model, weights, density):
    depth = np.einsum("lcqr,...lr->...lcq", model, density)
    return np.exp(-depth) @ weights


def initial(model, weights, observed, sigma):
    profiles = np.empty((*observed.shape[:2], model.shape[-1]))
    for shot in range(len(profiles)):
        for line in range(model.shape[0]):
            kernel = model[line]
            path = np.einsum("cqr,q->cr", kernel, weights)
            guess = np.linalg.lstsq(path, -np.log(observed[shot, line]), rcond=None)[0]

            def residual(n):
                return (np.exp(-kernel @ n) @ weights - observed[shot, line]) / sigma[shot, line]

            profiles[shot, line] = least_squares(
                residual, np.maximum(guess, 0.02), bounds=(0, np.inf)).x
    return profiles


def transport(info, parameters):
    edges = np.asarray(info["radial_edges_mm"])
    center = (edges[:-1] + edges[1:]) / 2
    areas = np.diff(edges**2) / 2
    diffusion, speed, bulk, transfer, cross, ternary = parameters
    pressure = pressure_ratio(info)
    diff = diffusion / pressure[:, None] * np.array([1.0, info["diffusion_ratio"]])
    flow = diff[..., None] * edges[1:-1] / np.diff(center)
    wall = speed / (1 + speed * (edges[-1] - center[-1]) / diff)
    backward = info["reverse_transfer_ratio"]
    bulk_at_pressure = bulk * (1 - ternary) * pressure + bulk * ternary * pressure**2

    def derivative(t, n):
        loss = cross * n[:, 0] * n[:, 1]
        swap = (transfer * pressure)[:, None] * (n[:, 1] - backward * n[:, 0])
        dn = -bulk_at_pressure[:, None, None] * n - loss[:, None]
        dn[:, 0] += swap
        dn[:, 1] -= swap
        flux = flow * np.diff(n, axis=-1)
        dn[:, :, :-1] += flux / areas[:-1]
        dn[:, :, 1:] -= flux / areas[1:]
        dn[:, :, -1] -= (edges[-1] / areas[-1]) * wall * n[:, :, -1]
        return dn

    return derivative, edges[-1] * wall, areas


def sampling(info):
    gates = np.asarray(info["gates_ms"])
    nodes, weights = leggauss(8)
    times = gates.mean(axis=1)[:, None] + np.diff(gates, axis=1) * nodes / 2
    return times, weights / 2


def simulate(info, p, n0, time):
    derivative, wall, areas = transport(info, p)
    horizon = info["wall_report_ms"]
    shape = n0.shape

    def ode(t, y):
        n = y[:-1].reshape(shape)
        report = info["report_preparation"]
        wall_rate = np.dot(n[report, :, -1], wall[report])
        return np.r_[derivative(t, n).ravel(), wall_rate]

    sol = solve_ivp(ode, (0, max(time.max(), horizon)), np.r_[n0.ravel(), 0.0],
                    dense_output=True, method="DOP853", rtol=2e-8, atol=1e-9)
    values = sol.sol(time.ravel())[:-1].T.reshape((*time.shape, *shape))
    population = (n0[info["report_preparation"]] * areas).sum()
    return values, float(sol.sol(horizon)[-1] / population)


def fit(info, observed, sigma, free=None, fixed=None):
    optics = geometry(info)
    weight = np.asarray(info["ray_weights"])
    n0 = initial(optics, weight, observed[:, 0], sigma[:, 0])
    times, gate_weights = sampling(info)
    center = np.array([2.1, 0.55, 0.12, 0.36, 0.38, 0.45])
    if free is None:
        free = np.arange(len(center))
    if fixed:
        for i, value in fixed.items():
            center[i] = value

    def unpack(x):
        p = center.copy()
        p[free] = x
        return p

    def residual(x):
        profiles, _ = simulate(info, unpack(x), n0, times)
        signal = light(optics, weight, profiles)
        prediction = np.einsum("gqslc,q->sglc", signal, gate_weights)
        return ((prediction - observed[:, 1:]) / sigma[:, 1:]).ravel()

    bounds = np.array([[0.05, 0.01, 0, 0, 0, 0], [10, 6, 2, 3, 3, 1]])
    opt = least_squares(residual, center[free], bounds=bounds[:, free],
                        x_scale="jac", ftol=1e-8)
    if not opt.success:
        raise RuntimeError(opt.message)
    p = unpack(opt.x)
    _, fraction = simulate(info, p, n0, times)
    return answer(p, fraction)


def answer(p, fraction):
    return dict(zip(NAMES, [*map(float, p[:4]), float(p[4]) * 1e-15, fraction]))


if __name__ == "__main__":
    source, destination = map(Path, sys.argv[1:])
    destination.write_text(json.dumps(fit(*read(source)), allow_nan=False) + "\n",
                           encoding="utf-8")
