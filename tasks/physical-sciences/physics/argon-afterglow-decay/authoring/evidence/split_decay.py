"""Independent reduction using exact reaction substeps and matrix exponentials."""
import json
import os
import sys
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.linalg import expm
from scipy.optimize import least_squares, nnls
from scipy.special import logsumexp


def assemble(meta):
    r = np.array(meta["radial_edges_mm"])
    volume = np.diff(r * r) / 2
    centers = 0.5 * (r[1:] + r[:-1])
    h = r[-1] - centers[-1]
    face = r[1:-1] / np.diff(centers)
    count = len(centers)
    jump = np.diff(np.eye(count), axis=0)
    lap = -(jump.T @ (face[:, None] * jump)) / volume[:, None]
    return lap, volume, h, r[-1]


def optics(meta):
    r = np.array(meta["radial_edges_mm"])
    y = np.add.outer(meta["impact_mm"], meta["ray_offsets_mm"])
    cap = np.sqrt(np.clip(np.subtract.outer(r * r, y * y), 0, None))
    path = 2 * np.diff(cap, axis=0).transpose(1, 2, 0)
    return np.array(meta["opacity_per_mm"])[:, None, None, None] * path


def reduce_case(folder):
    meta = json.loads((folder / "cell.json").read_text())
    with np.load(folder / "decays.npz", allow_pickle=False) as data:
        blank = data["reference"] - data["dark"]
        target = (data["probe"] - data["dark"]) / blank
    ray_y = np.add.outer(meta["impact_mm"], meta["ray_offsets_mm"])
    clear = np.all(np.abs(ray_y) > meta["radial_edges_mm"][-1], axis=1)
    target = target / np.median(target[..., clear], axis=-1)[..., None]
    rays = np.array(meta["ray_weights"])
    kernel = optics(meta)
    lap, volume, h, wall_radius = assemble(meta)
    nring = len(volume)
    reverse = meta["reverse_transfer_ratio"]
    report = meta["report_preparation"]
    pressure = (np.array(meta["gauge_readings_torr"]) - meta["evacuated_reading_torr"]) / meta["reference_pressure_torr"]
    end = meta["wall_report_ms"]
    gates = np.asarray(meta["gates_ms"])
    ticks = np.linspace(0, max(gates.max(), end), 1 + int(np.ceil(max(gates.max(), end) / 0.035)))
    dt = ticks[1]
    fractions = np.linspace(0, 1, 13)
    window = gates[:, :1] + (gates[:, 1:] - gates[:, :1]) * fractions

    def transmitted(n):
        return np.exp(-np.einsum("lcqr,...lr->...lcq", kernel, n)) @ rays

    initial = np.empty((target.shape[0], 2, nring))
    for s in range(target.shape[0]):
        for j, paths in enumerate(kernel):
            effective = rays @ paths
            target_depth = -np.log(target[s, 0, j])
            profile = nnls(effective, target_depth)[0]
            log_weights = np.log(rays)

            def optical_depth(n):
                return -logsumexp(log_weights - paths @ n, axis=1)

            def slope(n):
                log_flux = log_weights - paths @ n
                fractions = np.exp(log_flux - logsumexp(log_flux, axis=1, keepdims=True))
                return np.einsum("cq,cqr->cr", fractions, paths)

            initial[s, j] = least_squares(lambda n: optical_depth(n) - target_depth,
                                          np.maximum(profile, 0.01), jac=slope,
                                          bounds=(0, np.inf), x_scale="jac").x

    def forward(p):
        d, speed, bulk, forward_rate, cross, ternary = p
        ds = d / pressure[:, None] * np.array([1.0, meta["diffusion_ratio"]])
        boundary = wall_radius / (1 / speed + h / ds)
        swap = np.kron(np.array([[-reverse, 1], [reverse, -1]]) * forward_rate,
                       np.eye(nring))
        steps = {}
        for s, density_ratio in enumerate(pressure):
            if density_ratio in steps:
                continue
            op = np.kron(np.diag(ds[s]), lap) + swap * density_ratio
            op -= np.eye(2 * nring) * bulk * ((1 - ternary) * density_ratio + ternary * density_ratio**2)
            for j in range(2):
                i = (j + 1) * nring - 1
                op[i, i] -= boundary[s, j] / volume[-1]
            steps[density_ratio] = expm(op * (dt / 2))
        linear = np.stack([steps[q] for q in pressure])
        profiles = np.empty((len(ticks), *initial.shape))
        profiles[0] = initial
        for k in range(1, len(ticks)):
            inter = np.einsum("sij,sj->si", linear, profiles[k - 1].reshape(-1, 2 * nring)).reshape(initial.shape)
            high, low = inter.max(axis=1), inter.min(axis=1)
            gap = high - low
            phi = np.full_like(gap, cross * dt)
            np.divide(-np.expm1(-gap * cross * dt), gap, out=phi, where=gap > 0)
            consumed = high * low * phi / (1 + low * phi)
            inter -= consumed[:, None, :]
            profiles[k] = np.einsum("sij,sj->si", linear, inter.reshape(-1, 2 * nring)).reshape(initial.shape)
        sampled = CubicSpline(ticks, profiles)(window)
        signal = simpson(transmitted(sampled), x=fractions, axis=1).transpose(1, 0, 2, 3)
        rate = profiles[:, report, :, -1] @ boundary[report]
        fraction = CubicSpline(ticks, rate).integrate(0, end) / (initial[report] * volume).sum()
        return signal, float(fraction)

    def residual(p):
        predicted, _ = forward(p)
        return ((predicted - target[:, 1:]) * np.sqrt(blank[:, 1:] / (1 + target[:, 1:]))).ravel()

    result = least_squares(residual, [3.2, 0.95, 0.22, 0.6, 0.25, 0.55],
                           bounds=([0.05, 0.01, 0, 0, 0, 0], [10, 6, 2, 3, 3, 1]),
                           x_scale="jac", ftol=1e-8)
    d, v, b, t, c, _ = result.x.tolist()
    return dict(diffusion_mm2_ms=d, wall_speed_m_s=v, bulk_loss_ms_inv=b,
                transfer_ms_inv=t, cross_loss_m3_s=c * 1e-15,
                wall_fraction=forward(result.x)[1])


if __name__ == "__main__":
    source, destination = map(Path, sys.argv[1:])
    destination.write_text(json.dumps(reduce_case(source), allow_nan=False), encoding="utf-8")
