"""Recover two complex shear moduli from stepped-coupon transmission."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import j0, j1


def calibrated_response(records, label, frequency, use_pilot):
    data = records[label]
    spectrum = np.fft.rfft(data, axis=2)
    if use_pilot:
        zero = records[f"{label}_zero_us"]
        polarity = records[f"{label}_polarity"]
        correction = polarity[:, :, None] * np.exp(2j * np.pi * frequency[None, None, :, None]
                                                  * zero[:, :, None])
        spectrum *= correction[..., None]
    return spectrum.mean(axis=1)


def modal_estimates(response, deltas, frequency, slowness_grid, base=0.0, radius=0.0):
    """Use a well-separated ratio to establish the common receiver-space modes."""
    ratios = np.linalg.solve(response[0].T, response[1:].transpose(0, 2, 1)).transpose(0, 2, 1)
    eigvals, eigvecs = np.linalg.eig(ratios)
    separation = np.abs(eigvals[:, 0] - eigvals[:, 1])
    basis = eigvecs[np.argmax(separation)]
    multipliers = np.linalg.solve(basis, ratios @ basis).diagonal(axis1=1, axis2=2)
    expected = np.exp(-2j * np.pi * frequency * slowness_grid[:, None] * deltas)
    measured = multipliers / np.maximum(abs(multipliers), 1e-15)
    best = slowness_grid[np.argmax((expected.conj() @ measured).real, axis=0)]
    raw = multipliers
    for _ in range(3 if radius else 1):
        cycles = -np.angle(multipliers) / (2 * np.pi)
        cycles += np.rint(frequency * deltas[:, None] * best - cycles)
        slowness = deltas @ cycles / (frequency * (deltas @ deltas))
        attenuation = -deltas @ np.log(np.maximum(abs(multipliers), 1e-15)) / (deltas @ deltas)
        if radius:            # divide out the beam spreading of each mode between the thin region and region j
            best = slowness
            ratio = (spreading(frequency, slowness[None, :], (base + deltas)[:, None], radius)
                     / spreading(frequency, slowness[None, :], np.array([[base]]), radius))
            multipliers = raw / ratio
    return np.column_stack((slowness, attenuation)), basis


def placement_eigenvalues(direct, echo, basis):
    """Within one placement the triple transit and the direct pass share the couplant films, so
    echo x direct^-1 holds only twice the material path of that region."""
    inside = np.linalg.solve(direct.transpose(0, 2, 1), echo.transpose(0, 2, 1)).transpose(0, 2, 1)
    return np.linalg.solve(basis[None], inside @ basis[None]).diagonal(axis1=1, axis2=2)


def film_free_fit(lam, frequency, d, radius, f_ref, start):
    """Fit s = s0 + slope (f/f_ref - 1) and alpha = a_ref f/f_ref of each mode to the region-to-region
    change of the in-placement eigenvalues lam[f, region, mode]."""
    x = frequency / f_ref
    paths = 2 * (d[1:] - d[0])
    result = []
    for m in range(lam.shape[2]):
        s0, slope = start[m]
        for it in range(4):
            s = s0 + slope * (x - 1)
            corrected = lam[:, :, m]
            if radius:
                corrected = corrected / (spreading(frequency[:, None], s[:, None], 3 * d[None], radius)
                                         / spreading(frequency[:, None], s[:, None], d[None], radius))
            ratio = corrected[:, 1:] / corrected[:, :1]
            weight = abs(corrected[:, 1:] * corrected[:, :1]) / (abs(corrected[:, 1:])**2 + abs(corrected[:, :1])**2)
            if it == 0:
                unit = weight * ratio / np.maximum(abs(ratio), 1e-15)
                a_grid = np.arange(s0 - .06, s0 + .06, .0002)
                b_grid = np.arange(-.10, .1001, .002)
                trial = a_grid[:, None, None] + b_grid[None, :, None] * (x - 1)[None, None]
                score = np.zeros(trial.shape[:2])
                for r, path in enumerate(paths):
                    score += (np.exp(2j * np.pi * frequency * path * trial) * unit[None, None, :, r]).real.sum(-1)
                i, k = np.unravel_index(np.argmax(score), score.shape)
                s0, slope = a_grid[i], b_grid[k]
                s = s0 + slope * (x - 1)
            cycles = -np.angle(ratio) / (2 * np.pi)
            cycles += np.rint(frequency[:, None] * paths[None] * s[:, None] - cycles)
            per_f = (cycles * paths).sum(1) / (frequency * (paths @ paths))
            w = weight.sum(1)
            design = np.column_stack((np.ones_like(x), x - 1)) * np.sqrt(w)[:, None]
            s0, slope = np.linalg.lstsq(design, per_f * np.sqrt(w), rcond=None)[0]
        loss = -np.log(np.maximum(abs(ratio), 1e-15))
        a_ref = (weight * loss * x[:, None] * paths).sum() / (weight * (x[:, None] * paths)**2).sum()
        result.append((s0, a_ref))
    return result


def gated(response, start, end, dt):
    n = 2 * (response.shape[1] - 1)
    time = np.arange(n) * dt
    window = .5 * (np.tanh((time - start[:, None]) / .7) - np.tanh((time - end[:, None]) / .7))
    return np.fft.rfft(np.fft.irfft(response, n=n, axis=1) * window[..., None, None], axis=1)


def spreading(f, s, d, a):
    """Lommel factor of a coaxial piston pair (radius a mm) over d mm, in the exp(i w t - i k d) convention."""
    x = 2 * np.pi * f * s * a * a / d
    return np.conj(1 - np.exp(-1j * x) * (j0(x) + 1j * j1(x)))


def measure(folder, *, gate=True, subtract=True, use_pilot=True, longest_only=False, track=True, beam=True,
            films=True):
    geometry = json.loads((folder / "geometry.json").read_text())
    with np.load(folder / "waveforms.npz", allow_pickle=False) as arrays:
        n = arrays["shot"].shape[2]
        dt = geometry["sample_interval_us"]
        frequency = np.fft.rfftfreq(n, dt)
        response = calibrated_response(arrays, "shot", frequency, use_pilot)
        if subtract:
            response -= calibrated_response(arrays, "lift_off", frequency, use_pilot)

    start, end = np.asarray(geometry["direct_gate_us"]).T
    full = response
    if gate:
        # A smooth shoulder avoids spectral leakage at the supplied gate edges.
        response = gated(full, start, end, dt)

    d = np.asarray(geometry["thickness_mm"])
    delta = d[1:] - d[0]
    lo, hi = geometry["useful_band_mhz"]
    idx = np.flatnonzero((frequency >= lo) & (frequency <= hi))[::2]
    grid = np.linspace(*geometry["phase_slowness_bounds_us_per_mm"], 1601)
    if longest_only:
        # Linear phase of the longest baseline alone has multiple compatible branches.
        response = response[[0, -1]]
        delta = delta[-1:]
    estimates, bases = [], []
    previous = None
    for j in idx:
        modal, basis = modal_estimates(response[:, j], delta, frequency[j], grid, d[0],
                                       geometry.get("transducer_radius_mm", 0.0) if beam else 0.0)
        if track and previous is not None:
            similarity = abs(previous.conj().T @ basis)**2
            swap = similarity[0, 1] + similarity[1, 0] > np.trace(similarity)
        else:
            swap = modal[0, 0] > modal[1, 0]
        if swap:
            modal = modal[::-1]
            basis = basis[:, ::-1]
        previous = basis
        estimates.append(modal)
        bases.append(basis)
    estimates = np.asarray(estimates)
    radius = geometry.get("transducer_radius_mm", 0.0) if beam else 0.0
    if films and not longest_only:
        # Place the triple-transit gate from the direct-pass slownesses, then measure inside each placement.
        s_lo, s_hi = estimates[:, :, 0].min(), estimates[:, :, 0].max()
        echo = gated(full, start + 2 * d * s_lo * .94, end + 2 * d * s_hi * 1.06, dt)
        every = np.flatnonzero((frequency >= lo) & (frequency <= hi))
        lam = np.asarray([placement_eigenvalues(response[:, j], echo[:, j], bases[np.argmin(abs(idx - j))])
                          for j in every])
        f_ref = geometry["reference_frequency_mhz"]
        first = [np.polyfit(frequency[idx] / f_ref - 1, estimates[:, m, 0], 1)[::-1] for m in range(2)]
        values = film_free_fit(lam, frequency[every], d, radius, f_ref, first)
        return report(geometry, values)
    f_ref = geometry["reference_frequency_mhz"]
    centered = frequency[idx] / f_ref - 1
    values = []
    for mode in estimates.transpose(1, 0, 2):
        s, attenuation = mode.T
        design = np.column_stack((np.ones_like(centered), centered))
        coefficients = np.linalg.lstsq(design, s, rcond=None)[0]
        residual = s - design @ coefficients
        good = np.abs(residual - np.median(residual)) < .006
        if good.sum() >= len(good) // 2:
            coefficients = np.linalg.lstsq(design[good], s[good], rcond=None)[0]
        values.append((coefficients[0], np.median(attenuation / (1 + centered))))
    return report(geometry, values)


def report(geometry, values):
    f_ref = geometry["reference_frequency_mhz"]
    moduli = [(geometry["density_kg_m3"] / 1e6) *
              (1000 / complex(s, -a / (2 * np.pi * f_ref)))**2 for s, a in values]
    return {f"{mode}_{component}_mpa": float(val)
            for mode, modulus in zip(("fast", "slow"), [moduli[i] for i in np.argsort(np.asarray(values)[:, 0])])
            for component, val in (("storage", modulus.real), ("loss", modulus.imag))}


if __name__ == "__main__":
    result = measure(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result, allow_nan=False) + "\n", encoding="utf-8")
