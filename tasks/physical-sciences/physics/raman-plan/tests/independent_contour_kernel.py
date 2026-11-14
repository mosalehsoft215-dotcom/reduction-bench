"""Independent contour inversion for three-state Markov-modulated counts.

For row generators, the matrix count generating function is

    G(z) = exp(t * (Q + (z - 1) * diag(photon_rates))).

This module integrates G(z) on r=exp(-8/N), at half-step angles, using full
three-by-three eigendecompositions. It is independent of a unit-circle FFT or
an analytic triangular-block exponential. The sparse overflow implementation
in independent_count_kernel.py remains a separate small-case check.

Let B = P(Poisson(max_rate*t) >= N). Half-step aliasing adds (-r**N)**j
times coefficient n+jN to coefficient n. Therefore the exact-arithmetic K0
alias error is at most r**N * B per entry; the direct finite-sum K1 error is
at most (1+r**N) * B, including the omitted counts >= N. These bounds do not
cover cancellation/eigendecomposition roundoff. Direct tail integration
amplifies numerical error by up to about exp(8), so precision is absolute,
not relative in tiny probabilities. Raw closure and negative residuals are
reported before projection, and a requested 2N run checks quadrature stability.
"""

from __future__ import annotations

import numbers

import numpy as np
from scipy.linalg import expm
from scipy.stats import poisson


def _matrix_exponential_batched(matrices, condition_limit=1e5, chunk_size=8192):
    """Full complex eigendecomposition, with robust defective-matrix fallback."""
    matrices = np.asarray(matrices, dtype=complex)
    original_shape = matrices.shape
    flat = matrices.reshape(-1, 3, 3)
    result = np.empty_like(flat)
    fallback_count = 0
    for start in range(0, len(flat), chunk_size):
        stop = min(start + chunk_size, len(flat))
        block = flat[start:stop]
        try:
            values, vectors = np.linalg.eig(block)
            condition = np.linalg.cond(vectors)
            good = np.isfinite(condition) & (condition <= condition_limit)
            block_result = np.empty_like(block)
            if np.any(good):
                scaled_vectors = vectors[good] * np.exp(values[good])[:, None, :]
                block_result[good] = scaled_vectors @ np.linalg.inv(vectors[good])
                good[good] &= np.all(np.isfinite(block_result[good]), axis=(-2, -1))
        except np.linalg.LinAlgError:
            good = np.zeros(len(block), dtype=bool)
            block_result = np.empty_like(block)
        for index in np.flatnonzero(~good):
            block_result[index] = expm(block[index])
            fallback_count += 1
        result[start:stop] = block_result
    return result.reshape(original_shape), fallback_count


def count_comparator_kernel(Q, photon_rates, duration, threshold, *,
                            quadrature_size=None, tail_tolerance=1e-15,
                            absolute_tolerance=1e-10, return_info=False):
    """Return K0,K1, each shaped broadcast_batch_shape+(3,3).

    Q is a continuous-time row generator (...,3,3); photon_rates has shape
    (...,3), in counts/second. Duration is a positive scalar in seconds and
    threshold is a nonnegative integer. K0 includes counts <= threshold.
    State 2 may be absorbing; general real three-state generators also work.

    quadrature_size can request a larger power-of-two N for refinement. It
    cannot weaken the Poisson-tail or threshold coverage requirements.
    return_info=True appends a dictionary of numerical diagnostics.
    """
    if isinstance(threshold, (bool, np.bool_)) or not isinstance(threshold, numbers.Integral) or threshold < 0:
        raise ValueError("threshold must be a nonnegative integer, excluding bool")
    threshold = int(threshold)
    if isinstance(duration, (bool, np.bool_)) or not np.isscalar(duration):
        raise ValueError("duration must be a positive finite scalar")
    duration = float(duration)
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("duration must be a positive finite scalar")
    if not 0 < tail_tolerance < 1 or not np.isfinite(absolute_tolerance) or absolute_tolerance <= 0:
        raise ValueError("tail and absolute tolerances must be positive and finite")
    Q = np.asarray(Q, dtype=float)
    rates = np.asarray(photon_rates, dtype=float)
    if Q.ndim < 2 or Q.shape[-2:] != (3, 3) or rates.ndim < 1 or rates.shape[-1] != 3:
        raise ValueError("expected Q (...,3,3) and photon_rates (...,3)")
    if not np.all(np.isfinite(Q)) or not np.all(np.isfinite(rates)) or np.any(rates < 0):
        raise ValueError("Q must be finite and photon rates finite and nonnegative")
    off_diagonal = Q.copy()
    diagonal_index = np.arange(3)
    off_diagonal[..., diagonal_index, diagonal_index] = 0
    if np.any(off_diagonal < 0) or np.any(np.diagonal(Q, axis1=-2, axis2=-1) > 0):
        raise ValueError("Q must have nonnegative off-diagonal and nonpositive diagonal entries")
    generator_scale = np.maximum(1., np.max(np.abs(Q), axis=(-2, -1)))
    if np.any(np.max(np.abs(Q.sum(axis=-1)), axis=-1) > 1e-12 * generator_scale):
        raise ValueError("Q rows must sum to zero")
    try:
        batch_shape = np.broadcast_shapes(Q.shape[:-2], rates.shape[:-1])
    except ValueError as error:
        raise ValueError("Q and photon_rates batch dimensions must broadcast") from error
    generators = np.broadcast_to(Q, batch_shape + (3, 3)).reshape(-1, 3, 3)
    rates = np.broadcast_to(rates, batch_shape + (3,)).reshape(-1, 3)
    if not len(generators):
        raise ValueError("empty batches are not supported")
    max_count = float(np.max(rates) * duration)
    if not np.isfinite(max_count):
        raise ValueError("expected photon counts must be finite")
    # Conditional on any material path, the integrated intensity is no larger
    # than max_count. Poisson domination therefore bounds P(N > threshold)
    # for every starting state, including the joint end-state entries. Avoid
    # quadrature proportional to a huge, legal comparator threshold when that
    # entire response has mass below 1e-18. With at most 1023 internal nodes,
    # the union bound is below 1.023e-15, far below all artifact allowances.
    threshold_tail_bound = float(poisson.sf(threshold, max_count))
    if quadrature_size is None and threshold_tail_bound <= min(tail_tolerance, 1e-18):
        transition = np.asarray([expm(q * duration) for q in generators])
        np.maximum(transition, 0., out=transition)
        zero = transition.reshape(batch_shape + (3, 3))
        one = np.zeros_like(zero)
        info = dict(quadrature_size=0, radius=None,
                    max_expected_count=max_count,
                    poisson_tail_bound=threshold_tail_bound,
                    zero_alias_bound=0., one_alias_bound=0.,
                    threshold_omission_bound=threshold_tail_bound,
                    raw_closure_residual=0., raw_negative_residual=0.,
                    eigen_expm_fallbacks=0,
                    method='bounded-threshold-tail')
        return (zero, one, info) if return_info else (zero, one)
    cutoff = int(poisson.isf(tail_tolerance, max_count)) if max_count else 0
    minimum = max(cutoff + 1, threshold + 2)
    size = 1 << int(minimum).bit_length()
    while poisson.sf(size - 1, max_count) > tail_tolerance:
        size *= 2
    if quadrature_size is not None:
        if (isinstance(quadrature_size, (bool, np.bool_))
                or not isinstance(quadrature_size, numbers.Integral)
                or quadrature_size < size or quadrature_size & (quadrature_size - 1)):
            raise ValueError(f"quadrature_size must be a power of two at least {size}")
        size = int(quadrature_size)
    radius = float(np.exp(-8. / size))
    tail_bound = float(poisson.sf(size - 1, max_count))
    hidden_transition = np.asarray([expm(q * duration) for q in generators])
    zero = np.zeros_like(hidden_transition)
    one = np.zeros_like(hidden_transition)
    analytic_zero_q = np.all(generators == 0, axis=(-2, -1))
    analytic_zero_rates = np.all(rates == 0, axis=-1)
    for index in np.flatnonzero(analytic_zero_q):
        zero[index] = np.diag(poisson.cdf(threshold, rates[index] * duration))
        one[index] = np.diag(poisson.sf(threshold, rates[index] * duration))
    zero[analytic_zero_rates] = hidden_transition[analytic_zero_rates]
    one[analytic_zero_rates] = 0
    work_indices = np.flatnonzero(~(analytic_zero_q | analytic_zero_rates))

    # Conjugate nodes need not be evaluated: pair k with N-1-k.
    angles = 2 * np.pi * (np.arange(size // 2) + .5) / size
    log_inverse_z = 8. / size - 1j * angles
    z = np.exp(-log_inverse_z)
    denominator = np.expm1(log_inverse_z)
    zero_weights = np.expm1((threshold + 1) * log_inverse_z) / denominator
    # At half-step nodes (1/z)**N is exactly -exp(8). Using that identity
    # avoids roundoff from evaluating trigonometric functions at large N*angle.
    one_weights = (-np.exp(8.) - np.exp((threshold + 1) * log_inverse_z)) / denominator
    fallback_count = 0
    # Bound temporary eigendecomposition storage even for large state batches.
    state_chunk = max(1, 16384 // (size // 2))
    for start in range(0, len(work_indices), state_chunk):
        indices = work_indices[start:start + state_chunk]
        tilted = np.broadcast_to(generators[indices, None], (len(indices), size // 2, 3, 3)).astype(complex).copy()
        for state in range(3):
            tilted[..., state, state] += (z[None, :] - 1) * rates[indices, state, None]
        pgf, fallbacks = _matrix_exponential_batched(tilted * duration)
        fallback_count += fallbacks
        zero[indices] = (2. / size) * np.einsum("bkij,k->bij", pgf, zero_weights).real
        one[indices] = (2. / size) * np.einsum("bkij,k->bij", pgf, one_weights).real

    raw_closure = float(np.max(np.abs(zero + one - hidden_transition)))
    raw_negative = float(max(0., -min(zero.min(), one.min(), hidden_transition.min())))
    if raw_closure > absolute_tolerance or raw_negative > absolute_tolerance:
        raise ArithmeticError(
            f"contour residual exceeds tolerance: closure={raw_closure:.3g}, negative={raw_negative:.3g}; "
            "retry with a larger quadrature_size and compare the kernels"
        )
    # Both tails were calculated directly. Reconcile only the recorded small
    # numerical residual, preserving their ratio and enforcing nonnegativity.
    np.maximum(zero, 0., out=zero)
    np.maximum(one, 0., out=one)
    np.maximum(hidden_transition, 0., out=hidden_transition)
    total = zero + one
    missing = total == 0
    if np.any(missing & (hidden_transition > absolute_tolerance)):
        raise ArithmeticError("nonzero hidden transition has zero contour mass")
    scale = np.ones_like(total)
    np.divide(hidden_transition, total, out=scale, where=~missing)
    zero *= scale
    one *= scale
    zero[missing] = hidden_transition[missing]
    residual = hidden_transition - (zero + one)
    larger_zero = zero >= one
    zero[larger_zero] += residual[larger_zero]
    one[~larger_zero] += residual[~larger_zero]
    info = dict(
        quadrature_size=size, radius=radius, max_expected_count=max_count,
        poisson_tail_bound=tail_bound,
        zero_alias_bound=float(np.exp(-8) * tail_bound),
        one_alias_bound=float((1 + np.exp(-8)) * tail_bound),
        raw_closure_residual=raw_closure, raw_negative_residual=raw_negative,
        eigen_expm_fallbacks=fallback_count,
    )
    answer = (zero.reshape(batch_shape + (3, 3)), one.reshape(batch_shape + (3, 3)))
    return (*answer, info) if return_info else answer


contour_comparator_kernel = count_comparator_kernel
comparator_kernel = count_comparator_kernel
