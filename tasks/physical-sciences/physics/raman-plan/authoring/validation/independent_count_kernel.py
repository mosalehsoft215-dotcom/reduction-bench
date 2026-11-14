"""Independent count-index comparator kernels for a three-state material CTMC.

All generators use row-vector conventions.  Photons arrive at the rate of the
current material state and do not themselves change that state.  A count bin
``threshold + 1`` collects all overflow events while the material CTMC keeps
evolving.  Thus both comparator tails are computed directly, without forming
one as a difference of nearly equal probabilities.

This is a validation operator, not a policy search implementation.  Its sparse
state space has ``3 * (threshold + 2)`` entries.  Like expm_multiply itself, its
accuracy is absolute rather than a guarantee of relative precision in tails
far below machine epsilon.
"""

from __future__ import annotations

import numbers

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply


def count_comparator_kernel(Q, photon_rates, duration, threshold):
    """Return ``(K0, K1)`` with shape ``broadcast_batch_shape + (3, 3)``.

    ``Q`` has shape ``(..., 3, 3)`` and contains a continuous-time row
    generator in inverse seconds. ``photon_rates`` has shape ``(..., 3)`` and
    contains nonnegative photon rates per second. Batch axes broadcast.
    ``duration`` is a positive finite scalar in seconds; ``threshold`` is a
    nonnegative integer, excluding Boolean values.

    ``Ky[i,j] = P(bit=y, final_state=j | initial_state=i)`` where bit zero
    means ``N <= threshold``. State 2 may be absorbing; the operator also
    supports other three-state generators. To roundoff, ``K0 + K1`` equals
    ``expm(Q * duration)``. A tiny closure residual is assigned to the larger
    comparator probability for each material transition, preserving the
    directly calculated small tail.
    """
    if isinstance(threshold, (bool, np.bool_)) or not isinstance(threshold, numbers.Integral):
        raise ValueError("threshold must be a nonnegative integer, excluding bool")
    threshold = int(threshold)
    if threshold < 0:
        raise ValueError("threshold must be nonnegative")
    if isinstance(duration, (bool, np.bool_)) or not np.isscalar(duration):
        raise ValueError("duration must be a positive finite scalar")
    duration = float(duration)
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("duration must be a positive finite scalar")

    Q = np.asarray(Q, dtype=float)
    photon_rates = np.asarray(photon_rates, dtype=float)
    if Q.ndim < 2 or Q.shape[-2:] != (3, 3):
        raise ValueError("Q must have final dimensions (3, 3)")
    if photon_rates.ndim < 1 or photon_rates.shape[-1] != 3:
        raise ValueError("photon_rates must have final dimension 3")
    if not np.all(np.isfinite(Q)) or not np.all(np.isfinite(photon_rates)):
        raise ValueError("Q and photon_rates must be finite")
    if np.any(photon_rates < 0):
        raise ValueError("photon_rates must be nonnegative")
    off_diagonal = Q.copy()
    diagonal_index = np.arange(3)
    off_diagonal[..., diagonal_index, diagonal_index] = 0
    if np.any(off_diagonal < 0) or np.any(np.diagonal(Q, axis1=-2, axis2=-1) > 0):
        raise ValueError("Q must have nonnegative off-diagonal and nonpositive diagonal entries")
    generator_scale = np.maximum(1.0, np.max(np.abs(Q), axis=(-2, -1)))
    if np.any(np.max(np.abs(Q.sum(axis=-1)), axis=-1) > 1e-12 * generator_scale):
        raise ValueError("rows of Q must sum to zero")

    try:
        batch_shape = np.broadcast_shapes(Q.shape[:-2], photon_rates.shape[:-1])
    except ValueError as error:
        raise ValueError("Q and photon_rates batch dimensions must broadcast") from error
    Q_batch = np.broadcast_to(Q, batch_shape + (3, 3)).reshape(-1, 3, 3)
    rate_batch = np.broadcast_to(photon_rates, batch_shape + (3,)).reshape(-1, 3)
    kernel_zero = np.empty((len(Q_batch), 3, 3), dtype=float)
    kernel_one = np.empty_like(kernel_zero)

    bin_count = threshold + 2
    # Row c -> c+1 at unit rate, except in the terminal overflow bin.
    count_diagonal = np.r_[-np.ones(bin_count - 1), 0.0]
    count_generator = sparse.diags(
        (count_diagonal, np.ones(bin_count - 1)), (0, 1),
        shape=(bin_count, bin_count), format="csr",
    )
    count_identity = sparse.eye(bin_count, format="csr")
    starts = np.zeros((3 * bin_count, 3), dtype=float)
    starts[:3, :] = np.eye(3)

    for index, (material_generator, rates) in enumerate(zip(Q_batch, rate_batch)):
        joint_generator = (
            sparse.kron(count_identity, sparse.csr_matrix(material_generator), format="csr")
            + sparse.kron(count_generator, sparse.diags(rates), format="csr")
        )
        # Column-action exponential propagates all three starting states at once.
        trace = bin_count * np.trace(material_generator) - (bin_count - 1) * rates.sum()
        propagated = expm_multiply(
            joint_generator.T * duration, starts, traceA=trace * duration,
        ).reshape(bin_count, 3, 3)
        zero = propagated[:-1].sum(axis=0).T.copy()
        one = propagated[-1].T.copy()
        hidden_transition = expm(material_generator * duration)
        residual = hidden_transition - (zero + one)
        larger_zero = zero >= one
        zero[larger_zero] += residual[larger_zero]
        one[~larger_zero] += residual[~larger_zero]
        kernel_zero[index], kernel_one[index] = zero, one

    return (
        kernel_zero.reshape(batch_shape + (3, 3)),
        kernel_one.reshape(batch_shape + (3, 3)),
    )


# Short alias for callers comparing several independently derived operators.
comparator_kernel = count_comparator_kernel
