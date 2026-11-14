"""Small analytic and structural checks; run this file directly."""

import time
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm
from scipy.stats import poisson

from independent_count_kernel import count_comparator_kernel


class IndependentCountKernelTests(unittest.TestCase):
    def test_zero_generator_poisson_diagonal(self):
        rates = np.array([0.0, 0.7, 13.0])
        for threshold in (0, 2, 12, 35):
            zero, one = count_comparator_kernel(np.zeros((3, 3)), rates, 1.3, threshold)
            assert_allclose(zero, np.diag(poisson.cdf(threshold, rates * 1.3)), atol=2e-14, rtol=2e-13)
            assert_allclose(one, np.diag(poisson.sf(threshold, rates * 1.3)), atol=2e-14, rtol=2e-13)

    def test_equal_rates_factorize_and_batches_broadcast(self):
        generators = np.array([
            [[-0.7, 0.5, 0.2], [0.3, -0.9, 0.6], [0, 0, 0]],
            [[-0.1, 0.1, 0], [0.2, -0.4, 0.2], [0, 0, 0]],
        ])[:, None, :, :]
        rate_levels = np.array([0.2, 2.0, 17.0])
        rates = np.repeat(rate_levels[:, None], 3, axis=-1)
        duration = 2.1
        zero, one = count_comparator_kernel(generators, rates, duration, 4)
        self.assertEqual(zero.shape, (2, 3, 3, 3))
        for i in range(2):
            transition = expm(generators[i, 0] * duration)
            for j, rate in enumerate(rate_levels):
                assert_allclose(zero[i, j], poisson.cdf(4, rate * duration) * transition, atol=3e-14, rtol=3e-13)
                assert_allclose(one[i, j], poisson.sf(4, rate * duration) * transition, atol=3e-14, rtol=3e-13)

    def test_direct_small_tail_is_not_a_complement(self):
        zero, one = count_comparator_kernel(np.zeros((3, 3)), np.full(3, 0.7), 1.0, 15)
        expected = poisson.sf(15, 0.7)
        self.assertLess(expected, 1e-16)
        self.assertGreater(one[0, 0], 0)
        # Sparse Taylor propagation guarantees absolute precision. At this
        # sub-epsilon tail its relative error can be larger; it still retains
        # the tail instead of rounding it to zero by subtraction.
        assert_allclose(np.diag(one), expected, atol=1e-19, rtol=0)
        # A subtraction implementation would give zero for this tail.
        self.assertEqual(1 - zero[0, 0], 0)

    def test_general_conservation_and_absorbing_row(self):
        generator = np.array([[-0.3, 0.2, 0.1], [0.4, -0.6, 0.2], [0, 0, 0]])
        rates = np.array([0.02, 9.0, 0.3])
        zero, one = count_comparator_kernel(generator, rates, 2.7, 8)
        self.assertGreaterEqual(min(zero.min(), one.min()), -1e-15)
        assert_allclose(zero + one, expm(generator * 2.7), atol=2e-15, rtol=2e-14)
        assert_allclose((zero + one).sum(axis=1), np.ones(3), atol=2e-15)
        assert_allclose(zero[2], [0, 0, poisson.cdf(8, 0.3 * 2.7)], atol=2e-14)
        assert_allclose(one[2], [0, 0, poisson.sf(8, 0.3 * 2.7)], atol=2e-14)

    def test_multistep_order_matters_and_memory_persists(self):
        generator = np.array([[-0.8, 0.7, 0.1], [0.2, -0.5, 0.3], [0, 0, 0]])
        a_zero, a_one = count_comparator_kernel(generator, [12, 0.3, 1.2], 0.8, 2)
        b_zero, b_one = count_comparator_kernel(generator, [0.4, 10, 1.2], 1.1, 3)
        initial = np.array([1.0, 0.0, 0.0])
        sequence_ab = initial @ a_one @ b_zero
        sequence_ba = initial @ b_zero @ a_one
        self.assertGreater(np.max(np.abs(sequence_ab - sequence_ba)), 0.01)
        # Resetting or independently marginalizing material state gives a wrong
        # two-bit likelihood. Sum all four outcomes to recover hidden evolution.
        reset_likelihood = (initial @ a_one).sum() * (initial @ b_zero).sum()
        self.assertGreater(abs(sequence_ab.sum() - reset_likelihood), 0.01)
        total = sum(initial @ a @ b for a in (a_zero, a_one) for b in (b_zero, b_one))
        assert_allclose(total, initial @ expm(generator * 1.9), atol=2e-14)

    def test_input_validation(self):
        Q = np.zeros((3, 3))
        for duration, threshold in ((0, 1), (-1, 1), (np.inf, 1), (1, -1), (1, True), (1, 1.2)):
            with self.assertRaises(ValueError):
                count_comparator_kernel(Q, [1, 2, 3], duration, threshold)
        with self.assertRaises(ValueError):
            count_comparator_kernel(np.ones((3, 3)), [1, 2, 3], 1, 1)
        with self.assertRaises(ValueError):
            count_comparator_kernel(Q, [-1, 2, 3], 1, 1)


if __name__ == "__main__":
    start = time.perf_counter()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(IndependentCountKernelTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print(f"Elapsed test time: {time.perf_counter() - start:.4f} s")
    raise SystemExit(0 if result.wasSuccessful() else 1)
