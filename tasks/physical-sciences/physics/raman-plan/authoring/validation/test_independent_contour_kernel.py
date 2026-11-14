"""Analytic, sparse, and sequence checks for independent contour inversion."""

import time
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm
from scipy.stats import poisson

from independent_contour_kernel import count_comparator_kernel as contour
from independent_count_kernel import count_comparator_kernel as sparse_kernel


Q = np.array([[-.8, .7, .1], [.2, -.5, .3], [0., 0., 0.]])


class IndependentContourTests(unittest.TestCase):
    def test_maximum_legal_threshold_uses_bounded_memory(self):
        zero, one, info = contour(Q, [12., .3, 1.2], 22.4,
                                  1000000, return_info=True)
        self.assertEqual(info['method'], 'bounded-threshold-tail')
        self.assertEqual(info['quadrature_size'], 0)
        self.assertLessEqual(info['threshold_omission_bound'], 1e-18)
        assert_allclose(zero, expm(Q * 22.4), atol=2e-15, rtol=0)
        assert_allclose(one, 0., atol=0., rtol=0)

    def test_tail_bound_agrees_with_analytic_equal_rate_law(self):
        zero, one, info = contour(Q, [1000.] * 3, 1., 1500,
                                  return_info=True)
        self.assertEqual(info['method'], 'bounded-threshold-tail')
        bound = poisson.sf(1500, 1000.)
        self.assertGreater(bound, 0.)
        self.assertEqual(info['threshold_omission_bound'], bound)
        transition = expm(Q)
        expected_one = bound * transition
        self.assertLessEqual(np.max(abs(one - expected_one)), bound)
        assert_allclose(zero + one, transition, atol=2e-15, rtol=0)

    def test_tail_fast_path_broadcasts_all_initial_states(self):
        generators = np.stack([Q, .4 * Q])[:, None]
        rates = np.array([[3., 1., 2.], [8., 4., 6.]])
        zero, one, info = contour(generators, rates, .7, 1000000,
                                  return_info=True)
        self.assertEqual(zero.shape, (2, 2, 3, 3))
        self.assertEqual(info['quadrature_size'], 0)
        for i in range(2):
            for j in range(2):
                assert_allclose(zero[i, j], expm(generators[i, 0] * .7),
                                atol=2e-15, rtol=0)
        assert_allclose(one, 0., atol=0., rtol=0)

    def test_zero_q_poisson_diagonal_at_high_flux(self):
        rates = np.array([0., 700., 1000.])
        for threshold in (0, 700, 1000, 1100):
            zero, one = contour(np.zeros((3, 3)), rates, 1., threshold)
            assert_allclose(zero, np.diag(poisson.cdf(threshold, rates)), atol=2e-15, rtol=0)
            assert_allclose(one, np.diag(poisson.sf(threshold, rates)), atol=2e-15, rtol=0)

    def test_equal_rates_analytic_at_mean_count_1000(self):
        zero, one, info = contour(Q, [1000.] * 3, 1., 1000, return_info=True)
        transition = expm(Q)
        assert_allclose(zero, poisson.cdf(1000, 1000.) * transition, atol=1e-10, rtol=0)
        assert_allclose(one, poisson.sf(1000, 1000.) * transition, atol=1e-10, rtol=0)
        self.assertLessEqual(info["poisson_tail_bound"], 1e-15)
        self.assertLess(info["raw_closure_residual"], 1e-10)

    def test_sparse_third_check_at_small_and_high_counts(self):
        for rates, duration, threshold in (([3., 9., 1.], 1., 4), ([220., 1000., 170.], 1., 700)):
            zero, one = contour(Q, rates, duration, threshold)
            sparse_zero, sparse_one = sparse_kernel(Q, rates, duration, threshold)
            assert_allclose(zero, sparse_zero, atol=1e-10, rtol=0)
            assert_allclose(one, sparse_one, atol=1e-10, rtol=0)
            assert_allclose(zero + one, expm(Q * duration), atol=3e-15, rtol=0)
            self.assertGreaterEqual(min(zero.min(), one.min()), 0.)

    def test_defective_equal_rate_generator_uses_expm_fallback(self):
        defective = np.array([[-1., 1., 0.], [0., -1., 1.], [0., 0., 0.]])
        zero, one, info = contour(defective, [9.] * 3, .8, 7, return_info=True)
        self.assertGreater(info["eigen_expm_fallbacks"], 0)
        transition = expm(defective * .8)
        assert_allclose(zero, poisson.cdf(7, 7.2) * transition, atol=1e-10, rtol=0)
        assert_allclose(one, poisson.sf(7, 7.2) * transition, atol=1e-10, rtol=0)

    def test_broadcasting_and_contour_refinement(self):
        generators = np.array([Q, Q * .3])[:, None]
        rates = np.array([[2., 8., .2], [20., 40., 2.], [100., 90., 10.]])
        zero, one, info = contour(generators, rates, .7, 25, return_info=True)
        self.assertEqual(zero.shape, (2, 3, 3, 3))
        refined_zero, refined_one = contour(generators, rates, .7, 25,
                                           quadrature_size=2 * info["quadrature_size"])
        assert_allclose(zero, refined_zero, atol=1e-10, rtol=0)
        assert_allclose(one, refined_one, atol=1e-10, rtol=0)
        for i in range(2):
            for j in range(3):
                expected_zero, expected_one = sparse_kernel(generators[i, 0], rates[j], .7, 25)
                assert_allclose(zero[i, j], expected_zero, atol=1e-10, rtol=0)
                assert_allclose(one[i, j], expected_one, atol=1e-10, rtol=0)

    def test_noncommuting_sequence_matches_sparse(self):
        a = contour(Q, [12., .3, 1.2], .8, 2)
        b = contour(Q, [.4, 10., 1.2], 1.1, 3)
        sparse_a = sparse_kernel(Q, [12., .3, 1.2], .8, 2)
        sparse_b = sparse_kernel(Q, [.4, 10., 1.2], 1.1, 3)
        initial = np.array([1., 0., 0.])
        ab = initial @ a[1] @ b[0]
        ba = initial @ b[0] @ a[1]
        self.assertGreater(np.max(np.abs(ab - ba)), .01)
        assert_allclose(ab, initial @ sparse_a[1] @ sparse_b[0], atol=1e-10, rtol=0)
        combined = sum(initial @ x @ y for x in a for y in b)
        assert_allclose(combined, initial @ expm(Q * 1.9), atol=3e-15, rtol=0)

    def test_invalid_inputs_and_underresolved_quadrature(self):
        for duration, threshold in ((0, 1), (1, -1), (1, True), (1, .5)):
            with self.assertRaises(ValueError):
                contour(Q, [1, 2, 3], duration, threshold)
        with self.assertRaises(ValueError):
            contour(Q, [1, 2, 3], 1, 4, quadrature_size=4)
        with self.assertRaises(ValueError):
            contour(np.ones((3, 3)), [1, 2, 3], 1, 4)


if __name__ == "__main__":
    started = time.perf_counter()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(IndependentContourTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print(f"Elapsed test time: {time.perf_counter() - started:.4f} s")
    raise SystemExit(0 if result.wasSuccessful() else 1)
