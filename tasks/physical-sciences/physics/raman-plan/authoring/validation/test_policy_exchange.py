"""Independent small-tree checks for full three-state policy exchanges."""

import copy
import io
import time
import unittest
from contextlib import redirect_stdout

import numpy as np
from numpy.testing import assert_allclose

from policy_exchange import TRUE_DOMAIN, _candidate_metrics, _local_candidate_indices, improve


def evaluate(node, actions, K0, K1):
    """Direct recursive continuation values, independent of exchange records."""
    if "domain" in node:
        error = np.broadcast_to((TRUE_DOMAIN != node["domain"])[:, None], (72, 3)).astype(float)
        return error, np.broadcast_to([0., 0., 1.], (72, 3)).copy()
    a = next(i for i, row in enumerate(actions.tolist()) if row == node["action"])
    zero = evaluate(node["children"][0], actions, K0, K1)
    one = evaluate(node["children"][1], actions, K0, K1)
    result = []
    for key in range(2):
        value = np.empty((72, 3))
        for s in range(72):
            value[s] = K0[a, s] @ zero[key][s] + K1[a, s] @ one[key][s]
        result.append(value)
    return tuple(result)


def config():
    return dict(exposure_unit_s=1., readout_s_per_spectrum=0., target_switch_s=0.,
                motion_cost_s=dict(azimuth_per_deg=1., tilt_per_deg_each_axis=0.,
                                   half_wave_plate_per_deg_mod90=0., analyzer_toggle=0., laser_switch=0.))


def limits(budget=20):
    return dict(worst_error_max=.99, mean_error_max=.8,
                worst_memory_loss_max=.4, mean_memory_loss_max=.2,
                time_budget_s=budget, max_policy_depth=3)


def objective(risks, losses, lim):
    return float(
        np.maximum(risks / lim["worst_error_max"] - .95, 0).dot(np.maximum(risks / lim["worst_error_max"] - .95, 0))
        + np.maximum(losses / lim["worst_memory_loss_max"] - .95, 0).dot(np.maximum(losses / lim["worst_memory_loss_max"] - .95, 0))
        + 72 * max(risks.mean() / lim["mean_error_max"] - .95, 0) ** 2
        + 72 * max(losses.mean() / lim["mean_memory_loss_max"] - .95, 0) ** 2
        + 1e-7 * (risks.mean() / lim["mean_error_max"] + losses.mean() / lim["mean_memory_loss_max"])
    )


class PolicyExchangeTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(7441)
        self.actions = np.array([[10, 0, 1, 0], [20, 0, 1, 0], [30, 1, 2, 3], [10, 1, 3, 4]])
        self.catalog = np.array([[1, 10, 532, 0, 0, 0, 0, 0],
                                 [1, 20, 532, 0, 0, 0, 0, 0],
                                 [1, 30, 532, 0, 0, 0, 0, 0]], float)
        transition = rng.random((4, 72, 3, 3))
        transition /= transition.sum(axis=-1, keepdims=True)
        transition[:, :, 2, :] = [0, 0, 1]
        split = rng.uniform(.03, .97, transition.shape)
        self.K0 = transition * split
        self.K1 = transition - self.K0
        self.tree = dict(action=self.actions[0].tolist(), children=[
            dict(action=self.actions[1].tolist(), children=[dict(domain=1), dict(domain=7)]),
            dict(domain=3),
        ])

    def test_all_action_deltas_match_whole_tree_recomputation(self):
        original_error, original_loss = evaluate(self.tree, self.actions, self.K0, self.K1)
        initial = np.zeros((72, 3))
        initial[:, 0] = 1.
        for path, alpha in (([], initial), ([0], np.einsum("sf,sfg->sg", initial, self.K0[0]))):
            current = self.tree if not path else self.tree["children"][0]
            current_error, current_loss = evaluate(current, self.actions, self.K0, self.K1)
            left = evaluate(current["children"][0], self.actions, self.K0, self.K1)
            right = evaluate(current["children"][1], self.actions, self.K0, self.K1)
            risk, loss = _candidate_metrics(alpha, current_error, current_loss,
                                            left[0], left[1], right[0], right[1],
                                            self.K0, self.K1, original_error[:, 0], original_loss[:, 0])
            for candidate in range(len(self.actions)):
                changed = copy.deepcopy(self.tree)
                node = changed if not path else changed["children"][0]
                node["action"] = self.actions[candidate].tolist()
                expected_error, expected_loss = evaluate(changed, self.actions, self.K0, self.K1)
                assert_allclose(risk[candidate], expected_error[:, 0], atol=7e-16, rtol=0)
                assert_allclose(loss[candidate], expected_loss[:, 0], atol=7e-16, rtol=0)

    def test_improvement_decreases_exact_objective_and_relabels(self):
        tree = copy.deepcopy(self.tree)
        before = evaluate(tree, self.actions, self.K0, self.K1)
        with redirect_stdout(io.StringIO()):
            result = improve(tree, self.actions, self.K0, self.K1, self.catalog, config(), limits(), sweeps=3)
        after = evaluate(result, self.actions, self.K0, self.K1)
        self.assertIs(result, tree)
        self.assertLess(objective(after[0][:, 0], after[1][:, 0], limits()),
                        objective(before[0][:, 0], before[1][:, 0], limits()))
        self.assertNotEqual(
            [result["children"][0]["children"][0]["domain"],
             result["children"][0]["children"][1]["domain"], result["children"][1]["domain"]],
            [1, 7, 3],
        )

    def test_material_state_one_is_not_final_loss(self):
        actions = self.actions[:2]
        zero = np.zeros((2, 72, 3, 3))
        zero[0] = [[0, 0, 1], [0, 1, 0], [0, 0, 1]]
        zero[1] = [[0, 1, 0], [0, 1, 0], [0, 0, 1]]
        one = np.zeros_like(zero)
        tree = dict(action=actions[0].tolist(), children=[dict(domain=0), dict(domain=0)])
        with redirect_stdout(io.StringIO()):
            improve(tree, actions, zero, one, self.catalog, config(), limits(), sweeps=1)
        self.assertEqual(tree["action"], actions[1].tolist())
        _, final_loss = evaluate(tree, actions, zero, one)
        assert_allclose(final_loss[:, 0], 0, atol=0)

    def test_rejects_zero_probability_branch_that_exceeds_clock(self):
        tree = copy.deepcopy(self.tree)
        catalog = self.catalog.copy()
        catalog[1, 3] = 180
        zero = self.K0 + self.K1
        one = np.zeros_like(zero)  # The long second acquisition is reached with certainty;
        tree["children"].reverse()  # reverse children so it is instead on the zero-mass bit.
        with self.assertRaisesRegex(ValueError, "clock"):
            improve(tree, self.actions, zero, one, catalog, config(), limits(5), sweeps=0)

    def test_action_exchange_respects_child_subtree_motion(self):
        # A remote action is a better observation but has an impossible return
        # leg to the preserved child subtree. Both root and child swaps must
        # reject it. IDs intentionally differ from catalog row indices.
        catalog = self.catalog.copy()
        catalog[2, 3] = 180
        actions = self.actions[:3].copy()
        actions[2, 2] = 1
        zero = np.zeros((3, 72, 3, 3))
        one = np.zeros_like(zero)
        probabilities = np.array([
            np.where(TRUE_DOMAIN == 0, .8, 0.),
            np.where(TRUE_DOMAIN == 1, .8, 0.),
            np.where(TRUE_DOMAIN == 0, 1., 0.),
        ])
        for a in range(3):
            for s in range(72):
                zero[a, s] = np.eye(3) * probabilities[a, s]
                one[a, s] = np.eye(3) * (1 - probabilities[a, s])
        tree = dict(action=actions[0].tolist(), children=[dict(domain=0),
                    dict(action=actions[1].tolist(), children=[dict(domain=1), dict(domain=2)])])
        with redirect_stdout(io.StringIO()):
            result = improve(tree, actions, zero, one, catalog, config(), limits(5), sweeps=2)
        self.assertNotEqual(result["action"][0], 30)
        self.assertNotEqual(result["children"][1]["action"][0], 30)

    def test_local_filter_keeps_setting_channel_and_integer_neighborhood(self):
        actions = np.array([[10, 0, 10, 50], [10, 0, 13, 71],
                            [10, 0, 14, 70], [10, 0, 13, 72],
                            [10, 1, 10, 50], [20, 0, 10, 50]])
        self.assertEqual(_local_candidate_indices(actions, actions[0]).tolist(), [0, 1])

    def test_local_exchange_rejects_better_nonlocal_comparison(self):
        actions = np.array([[10, 0, 10, 50], [10, 0, 12, 60],
                            [10, 0, 14, 70], [10, 1, 10, 50], [20, 0, 10, 50]])
        zero = np.zeros((len(actions), 72, 3, 3))
        for a, probability in enumerate((.9, .2, 0., 0., 0.)):
            zero[a] = [[1-probability, 0., probability], [0., 1., 0.], [0., 0., 1.]]
        tree = dict(action=actions[0].tolist(), children=[dict(domain=0), dict(domain=0)])
        with redirect_stdout(io.StringIO()):
            improve(tree, actions, zero, np.zeros_like(zero), self.catalog,
                    config(), limits(100), sweeps=1, local=True)
        self.assertEqual(tree["action"], actions[1].tolist())


if __name__ == "__main__":
    started = time.perf_counter()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PolicyExchangeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print(f"Elapsed test time: {time.perf_counter() - started:.4f} s")
    raise SystemExit(0 if result.wasSuccessful() else 1)
