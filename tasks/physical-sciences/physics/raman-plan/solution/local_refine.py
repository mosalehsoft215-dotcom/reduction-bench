"""Local integer dwell/threshold polishing from a caller-supplied policy.

No policy files or seeds are read. Each (setting,channel,dwell) group invokes
the supplied physical Fourier operator once for its union of threshold values.
"""
from __future__ import annotations

from collections import defaultdict
import numbers

import numpy as np

from policy_exchange import improve
from count_kernels import fourier_kernels


OFFSETS = (-6, -3, -1, 0, 1, 3, 6)


def neighbor_actions(tree, threshold_max=1000000):
    """Sorted unique union of integer neighbors of every current action."""
    bank = set()
    seen = set()

    def visit(node):
        if id(node) in seen:
            raise ValueError("policy must be a tree without cycles or shared nodes")
        seen.add(id(node))
        if "domain" in node:
            return
        action = node.get("action")
        if not isinstance(action, (list, tuple, np.ndarray)) or len(action) != 4:
            raise ValueError("internal nodes need a four-integer action")
        if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, numbers.Integral) for x in action):
            raise ValueError("actions must contain integers excluding Boolean values")
        setting, channel, old_dwell, old_threshold = map(int, action)
        if not 1 <= old_dwell <= 32 or not 0 <= old_threshold <= threshold_max:
            raise ValueError("initial action exceeds dwell or threshold limits")
        bank.add((setting, channel, old_dwell, old_threshold))
        for dwell in range(max(1, old_dwell - 3), min(32, old_dwell + 3) + 1):
            center = round(old_threshold * dwell / old_dwell)
            for delta in OFFSETS:
                threshold = min(threshold_max, max(0, center + delta))
                bank.add((setting, channel, dwell, threshold))
        if len(node.get("children", ())) != 2:
            raise ValueError("internal nodes need exactly two children")
        for child in node["children"]:
            visit(child)

    visit(tree)
    return np.asarray(sorted(bank), dtype=np.int64).reshape(-1, 4)


def build_local_bank(tree, catalog, photon_rates, generators, cfg, limits):
    """Return (actions,K0,K1), retaining every original tree action."""
    catalog = np.asarray(catalog)
    rates, q = np.asarray(photon_rates, float), np.asarray(generators, float)
    if catalog.ndim != 2 or catalog.shape[1] != 8:
        raise ValueError("catalog must have shape (K,8)")
    if rates.shape != (len(catalog), q.shape[1], 3, 4) or q.shape != (len(catalog), rates.shape[1], 3, 3):
        raise ValueError("expected rates (K,S,3,4), generators (K,S,3,3)")
    by_id = {int(row[1]): i for i, row in enumerate(catalog)}
    if len(by_id) != len(catalog):
        raise ValueError("catalog setting IDs must be unique")
    threshold_max = int(limits.get("threshold_max", cfg.get("threshold_max", 1000000)))
    actions = neighbor_actions(tree, threshold_max)
    zero = np.empty((len(actions), rates.shape[1], 3, 3))
    one = np.empty_like(zero)
    groups = defaultdict(list)
    for index, action in enumerate(actions):
        setting, channel, dwell, threshold = map(int, action)
        if setting not in by_id:
            raise ValueError("action refers to an unknown catalog setting")
        k = by_id[setting]
        channels = 3 if int(catalog[k, 0]) == 0 else 4
        if not 0 <= channel < channels:
            raise ValueError("action refers to an unavailable channel")
        groups[(setting, channel, dwell)].append(index)
    for (setting, channel, dwell), indices in groups.items():
        k = by_id[setting]
        thresholds = actions[indices, 3]
        lo, hi = fourier_kernels(q[k], rates[k, :, :, channel],
                                 dwell * float(cfg["exposure_unit_s"]), thresholds)
        zero[indices], one[indices] = lo, hi
    print("Local bank", len(actions), "actions", len(groups), "Fourier groups", flush=True)
    return actions, zero, one


def refine_with_bank(tree, catalog, photon_rates, generators, cfg, limits, sweeps=4):
    """Polish in place; return (tree, (actions,K0,K1)) for a final audit."""
    bank = build_local_bank(tree, catalog, photon_rates, generators, cfg, limits)
    actions, zero, one = bank
    if len(actions):
        tree = improve(tree, actions, zero, one, catalog, cfg, limits,
                       sweeps=sweeps, local=True)
    return tree, bank


def refine(tree, catalog, photon_rates, generators, cfg, limits, sweeps=4):
    """Polish in place and return the policy; see refine_with_bank for audit."""
    return refine_with_bank(tree, catalog, photon_rates, generators, cfg, limits, sweeps)[0]
