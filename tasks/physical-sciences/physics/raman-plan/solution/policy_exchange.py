"""Exact local exchanges for a finite-initial-state, three-material-state policy.

Policies are binary dictionaries with either ``domain`` or ``action`` and
``children``. Actions are [setting_id, channel, dwell_units, threshold]. Catalog
rows are [target, setting_id, laser_nm, azimuth, tilt_x, tilt_y, plate, crossed].
Every initial state starts in material state 0. The reported domain is its
original label (state_index // (state_count // 12)), and final material state 2 defines loss.
No policy, selected setting IDs, or candidate action seed is embedded here.
"""

from __future__ import annotations

import numbers

import numpy as np


TRUE_DOMAIN = np.repeat(np.arange(12), 6)
LIMIT_KEYS = (
    "worst_error_max", "mean_error_max",
    "worst_memory_loss_max", "mean_memory_loss_max",
)


def motion_seconds(a, b, cfg):
    """Serial commanded-motion clock, in seconds, for two catalog rows."""
    costs = cfg["motion_cost_s"]
    azimuth = abs(a[3] - b[3]) % 360
    plate = abs(a[6] - b[6]) % 90
    return float(
        costs["azimuth_per_deg"] * min(azimuth, 360 - azimuth)
        + costs["tilt_per_deg_each_axis"] * (abs(a[4] - b[4]) + abs(a[5] - b[5]))
        + costs["half_wave_plate_per_deg_mod90"] * min(plate, 90 - plate)
        + costs["analyzer_toggle"] * (a[7] != b[7])
        + costs["laser_switch"] * (a[2] != b[2])
        + cfg["target_switch_s"] * (a[0] != b[0])
    )


def _candidate_metrics(alpha, current_error, current_loss,
                       left_error, left_loss, right_error, right_loss,
                       K0, K1, risks, losses):
    """Exact whole-policy deltas, retaining only (actions,S) result arrays."""
    old_risk = np.einsum("sf,sf->s", alpha, current_error)
    old_loss = np.einsum("sf,sf->s", alpha, current_loss)
    # The state-local outer products are tiny and reused across all actions.
    # A direct two-operand contraction avoids moving the large action bank.
    next_risk = np.einsum("asfg,sfg->as", K0, alpha[:,:,None]*left_error[:,None,:], optimize=False)
    next_risk += np.einsum("asfg,sfg->as", K1, alpha[:,:,None]*right_error[:,None,:], optimize=False)
    next_loss = np.einsum("asfg,sfg->as", K0, alpha[:,:,None]*left_loss[:,None,:], optimize=False)
    next_loss += np.einsum("asfg,sfg->as", K1, alpha[:,:,None]*right_loss[:,None,:], optimize=False)
    next_risk += risks[None, :] - old_risk[None, :]
    next_loss += losses[None, :] - old_loss[None, :]
    return next_risk, next_loss


def _local_candidate_indices(actions, current_action):
    """Integer dwell/threshold neighbors, keeping setting and channel fixed."""
    setting, channel, dwell, threshold = current_action
    center = np.rint(threshold * actions[:, 2] / dwell)
    mask = ((actions[:, 0] == setting) & (actions[:, 1] == channel)
            & (abs(actions[:, 2] - dwell) <= 3)
            & (abs(actions[:, 3] - center) <= 6))
    return np.flatnonzero(mask)


def improve(tree, actions, K0, K1, catalog, cfg, limits, sweeps=8, local=False):
    """Improve ``tree`` in place and return it, using exact candidate deltas.

    K0 and K1 have shape (actions,S,3,3), in row-start/column-end convention.
    All candidate actions are considered at each internal node unless local
    is true, when only same-setting/channel integer neighbors are considered.
    Dwell differs by at most 3 and threshold by at most 6 from the current
    threshold scaled to the candidate dwell. Leaf labels
    are also exchanged. The score penalizes every conditional risk/loss and
    both means relative to the four limits, with a small safety margin.
    The initial tree must obey depth and clock limits; every accepted exchange
    preserves these limits for every branch, independently of its probability.
    """
    if isinstance(sweeps, (bool, np.bool_)) or not isinstance(sweeps, numbers.Integral) or sweeps < 0:
        raise ValueError("sweeps must be a nonnegative integer")
    if not isinstance(local, (bool, np.bool_)):
        raise ValueError("local must be Boolean")
    state_count=np.asarray(K0).shape[1]
    TRUE_DOMAIN=np.repeat(np.arange(12),state_count//12)
    actions = np.asarray(actions)
    catalog = np.asarray(catalog, dtype=float)
    K0, K1 = np.asarray(K0), np.asarray(K1)
    if actions.ndim != 2 or actions.shape[1] != 4 or not len(actions):
        raise ValueError("actions must have shape (A,4), with A positive")
    if not np.all(np.isfinite(actions)) or not np.all(actions == np.floor(actions)):
        raise ValueError("actions must contain finite integers")
    actions = actions.astype(np.int64, copy=False)
    if catalog.ndim != 2 or catalog.shape[1] != 8:
        raise ValueError("catalog must have shape (settings,8)")
    if K0.shape != (len(actions), state_count, 3, 3) or K1.shape != K0.shape:
        raise ValueError("kernels must have shape (actions,S,3,3)")
    lookup = {tuple(row): i for i, row in enumerate(actions.tolist())}
    if len(lookup) != len(actions):
        raise ValueError("candidate actions must be unique")
    setting_index = {int(row[1]): i for i, row in enumerate(catalog)}
    if len(setting_index) != len(catalog):
        raise ValueError("catalog setting IDs must be unique")
    try:
        action_settings = np.array([setting_index[int(a[0])] for a in actions], dtype=int)
    except KeyError as error:
        raise ValueError("candidate action references an unknown setting") from error
    targets = np.asarray([limits[k] for k in LIMIT_KEYS], dtype=float)
    if np.any(~np.isfinite(targets)) or np.any(targets <= 0):
        raise ValueError("the four risk/loss limits must be finite and positive")
    budget = float(limits["time_budget_s"] if "time_budget_s" in limits else cfg["time_budget_s"])
    if not np.isfinite(budget) or budget < 0:
        raise ValueError("time budget must be finite and nonnegative")
    depth_max = int(limits.get("max_policy_depth", cfg.get("max_policy_depth", 7)))
    movement = np.asarray([[motion_seconds(a, b, cfg) for b in catalog] for a in catalog])
    durations = actions[:, 2] * float(cfg["exposure_unit_s"]) + float(cfg["readout_s_per_spectrum"])
    if np.any(actions[:, 2] <= 0) or np.any(durations <= 0):
        raise ValueError("candidate dwell and acquisition duration must be positive")

    def score(risk, loss):
        risk = np.maximum(risk, 0)
        loss = np.maximum(loss, 0)
        penalty = np.maximum(risk / targets[0] - .95, 0) ** 2
        penalty += np.maximum(loss / targets[2] - .95, 0) ** 2
        risk_mean, loss_mean = risk.mean(axis=-1), loss.mean(axis=-1)
        return (
            penalty.sum(axis=-1)
            + state_count * np.maximum(risk_mean / targets[1] - .95, 0) ** 2
            + state_count * np.maximum(loss_mean / targets[3] - .95, 0) ** 2
            + 1e-7 * (risk_mean / targets[1] + loss_mean / targets[3])
        )

    def backward(node, info, depth):
        if id(node) in info:
            raise ValueError("policy must be a tree, without shared nodes or cycles")
        info[id(node)] = None  # Mark before descent to detect cycles.
        if depth > depth_max:
            raise ValueError("initial policy exceeds the acquisition-count limit")
        if "domain" in node:
            label = node["domain"]
            if isinstance(label, (bool, np.bool_)) or not isinstance(label, numbers.Integral) or not 0 <= label < 12:
                raise ValueError("terminal domain must be an integer from 0 to 11")
            error = np.repeat((TRUE_DOMAIN != label)[:, None], 3, axis=1).astype(float)
            loss = np.broadcast_to([0., 0., 1.], (state_count, 3))
            record = dict(node=node, error=error, loss=loss, depth=depth, remaining=0., setting=None)
        else:
            if depth == depth_max or "action" not in node or len(node.get("children", [])) != 2:
                raise ValueError("internal policy nodes need an action and two children within the depth limit")
            try:
                action = lookup[tuple(node["action"])]
            except KeyError as error:
                raise ValueError("tree action is absent from the candidate action bank") from error
            left = backward(node["children"][0], info, depth + 1)
            right = backward(node["children"][1], info, depth + 1)
            error = np.einsum("sfg,sg->sf", K0[action], left["error"])
            error += np.einsum("sfg,sg->sf", K1[action], right["error"])
            loss = np.einsum("sfg,sg->sf", K0[action], left["loss"])
            loss += np.einsum("sfg,sg->sf", K1[action], right["loss"])
            setting = int(action_settings[action])
            suffix = [0. if child["setting"] is None else movement[setting, child["setting"]] + child["remaining"]
                      for child in (left, right)]
            record = dict(node=node, error=error, loss=loss, depth=depth,
                          remaining=float(durations[action] + max(suffix)),
                          setting=setting, action=action)
        info[id(node)] = record
        return record

    def forward(node, alpha, previous, before, info):
        record = info[id(node)]
        record.update(alpha=alpha, previous=previous, before=before)
        if "domain" in node:
            return
        action, setting = record["action"], record["setting"]
        end = before + (movement[previous, setting] if previous is not None else 0.) + durations[action]
        for bit, kernel in enumerate((K0[action], K1[action])):
            following = np.einsum("sf,sfg->sg", alpha, kernel)
            forward(node["children"][bit], following, setting, end, info)

    def rebuild():
        info = {}
        root = backward(tree, info, 0)
        alpha = np.zeros((state_count, 3))
        alpha[:, 0] = 1.
        forward(tree, alpha, None, 0., info)
        if root["remaining"] > budget + 1e-9:
            raise ValueError("policy exceeds the hard clock budget")
        return info, root["error"][:, 0], root["loss"][:, 0]

    info, risks, losses = rebuild()
    terminal_zero = K0.sum(axis=-1)
    terminal_one = K1.sum(axis=-1)
    terminal_lost = K0[...,2] + K1[...,2]
    parents = {id(tree): None}
    def record_parents(node):
        for child in node.get('children', []):
            parents[id(child)] = node
            record_parents(child)
    record_parents(tree)

    def refresh(node):
        """Update changed ancestors and affected forward beliefs only."""
        changed = node
        while node is not None:
            record = info[id(node)]
            if 'domain' in node:
                record['error'] = np.repeat((TRUE_DOMAIN != node['domain'])[:, None], 3, axis=1).astype(float)
            else:
                action = lookup[tuple(node['action'])]
                left, right = [info[id(child)] for child in node['children']]
                for name in ('error', 'loss'):
                    value = np.einsum('sfg,sg->sf', K0[action], left[name])
                    value += np.einsum('sfg,sg->sf', K1[action], right[name])
                    record[name] = value
                setting = int(action_settings[action])
                suffix = [0. if child['setting'] is None else movement[setting, child['setting']] + child['remaining'] for child in (left,right)]
                record.update(action=action,setting=setting,remaining=float(durations[action]+max(suffix)))
            node = parents[id(node)]
        record = info[id(changed)]
        forward(changed, record['alpha'], record['previous'], record['before'], info)
        root = info[id(tree)]
        assert root['remaining'] <= budget + 1e-9
        return info, root['error'][:,0], root['loss'][:,0]

    best = float(score(risks, losses))
    print("Policy exchange initial", best, dict(worst_error=float(risks.max()), mean_error=float(risks.mean()), worst_memory_loss=float(losses.max()), mean_memory_loss=float(losses.mean())), flush=True)
    for sweep in range(int(sweeps)):
        changed = 0
        nodes = [record["node"] for record in sorted(info.values(), key=lambda record: -record["depth"])]
        for node_number, node in enumerate(nodes):
            if node_number and node_number % 128 == 0:
                print("Exchange progress", sweep, node_number, len(nodes), float(risks.mean()), float(risks.max()), flush=True)
            current = info[id(node)]
            alpha = current["alpha"]
            if "domain" in node:
                labels = np.arange(12)[:, None]
                old_error = (TRUE_DOMAIN != node["domain"]).astype(float)
                next_risks = risks[None, :] + alpha.sum(axis=1)[None, :] * (
                    (TRUE_DOMAIN[None, :] != labels).astype(float) - old_error[None, :]
                )
                scores = score(next_risks, np.broadcast_to(losses, next_risks.shape))
                choice = int(np.argmin(scores))
                if scores[choice] < best - 1e-10:
                    node["domain"] = choice
                    changed += 1
                    info, risks, losses = refresh(node)
                    best = float(score(risks, losses))
                continue

            left, right = [info[id(child)] for child in node["children"]]
            candidates = (_local_candidate_indices(actions, actions[current["action"]])
                          if local else np.arange(len(actions)))
            if not len(candidates):
                continue
            candidate_zero = K0[candidates] if local else K0
            candidate_one = K1[candidates] if local else K1
            if 'domain' in left['node'] and 'domain' in right['node']:
                z = terminal_zero[candidates] if local else terminal_zero
                u = terminal_one[candidates] if local else terminal_one
                dead = terminal_lost[candidates] if local else terminal_lost
                next_risks = np.einsum('asf,sf->as',z,alpha*left['error'][:,0,None])
                next_risks += np.einsum('asf,sf->as',u,alpha*right['error'][:,0,None])
                next_losses = np.einsum('asf,sf->as',dead,alpha)
                next_risks += risks[None,:]-np.einsum('sf,sf->s',alpha,current['error'])[None,:]
                next_losses += losses[None,:]-np.einsum('sf,sf->s',alpha,current['loss'])[None,:]
            else:
                next_risks, next_losses = _candidate_metrics(
                    alpha, current["error"], current["loss"],
                    left["error"], left["loss"], right["error"], right["loss"],
                    candidate_zero, candidate_one, risks, losses,
                )
            scores = score(next_risks, next_losses)
            candidate_settings = action_settings[candidates]
            end = current["before"] + durations[candidates]
            if current["previous"] is not None:
                end = end + movement[current["previous"], candidate_settings]
            suffixes = [np.zeros(len(candidates)) if child["setting"] is None
                        else movement[candidate_settings, child["setting"]] + child["remaining"]
                        for child in (left, right)]
            total = end + np.maximum(suffixes[0], suffixes[1])
            scores[total > budget + 1e-9] = np.inf
            choice_index = int(np.argmin(scores))
            choice = int(candidates[choice_index])
            if scores[choice_index] < best - 1e-10:
                node["action"] = actions[choice].tolist()
                changed += 1
                info, risks, losses = refresh(node)
                best = float(score(risks, losses))

        # A full recomputation audits the incremental update once per sweep.
        complete, checked_risks, checked_losses = rebuild()
        assert np.allclose(checked_risks, risks, atol=1e-12, rtol=1e-12)
        assert np.allclose(checked_losses, losses, atol=1e-12, rtol=1e-12)
        for key in info:
            assert np.allclose(info[key]['alpha'], complete[key]['alpha'], atol=1e-12, rtol=1e-12)
            assert abs(info[key]['before']-complete[key]['before']) < 1e-10
        info, risks, losses = complete, checked_risks, checked_losses
        best = float(score(risks, losses))
        print("Policy exchange", sweep, changed,
              dict(worst_error=float(risks.max()), mean_error=float(risks.mean()),
                   worst_memory_loss=float(losses.max()), mean_memory_loss=float(losses.mean())),
              flush=True)
        if changed == 0 or (
            risks.max() <= .98 * targets[0] and risks.mean() <= .98 * targets[1]
            and losses.max() <= .98 * targets[2] and losses.mean() <= .98 * targets[3]
        ):
            break
    return tree
