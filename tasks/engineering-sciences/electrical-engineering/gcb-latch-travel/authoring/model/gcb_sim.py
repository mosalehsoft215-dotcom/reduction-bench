"""Generator-breaker closings recorded by the unit fault recorder: the truth model.

Closing circuit: battery -> control cable (R_cable) -> [interposing relay on some units] -> closing coil (R_coil at 20 C,
copper alpha, inductance L_coil). The coil current rises as (V_coil/R)(1 - exp(-t/tau)), tau = L/R; the latch releases
when it reaches the pick-up current I_pu (t_latch = -tau ln(1 - I_pu R / V_coil)); the armature then moves (a dip in
the current) and the charged spring drives each pole to contact touch after its own travel time. Pole j touches at
coil energisation + t_latch + t_arm + t_travel_j.
Network: generator EMF behind L_g, R_g per phase, neutral earthed; GSU LV delta behind L_t, R_t, floating. Before
contact the two sides are independent; one closed pole shifts the floating side; two closed poles carry a loop
current; three closed poles a three-phase current with no zero sequence. VT voltages are node voltages.
"""
import math

import numpy as np
from scipy.signal import lfilter

FS_TRUE = 4000.0
SUB = 50  # integration steps per sample


def coil_timing(u, V_coil, temp_c):
    R = u["R_coil"] * (1 + 0.00393 * (temp_c - 20.0))
    tau = u["L_coil"] / R
    ratio = u["I_pu"] * R / V_coil
    t_latch = -tau * math.log(1 - ratio)
    return R, tau, t_latch


def pole_times(u, V_coil, temp_c):
    _, _, t_latch = coil_timing(u, V_coil, temp_c)
    return [t_latch + u["t_arm"] + tr for tr in u["t_travel"]]


def simulate(u, close, rng, fs_true=FS_TRUE):
    """One closing: returns time stamps (as the recorder labels them) and channels."""
    fs_rec = u.get("fs_label", fs_true)
    T0, T1 = -0.100, 0.160
    dt = 1.0 / (fs_true * SUB)
    n_fine = int(round((T1 - T0) / dt))
    tf = T0 + dt * np.arange(n_fine)
    # coil energisation relative to the recorded close command (binary input)
    t_e = u.get("relay_delay", 0.0)
    # in the field the whole loop (coil, cable, any extra resistance) sets both the final current and tau
    R_hot = u["R_coil"] * (1 + 0.00393 * (close["temp"] - 20.0))
    R = R_hot + u["R_cable"] + u.get("R_extra", 0.0)
    V_coil = close["V_batt"] * R_hot / R
    tau = u["L_coil"] / R
    t_latch = -tau * math.log(1 - u["I_pu"] * R / close["V_batt"])
    touches = sorted(((t_e + t_latch + u["t_arm"] + tr, k) for k, tr in enumerate(u["t_travel"])))
    # EMFs: generator at bus frequency plus slip, at the synchroniser's angle when the command is given
    phase_g = close["angle_at_cmd"] * math.pi / 180.0
    eg = np.array([math.sqrt(2) * u["V_ph"] * close["vg_pu"] * np.sin(2 * math.pi * (close["f_bus"] + close["slip"]) * tf + phase_g + close["phi0"] - 2 * math.pi * k / 3) for k in range(3)])
    eb = np.array([math.sqrt(2) * u["V_ph"] * np.sin(2 * math.pi * close["f_bus"] * tf + close["phi0"] - 2 * math.pi * k / 3) for k in range(3)])
    L = u["L_g"] + u["L_t"]
    Rr = u["R_g"] + u["R_t"]
    i = np.zeros((3, n_fine))
    vg = eg.copy()
    vb = eb.copy()
    closed = []
    cur = np.zeros(3)
    for n in range(n_fine):
        t = tf[n]
        while len(closed) < 3 and t >= touches[len(closed)][0]:
            closed.append(touches[len(closed)][1])
        d = eg[:, n] - eb[:, n]
        if len(closed) == 2:
            a, b = closed
            loop = d[a] - d[b]
            di = (loop - 2 * Rr * cur[a]) / (2 * L)
            cur[a] += di * dt
            cur[b] = -cur[a]
            dcur = np.zeros(3)
            dcur[a], dcur[b] = di, -di
        elif len(closed) == 3:
            drive = d - d.mean()
            dcur = (drive - Rr * cur) / L
            cur = cur + dcur * dt
        else:
            dcur = np.zeros(3)
        i[:, n] = cur
        vg[:, n] = eg[:, n] - u["L_g"] * dcur - u["R_g"] * cur
        if closed:
            # floating side follows the closed pole(s)
            a = closed[0]
            vb_node = eb[:, n] + u["L_t"] * dcur + u["R_t"] * cur
            shift = vg[a, n] - vb_node[a]
            vb[:, n] = vb_node + shift
        # else vb = eb
    # coil current: exponential rise, armature dip after latch release, cut-off by the auxiliary switch
    tc = tf - t_e
    plateau = close["V_batt"] / R
    ic = np.where(tc > 0, plateau * (1 - np.exp(-np.clip(tc, 0, None) / tau)), 0.0)
    dip = np.where(tc > t_latch, 0.18 * plateau * np.exp(-((tc - t_latch - 0.004) / 0.003) ** 2), 0.0)
    ic = ic - dip
    cut = touches[-1][0] + u["aux_cut"] - t_e
    ic = np.where(tc > cut, ic * np.exp(-np.clip(tc - cut, 0, None) / 0.0015), ic)
    if u.get("clamp_lag"):
        # the clamp's output follows the current through a first-order lag (its low-bandwidth filter setting)
        ic = lfilter([1 - math.exp(-dt / u["clamp_lag"])], [1, -math.exp(-dt / u["clamp_lag"])], ic)
    # sampling: the recorder samples at fs_true but labels samples at fs_rec
    idx = np.arange(0, n_fine, SUB)
    n_rec = idx.size
    n_trig = int(np.flatnonzero(tf[idx] >= 0.0)[0])  # the export counts time from the trigger sample at its own rate
    t_lab = (np.arange(n_rec) - n_trig) / fs_rec
    out = {"t_ms": 1000.0 * t_lab}
    names = "ABC"
    for k in range(3):
        out["VG_" + names[k]] = vg[k, idx] / 1000.0 + rng.normal(0, 0.008, n_rec)
        out["VB_" + names[k]] = vb[k, idx] / 1000.0 + rng.normal(0, 0.008, n_rec)
        out["IG_" + names[k]] = i[k, idx] / 1000.0 + rng.normal(0, 0.003, n_rec)
    out["ICOIL"] = ic[idx] + rng.normal(0, 0.02, n_rec)
    out["BI_CLOSE"] = (tf[idx] >= 0.0).astype(int)
    out["BI_52A"] = (tf[idx] >= touches[-1][0] + u["aux_a"]).astype(int)
    truth = {"t_e": t_e, "V_coil": V_coil, "touch": {names[k]: tk for tk, k in touches}, "t_latch": t_latch}
    return out, truth
