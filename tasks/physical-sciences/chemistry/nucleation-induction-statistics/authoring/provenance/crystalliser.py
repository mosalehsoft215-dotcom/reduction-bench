"""Physical model of a parallel-vial crystalliser run (four blocks of eight stirred vials, repeated cycles).

Each cycle heats the blocks at a fixed rate to a dissolution temperature, holds it, cools to the hold temperature and
holds there. In every vial:

- primary nucleation is a Poisson process with hazard V J(S(t)), where S(t) = c / c*(T(t)) follows the vial
  temperature through the cooling ramp and the hold, c*(T) = exp(a - b / T) in mg per g of solvent, and
  J(S) = A S exp(-B / ln^2 S) with B = 16 pi gamma^3 v0^2 / (3 (k T)^3) (classical nucleation theory);
- the first nucleus becomes visible after a delay theta / (S - 1) evaluated at the nucleation instant, scattered
  log-normally from vial to vial; transmission then falls within a few seconds to a floor of 2-15 %;
- on the next heating the crystals dissolve; the transmission comes back half-way 20 s of heating after the vial
  passed its saturation temperature (logistic width 0.05 K); a vial that never nucleated stays clear.

Recorded: the transmission of every vial (0.01 % steps) and each block probe's reading, every 10 s.

Effects that only some runs contain (all physical; none is announced by a field of its own):
- `bubbles`: stirring fast enough to entrain air; each bubble crossing the beam is a dip of a few seconds;
- `loss_block`: a block whose vials lose solvent through their caps, raising c by `loss_per_cycle` each cycle;
- `dew`: condensation on the vials at a cold hold, a slow loss of transmission with no crystals in the vial;
- `pause`: (cycle, seconds) - the instrument stayed in that cycle's hold longer, shifting every later cycle;
- `block_volume_ml`: {block index: mL} - vials of a block filled to another volume than the protocol's;
- `offset_block`: a block whose probe reads `offset_K` below the block's true temperature (not used by any run).
"""
import math

import numpy as np

KB = 1.380649e-23
NA = 6.02214076e23
DISSOLVE_S = 20.0          # clear point = saturation temperature + heating rate x 20 s


def c_star(T, a, b):
    return np.exp(a - b / np.asarray(T, float))


def t_sat(c, a, b):
    return b / (a - math.log(c))


def barrier(gamma, v0, T):
    return 16 * math.pi * gamma ** 3 * v0 ** 2 / (3 * (KB * np.asarray(T, float)) ** 3)


def rate(S, A, B):
    S = np.asarray(S, float)
    B = np.broadcast_to(np.asarray(B, float), S.shape)
    out = np.zeros_like(S)
    m = S > 1.0 + 1e-12
    out[m] = A * S[m] * np.exp(-B[m] / np.log(S[m]) ** 2)
    return out


def molecular_volume(run):
    return run["molar_mass_g_mol"] * 1e-3 / (run["crystal_density_kg_m3"] * NA)


def solve_s(target_wait, A, B, V):
    """supersaturation at which the mean waiting time 1 / (V J) equals target_wait."""
    lo, hi = 1.0005, 6.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rate([mid], A, B)[0] * V * target_wait < 1:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def programme(run):
    """sample times, setpoint (K) and the timing of every cycle as the instrument ran it.

    A cycle heats from the hold temperature, holds at the top, cools and holds. `pause` = (cycle, seconds) keeps
    that cycle's hold going for the extra seconds, so every later cycle starts that much later.
    """
    heat = (run["T_top"] - run["T_hold"]) / (run["heat_K_min"] / 60)
    cool = (run["T_top"] - run["T_hold"]) / (run["cool_K_min"] / 60)
    pause_cycle, pause_s = run.get("pause", (None, 0.0))
    cycles, start = [], 0.0
    for k in range(run["cycles"]):
        hold = run["hold_s"] + (pause_s if k == pause_cycle else 0.0)
        c0 = start + heat + run["top_hold_s"]
        cycles.append(dict(start=start, cool=c0, hold=c0 + cool, end=c0 + cool + hold))
        start = c0 + cool + hold
    for k, cy in enumerate(cycles):
        cy["next_heat"] = cy["end"]
        cy["next_top"] = cy["end"] + heat
    n = int(round(start / run["dt_s"])) + 1
    t = np.arange(n) * run["dt_s"]
    sp = np.full(n, run["T_hold"])
    for cy in cycles:
        m = (t >= cy["start"]) & (t < cy["start"] + heat)
        sp[m] = run["T_hold"] + run["heat_K_min"] / 60 * (t[m] - cy["start"])
        sp[(t >= cy["start"] + heat) & (t < cy["cool"])] = run["T_top"]
        m = (t >= cy["cool"]) & (t < cy["hold"])
        sp[m] = run["T_top"] - run["cool_K_min"] / 60 * (t[m] - cy["cool"])
    return t, sp, cycles


def dew_factor(t, temp, run, rng):
    """water condensing on the outside of a vial colder than the room's dew point: transmission sinks slowly
    towards (1 - depth) while the block is below `above` K over the hold temperature, and clears within a minute
    once the block is warmer than that."""
    depth = run["dew"]["depth"] * rng.uniform(0.9, 1.1)
    tau, line = run["dew"]["tau_s"], run["T_hold"] + run["dew"]["above_K"]
    film = np.zeros(len(t))
    dt = np.diff(t, prepend=t[0])
    for i in range(1, len(t)):
        if temp[i] < line:
            film[i] = film[i - 1] + (1 - film[i - 1]) * dt[i] / tau
        else:
            film[i] = film[i - 1] * np.exp(-dt[i] / 40.0)
    return 1 - depth * film


def truths(run):
    v0 = molecular_volume(run)
    A = 10 ** run["log10_A"]
    B = barrier(run["gamma"], v0, run["T_hold"])
    s = run["S_ref"]
    return {"interfacial_energy_mj_m2": round(run["gamma"] * 1e3, 4),
            "log10_kinetic_factor": round(run["log10_A"], 4),
            "rate_at_reference_m3_s": float(f"{rate([s], A, B)[0]:.5g}"),
            "detection_delay_at_reference_s": round(run["theta"] / (s - 1), 3)}


def simulate(run, seed):
    rng = np.random.default_rng(seed)
    t, sp, cycles = programme(run)
    lag = run["heat_K_min"] / 60 * DISSOLVE_S
    nb, npos = 4, run["positions"]
    A = 10 ** run["log10_A"]
    v0 = molecular_volume(run)
    a, b = run["sol_a"], run["sol_b"]
    trans = np.zeros((len(t), nb * npos))
    probes = np.zeros((len(t), nb))
    tally = {"censored": 0, "before_hold": 0, "vial_cycles": 0}
    for j in range(nb):
        off = run.get("offset_K", 0.0) if run.get("offset_block") == j else 0.0
        true_T = sp + off
        V = run.get("block_volume_ml", {}).get(j, run["vial_ml"]) * 1e-6
        probes[:, j] = true_T - off + rng.normal(0, 0.02, len(t))
        for p in range(npos):
            base = 100.0 + np.cumsum(rng.normal(0, 0.004, len(t)))
            base -= np.linspace(0, base[-1] - 100.0, len(t))
            level = np.ones(len(t))
            for k in range(run["cycles"]):
                grow = run.get("loss_per_cycle", 0.0) if run.get("loss_block") == j else 0.0
                c = run["conc_mg_g"][j] * (1 + grow) ** k
                cy = cycles[k]
                c0, h0, h1 = cy["cool"], cy["hold"], cy["end"]
                fine = np.arange(c0, h1, 1.0)
                Tf = np.interp(fine, t, true_T)
                S = c / c_star(Tf, a, b)
                hz = rate(S, A, barrier(run["gamma"], v0, Tf)) * V
                H = np.cumsum(hz)
                tally["vial_cycles"] += 1
                draw = rng.exponential()
                if draw >= H[-1]:
                    tally["censored"] += 1
                    continue
                i = int(np.searchsorted(H, draw))
                td = fine[i] + run["theta"] / (S[i] - 1) * math.exp(rng.normal(0, run["delay_scatter"]))
                if td >= h1:
                    tally["censored"] += 1
                if td < h0:
                    tally["before_hold"] += 1
                end_heat = cy["next_top"]                    # end of the next heating ramp (or the record end)
                if td >= end_heat:
                    continue
                floor = rng.uniform(0.02, 0.15)
                turbid = 1.0 / (1 + np.exp(np.clip(-(t - td) / 2.0, -60, 60)))
                mid = t_sat(c, a, b) + lag
                dissolved = 1.0 / (1 + np.exp(np.clip(-(true_T - mid) / 0.05, -60, 60)))
                after = t >= cy["next_heat"]
                turbid = np.where(after, turbid * (1 - dissolved), turbid)
                turbid[t >= end_heat] = 0.0
                level = np.minimum(level, 1 - (1 - floor) * turbid)
            x = base * level
            if run.get("dew"):
                x = x * dew_factor(t, true_T, run, rng)
            x = x + rng.normal(0, 0.3, len(t))
            if run.get("bubbles", 0) > 0:
                for tb in rng.uniform(0, t[-1], rng.poisson(run["bubbles"] * t[-1])):
                    w = rng.uniform(2.0, 7.0)
                    x -= rng.uniform(15, 75) * np.exp(-0.5 * ((t - tb) / w) ** 2) * x / 100
            trans[:, j * npos + p] = x
    return t, np.clip(np.round(trans * 100), 0, 65535).astype(np.uint16), probes, tally
