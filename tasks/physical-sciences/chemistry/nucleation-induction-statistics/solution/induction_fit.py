"""Reference reduction of one parallel-crystalliser run.

    python3 induction_fit.py <run folder> <answer file>

1. Cycle timing comes from the probe trace (the median of the four blocks), not from the programme entries: each
   stretch at the dissolution temperature opens a cycle, and the following stretch at the hold temperature closes it.
   In every vial and cycle, the detection is the first sample after which the transmission stays, until the hold
   ends, below half of the highest reading of the 60 s before it: a fall within seconds, not a slow drift or a
   passing dip. A vial without such a fall by the end of the hold is censored there.
2. The vial contents come from the clear points. On each heating, the probe temperature where a turbid vial is back
   half-way, less (heating rate x 20 s), is its saturation temperature and gives c = c*(T). Per block, ln c is
   fitted linearly against the cycle number, so solvent creeping out of a block enters S(t) = c / c*(T_probe(t))
   through the vials' own saturation temperatures rather than through the nominal concentration.
3. Vial volume: the protocol's value, unless the operator's log records that a block was filled to another volume.
4. Likelihood. For each block and cycle the hazard V J(S(t)) runs on a 2 s grid from the top of the cycle to the end
   of the hold. Each nucleation instant carries a log-normal detection delay with median theta / (S - 1) and spread
   sigma. The probability of a detection inside its sample interval, or of none before the hold ends, gives the
   likelihood. The parameters ln J(S_ref), B, ln theta and ln sigma are found by Nelder-Mead with restarts, and
   B scales as (T_hold / T)^3 off the hold temperature.
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr

KB = 1.380649e-23
NA = 6.02214076e23
GRID = 2.0
DISSOLVE_S = 20.0


def load(folder):
    proto = json.loads((folder / "protocol.json").read_text())
    trans = np.load(folder / "transmission.npy").astype(float) / 100.0
    probes = np.loadtxt(folder / "block_probes.csv", delimiter=",", skiprows=1)
    return proto, trans, probes[:, 0], probes[:, 1:] + 273.15


def stretches(flag):
    """(first, last) index of every run of True in a boolean array."""
    edges = np.diff(np.concatenate([[0], flag.astype(int), [0]]))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1) - 1))


def segments(proto, t, temps):
    pg = proto["programme"]
    T_med = np.median(temps, axis=1)
    top = stretches(T_med > pg["dissolution_temperature_C"] + 273.15 - 0.3)
    hold = stretches(T_med < pg["hold_temperature_C"] + 273.15 + 0.1)
    out = []
    for k, (a0, a1) in enumerate(top):
        later = [h for h in hold if h[0] > a1]
        if not later:
            break
        h0, h1 = later[0]
        nxt = top[k + 1][0] if k + 1 < len(top) else len(t) - 1
        out.append(dict(cool=t[a1], hold=t[h0], end=t[h1], next_heat=t[h1], next_top=t[nxt]))
    return out


def nominal_segments(proto):
    """evidence switch only: cycle timing taken from the programme entries, as if no cycle ever ran long."""
    pg = proto["programme"]
    span = pg["dissolution_temperature_C"] - pg["hold_temperature_C"]
    heat, cool = span / (pg["heating_rate_K_min"] / 60), span / (pg["cooling_rate_K_min"] / 60)
    cyc = heat + pg["dissolution_hold_s"] + cool + pg["hold_s"]
    return [dict(cool=k * cyc + heat + pg["dissolution_hold_s"], hold=k * cyc + heat + pg["dissolution_hold_s"] + cool,
                 end=(k + 1) * cyc, next_heat=(k + 1) * cyc, next_top=(k + 1) * cyc + heat) for k in range(pg["cycles"])]


def observe(proto, trans, t, temps, persist=True, timing="probes", detector="relative"):
    segs = segments(proto, t, temps) if timing == "probes" else nominal_segments(proto)
    names = sorted(proto["blocks"])
    events = {}              # (block, cycle) -> list of (detected?, lo, hi)
    clear = {}               # (block, cycle) -> list of clear-point temperatures (K)
    for j, name in enumerate(names):
        for p in range(len(proto["blocks"][name]["positions"])):
            x = trans[:, 8 * j + p]
            prior = np.concatenate([np.full(6, x[0]), np.lib.stride_tricks.sliding_window_view(x, 6).max(axis=1)])
            for k, s in enumerate(segs):
                idx = np.flatnonzero((t >= s["cool"]) & (t <= s["end"]))
                ahead = np.flip(np.maximum.accumulate(np.flip(x[idx])))
                level = 0.5 * prior[idx] if detector == "relative" else 50.0
                drop = (ahead if persist else x[idx]) < level
                hit = np.flatnonzero(drop)
                if hit.size:
                    i = idx[hit[0]]
                    events.setdefault((j, k), []).append((1, t[i - 1], t[i]))
                else:
                    events.setdefault((j, k), []).append((0, s["end"], s["end"]))
                    continue
                heat = np.flatnonzero((t >= s["next_heat"]) & (t < s["next_top"]))
                if heat.size < 5:
                    continue
                turbid = np.median(x[idx[-5:]])
                clear_level = np.median(x[heat[-5:]])
                mid = 0.5 * (turbid + clear_level)
                up = np.flatnonzero(x[heat] > mid)
                if up.size and up[0] > 0:
                    i = heat[up[0]]
                    f = (mid - x[i - 1]) / (x[i] - x[i - 1])
                    clear.setdefault((j, k), []).append(temps[i - 1, j] + f * (temps[i, j] - temps[i - 1, j]))
    return segs, names, events, clear


def contents(proto, names, clear, n, use_clear=True, lag_corrected=True):
    a, b = proto["solubility"]["a"], proto["solubility"]["b_K"]
    lag = proto["programme"]["heating_rate_K_min"] / 60 * DISSOLVE_S if lag_corrected else 0.0
    conc = np.zeros((len(names), n))
    for j, name in enumerate(names):
        ks = [k for k in range(n) if (j, k) in clear]
        if use_clear and len(ks) >= 3:
            ln_c = [a - b / (np.median(clear[(j, k)]) - lag) for k in ks]
            slope, icept = np.polyfit(ks, ln_c, 1)
            conc[j] = np.exp(icept + slope * np.arange(n))
        else:
            conc[j] = proto["blocks"][name]["concentration_mg_g"]
    return conc


def volumes(proto, folder, names):
    """m^3 per block: the protocol's vial volume unless the log names another fill for a block."""
    vol = {name: proto["vial_volume_ml"] for name in names}
    log = folder / "run_log.txt"
    if log.exists():
        for name, ml in re.findall(r"block ([A-D]) vials filled to ([0-9.]+) mL", log.read_text()):
            vol[name] = float(ml)
    return [vol[name] * 1e-6 for name in names]


def reduce(folder, persist=True, censored="censor", use_clear=True, clock="cooling", timing="probes",
           log_volumes=True, lag_corrected=True, detector="relative"):
    """The keyword switches exist for the evidence scripts; the defaults are the reduction described above."""
    proto, trans, t, temps = load(folder)
    segs, names, events, clear = observe(proto, trans, t, temps, persist, timing, detector)
    conc = contents(proto, names, clear, len(segs), use_clear, lag_corrected)
    vols = volumes(proto, folder, names) if log_volumes else [proto["vial_volume_ml"] * 1e-6] * len(names)
    a, b = proto["solubility"]["a"], proto["solubility"]["b_K"]
    T_hold = proto["programme"]["hold_temperature_C"] + 273.15
    s_ref = proto["reference_supersaturation"]
    v0 = proto["molar_mass_g_mol"] * 1e-3 / (proto["crystal_density_kg_m3"] * NA)

    pieces = []
    for (j, k), rows in sorted(events.items()):
        s = segs[k]
        if censored == "drop":
            rows = [r for r in rows if r[0] == 1]
        if clock == "hold":
            rows = [(r[0], max(r[1], s["hold"]), max(r[2], s["hold"] + GRID)) if r[0] else r for r in rows]
        grid = np.arange(s["cool"] if clock == "cooling" else s["hold"], s["end"], GRID)
        T = np.interp(grid, t, temps[:, j])
        S = conc[j, k] / np.exp(a - b / T)
        live = S > 1.0005
        if not live.any():
            continue
        grid, T, S = grid[live], T[live], S[live]
        kind = np.array([r[0] for r in rows])
        lo = np.array([r[1] for r in rows])
        hi = np.array([r[2] for r in rows])
        query = np.concatenate([lo[kind == 1], hi[kind == 1], [s["end"]]])
        lag = query[:, None] - grid[None, :]
        pieces.append(dict(V=vols[j], lnS=np.log(S), S=S, cube=(T_hold / T) ** 3, n_det=int(kind.sum()),
                           n_cens=int((kind == 0).sum()), loglag=np.log(np.maximum(lag, 1e-9)), pos=lag > 0))

    def nll(q):
        u, B, ln_th, ln_sig = q
        if B <= 0 or not -6.0 < ln_sig < -0.5:
            return 1e12
        lnA = u - math.log(s_ref) + B / math.log(s_ref) ** 2
        sig = math.exp(ln_sig)
        total = 0.0
        for pc in pieces:
            hz = pc["V"] * np.exp(lnA + pc["lnS"] - B * pc["cube"] / pc["lnS"] ** 2)
            H = np.concatenate([[0.0], np.cumsum(hz[:-1]) * GRID])
            mass = np.exp(-H) * -np.expm1(-hz * GRID)
            z = (pc["loglag"] - (ln_th - np.log(pc["S"] - 1))[None, :]) / sig
            F = np.where(pc["pos"], ndtr(z), 0.0) @ mass
            n = pc["n_det"]
            p_det = F[n:2 * n] - F[:n]
            p_cens = 1.0 - F[-1]
            total -= np.log(np.maximum(p_det, 1e-300)).sum() + pc["n_cens"] * math.log(max(p_cens, 1e-300))
        return total

    best = None
    for start in ([math.log(3e3), 0.4, math.log(80.0), -3.3], [math.log(1e4), 0.25, math.log(150.0), -3.0],
                  [math.log(1e3), 0.6, math.log(50.0), -3.6]):
        r = minimize(nll, start, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=3000))
        for _ in range(2):
            r = minimize(nll, r.x, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=3000))
        if best is None or r.fun < best.fun:
            best = r
    u, B, ln_th, ln_sig = best.x
    lnA = u - math.log(s_ref) + B / math.log(s_ref) ** 2
    gamma = (3 * B * (KB * T_hold) ** 3 / (16 * math.pi * v0 ** 2)) ** (1 / 3)
    return {"interfacial_energy_mj_m2": gamma * 1e3, "log10_kinetic_factor": lnA / math.log(10),
            "rate_at_reference_m3_s": math.exp(u), "detection_delay_at_reference_s": math.exp(ln_th) / (s_ref - 1)}


if __name__ == "__main__":
    answer = reduce(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(answer, indent=1) + "\n")
