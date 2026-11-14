"""Reduces one mass-comparator weighing to the test weight's true-mass and conventional-mass deviations.

    python3 mass_reduce.py <weighing directory> <result.json>

0. Sensitivity: the readings come in blocks of four (R, R+S, R+S, R three times; R T T R twelve times; R, R+S, R+S, R
   three times; the artefacts).  Each sensitivity block gives k = (mean I_R+S - mean I_R) / (m_S - rho_a V_S); k_start
   and k_end are the means of the blocks before and after the comparisons.  The comparator can recalibrate itself
   between two loadings, as often as its temperature calls for it, without a mark in the file.  Plans tried: no
   recalibration (only where k_start and k_end agree within 3e-4; k is then their mean), one recalibration between two
   R T T R blocks (only where they differ; k_start before it, k_end after it), and two recalibrations (any two of the
   11 places; k_start before the first, k_end after the second, the cycles between them, read at a sensitivity no
   block measures, left out).  Each plan is reduced by steps 1-4 and scored by its weighted sum of squares (the
   residuals of step 4 over (1.5 ug)^2) plus 60 per cycle left out and 30 per recalibration; the lowest score wins.
   Every indication is divided by the sensitivity in force at it (the artefacts and the end blocks: k_end).
1. Air density through the session: the CIPM formula on every line of the sensor log, interpolated to each reading's
   minute.  The sensors follow the changes of the room air but may be off as a whole, so the log is tied to the
   buoyancy artefacts: each artefact cycle (solid hollow hollow solid) gives d = I_hollow - I_solid =
   -dm - rho_h V_h + rho_s V_s, rho_h and rho_s the log's densities at the hollow and solid readings plus one common
   offset b; b is their least-squares value over the four cycles.  Every reading's air density is the log's plus b.
2. The reference as it is now: the check standard C keeps its certificate mass, so m_R = m_C + (R - C), the mean of
   the last four entries of the control chart; the reference's own certificate is not used (a cleaning after it
   changes the mass).
3. ABBA cycles (R T T R): per cycle y = mean(I_T) - mean(I_R) + mean(rho at the T readings) V_T -
   mean(rho at the R readings) V_R, which removes a linear comparator drift and the buoyancy of the air at the minutes
   the loads were on.
4. Temperature of the test weight: the probe reads dT = T_weight - T_air (0.05 K scatter) every 0.75 min; per
   cycle its readings interpolated to the two T readings are averaged, and the cycles y are fitted by linear least
   squares as y = d - c dT (the convection of a weight off the air's temperature), d the result; the residuals of
   this fit are the weighted sum of squares that places a recalibration (step 0).  (physics=True instead fits the
   first-order law of the weight's temperature with the logged air temperature to the cycles and the probe.)
5. The test weight's density rho_T from the laboratory's hydrostatic weighing (the sheet's value is the maker's
   nominal one): the readings are mass less the buoyancy of the medium, so V = (I_air - (I_water - I_pan)) /
   (rho_w - rho_a), m = I_air + rho_a V and rho_T = m / V, rho_w from the water's temperature (Tanaka formula) and
   rho_a from the air's conditions then (CIPM formula).
6. m_T = m_R + d, V_T = m_T / rho_T (iterated, since V_T enters step 3); the conventional mass is
   m_T (1 - 1.2 / rho_T) / (1 - 1.2 / 8000).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def cipm(t, p, rh):
    return (0.34848 * p - 0.009 * rh * np.exp(0.061 * t)) / (273.15 + t)


def weight_excess(log, tau, minutes, grid=0.05):
    """(exp(-t / tau), q(t)) at the given minutes: q is the excess of the weight over the air produced by the air's own
    changes alone, from d(dT)/dt = -dT / tau - dT_air/dt with dT(0) = 0, the air temperature taken from the log."""
    end = max(float(np.max(minutes)), 1.0) + grid
    t = np.arange(0.0, end + grid, grid)
    t_air = np.interp(t, log[:, 0], log[:, 1])
    q = np.zeros_like(t)
    decay = np.exp(-grid / tau)
    for i in range(1, t.size):                     # exact for an air temperature linear over each step
        slope = (t_air[i] - t_air[i - 1]) / grid
        q[i] = q[i - 1] * decay - slope * tau * (1 - decay)
    minutes = np.asarray(minutes, float)
    return np.exp(-minutes / tau), np.interp(minutes, t, q)


def water_density(t):
    """Tanaka et al. (2001), kg/m3, t in C."""
    a1, a2, a3, a4, a5 = -3.983035, 301.797, 522528.9, 69.34881, 999.974950
    return a5 * (1 - (t + a1) ** 2 * (t + a2) / (a3 * (t + a4)))


def hydrostatic_density(folder):
    hy = json.loads((Path(folder) / "hydrostatic.json").read_text())
    rho_a = cipm(hy["air"]["temperature_c"], hy["air"]["pressure_hpa"], hy["air"]["relative_humidity_pct"])
    r = hy["readings_g"]
    vol = (r["weight_in_air"] - (r["weight_on_pan_in_water"] - r["pan_in_water"])) / (
        (water_density(hy["water_temperature_c"]) - rho_a) / 1e3)                  # cm3
    return (r["weight_in_air"] + rho_a / 1e3 * vol) / vol * 1e3                  # kg/m3


def reduce_weighing(folder, air="anchored", chart=True, probe=True, gain="step", physics=False, plan=None,
                    density="hydrostatic"):
    """air: 'anchored' (the log tied to the artefacts), 'artefacts' (the artefacts' density for the whole weighing) or
    'sensors' (the log as it reads); chart=False takes the reference's certificate; probe=False ignores the weight's
    temperature; gain=False takes the comparator's sensitivity as 1, gain="mean" the mean of k_start and k_end,
    gain="linear" a sensitivity changing linearly in time between them, gain="single" at most one recalibration,
    searched only where k_start and k_end differ; physics=False uses the probe's readings,
    interpolated, as the regressor of a linear fit y = d - c dT; density='sheet' takes the sheet's density.  Options other than the defaults are the shortcuts
    tabulated by authoring/evidence/bench.py."""
    folder = Path(folder)
    sh = json.loads((folder / "weighing.json").read_text())
    with open(folder / "readings.csv", newline="") as fh:
        rows = [(float(r["minute"]), r["load"], float(r["indication_mg"])) for r in csv.DictReader(fh)]
    log = np.genfromtxt(folder / "air_log.csv", delimiter=",", skip_header=1)
    rho_log = cipm(log[:, 1], log[:, 2], log[:, 3])
    at = lambda minute: np.interp(minute, log[:, 0], rho_log)

    # 0. the comparator's sensitivity: blocks of four readings; the sensitivity blocks at the start and at the end
    blocks = [rows[i:i + 4] for i in range(0, len(rows), 4)]
    sw = sh["sensitivity_weight"]
    v_sw = sw["mass_mg"] / 1e3 / (sw["density_kg_m3"] / 1e3)

    def k_of(block):
        plus = [x for x in block if x[1] == "R+S"]
        base = [x for x in block if x[1] == "R"]
        return ((np.mean([x[2] for x in plus]) - np.mean([x[2] for x in base]))
                / (sw["mass_mg"] - np.mean([at(x[0]) for x in plus]) * v_sw))

    kind = ["S" if any(x[1] == "R+S" for x in bl) else "T" if any(x[1] == "T" for x in bl) else "A" for bl in blocks]
    first_t = kind.index("T")
    k_start = float(np.mean([k_of(bl) for bl, kd in zip(blocks[:first_t], kind[:first_t]) if kd == "S"]))
    k_end = float(np.mean([k_of(bl) for bl, kd in zip(blocks[first_t:], kind[first_t:]) if kd == "S"]))
    t_blocks = [i for i, kd in enumerate(kind) if kd == "T"]
    same = abs(k_start - k_end) <= 3e-4
    if plan is None and gain in ("step", "single"):
        n = len(t_blocks)
        plans = [(0, 0)] if same else [(j, j) for j in range(1, n)]
        if gain == "step":
            plans += [(a, b) for a in range(1, n) for b in range(a + 1, n)]
        best = None
        for pl in plans:
            out = reduce_weighing(folder, air, chart, probe, gain, physics, plan=pl, density=density)
            score = out.pop("_cost") + 60.0 * (pl[1] - pl[0]) + 30.0 * ((pl[0] > 0) + (pl[1] > pl[0]))
            if best is None or score < best[0]:
                best = (score, out)
        return best[1]
    a_cut, b_cut = plan if plan else (0, 0)

    def k_block(i):
        if not gain:
            return 1.0
        if gain == "mean" or (gain in ("step", "single") and same):
            return 0.5 * (k_start + k_end)
        if i < t_blocks[0]:
            return k_start
        if i > t_blocks[-1]:
            return k_end
        if gain == "linear":
            u = (i - t_blocks[0] + 0.5) / len(t_blocks)
            return k_start + u * (k_end - k_start)
        return k_start if t_blocks.index(i) < a_cut else k_end

    rows = [(m, lab, ind / k_block(i // 4)) for i, (m, lab, ind) in enumerate(rows)]

    art = sh["artefacts"]
    v_s = art["solid_mass_g"] / (art["solid_density_kg_m3"] / 1e3)
    v_h, dm = art["hollow_volume_cm3"], art["true_mass_difference_mg"]
    arts = [r for r in rows if r[1].startswith("art-")]
    lhs, rhs = [], []
    for c in range(0, len(arts), 4):
        block = arts[c:c + 4]
        h = [x for x in block if x[1] == "art-hollow"]
        s = [x for x in block if x[1] == "art-solid"]
        d = np.mean([x[2] for x in h]) - np.mean([x[2] for x in s])
        rh_, rs_ = np.mean([at(x[0]) for x in h]), np.mean([at(x[0]) for x in s])
        # d = -dm - (rh_ + b) V_h + (rs_ + b) V_s  ->  b (V_s - V_h) = d + dm + rh_ V_h - rs_ V_s
        lhs.append(v_s - v_h)
        rhs.append(d + dm + rh_ * v_h - rs_ * v_s)
    lhs, rhs = np.array(lhs), np.array(rhs)
    b = float(lhs @ rhs / (lhs @ lhs))
    if air == "anchored":
        rho = lambda minute: at(minute) + b
    elif air == "artefacts":
        rho_end = float(np.mean([at(x[0]) for x in arts])) + b
        rho = lambda minute: rho_end
    else:
        rho = at

    nom = sh["nominal_g"] * 1e3
    if chart:
        with open(folder / "control_chart.csv", newline="") as fh:
            r_minus_c = [float(r["r_minus_c_mg"]) for r in csv.DictReader(fh)]
        m_r = nom + sh["check_standard"]["mass_deviation_mg"] + float(np.mean(r_minus_c[-4:]))
    else:
        m_r = nom + sh["reference"]["mass_deviation_mg"]
    temp_log = np.genfromtxt(folder / "weight_temperature.csv", delimiter=",", skip_header=1)
    v_r = (m_r / 1e3) / (sh["reference"]["density_kg_m3"] / 1e3)
    rho_t = sh["test_weight"]["density_kg_m3"] if density == "sheet" else hydrostatic_density(folder)
    main = [r for i, r in enumerate(rows) if kind[i // 4] == "T"   # the R T T R cycles outside a pair of steps
            and not a_cut <= t_blocks.index(i // 4) < b_cut]
    m_t = m_r
    fit_state = None
    for _ in range(4):
        v_t = (m_t / 1e3) / (rho_t / 1e3)
        ys, t_minutes = [], []
        for c in range(0, len(main), 4):
            block = main[c:c + 4]
            tr = [x for x in block if x[1] == "T"]
            rr = [x for x in block if x[1] == "R"]
            ys.append(np.mean([x[2] for x in tr]) - np.mean([x[2] for x in rr])
                      + np.mean([rho(x[0]) for x in tr]) * v_t - np.mean([rho(x[0]) for x in rr]) * v_r)
            t_minutes.append([x[0] for x in tr])
        ys = np.array(ys)
        d_rest = float(np.mean(ys))
        cost = float(np.sum((ys - d_rest) ** 2) / 0.0015 ** 2)
        if probe and not physics:
            dts = np.array([np.mean(np.interp(tm, temp_log[:, 0], temp_log[:, 1])) for tm in t_minutes])
            design = np.column_stack([np.ones_like(dts), -dts])
            coef, *_ = np.linalg.lstsq(design, ys, rcond=None)
            d_rest = float(coef[0])
            cost = float(np.sum((ys - design @ coef) ** 2) / 0.0015 ** 2)
        elif probe:
            flat_t = np.array(t_minutes).ravel()

            def parts(tau):
                e_t, q_t = weight_excess(log, tau, flat_t)
                e_p, q_p = weight_excess(log, tau, temp_log[:, 0])
                return e_t.reshape(-1, 2).mean(axis=1), q_t.reshape(-1, 2).mean(axis=1), e_p, q_p

            def resid(x):
                d, c, dt0, tau = x
                e_c, q_c, e_p, q_p = parts(tau)
                r_cycles = (ys - (d - c * (dt0 * e_c + q_c))) / 0.0015
                r_probe = (temp_log[:, 1] - (dt0 * e_p + q_p)) / 1.5
                return np.r_[r_cycles, r_probe]

            if fit_state is None:                  # start: the best of a grid of time constants
                best = None
                for tau0 in (10.0, 20.0, 35.0, 60.0, 100.0):
                    e_c, q_c, e_p, q_p = parts(tau0)
                    dt0 = float(np.mean(temp_log[:3, 1]))
                    a = np.column_stack([np.ones_like(ys), -(dt0 * e_c + q_c)])
                    coef, *_ = np.linalg.lstsq(a, ys, rcond=None)
                    x0 = np.r_[coef[0], np.clip(coef[1], 0.0, 0.2), dt0, tau0]
                    cost = float(np.sum(resid(x0) ** 2))
                    if best is None or cost < best[0]:
                        best = (cost, x0)
                fit_state = best[1]
            sol = least_squares(resid, fit_state, bounds=([-np.inf, 0.0, -15.0, 3.0], [np.inf, 0.3, 15.0, 300.0]),
                                x_scale="jac")
            fit_state = sol.x
            d_rest = float(sol.x[0])
            cost = float(2 * sol.cost)
        m_t = m_r + d_rest
    conv = m_t * (1 - 1.2 / rho_t) / (1 - 1.2 / 8000.0)
    out = {"true_mass_deviation_mg": float(m_t - nom), "conventional_mass_deviation_mg": float(conv - nom)}
    if plan is not None:
        out["_cost"] = cost
    return out


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_weighing(sys.argv[1])))
