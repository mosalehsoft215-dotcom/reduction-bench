"""Physics of one bipolar membrane chlor-alkali electrolyser line during a guarantee run.

Everything the records show is computed from a single time-stepped plant state (10 s steps):
the rectifier current, the catholyte loop mass and NaOH fraction under level and concentration
control, the element voltages, the busbar, the rectifier and its AC meter, and the hydrogen
leaving the cooler. Instruments then read that state with their own errors.
"""
import math

import numpy as np

FARADAY = 96485.33212        # C/mol
M_NAOH = 39.99711e-3         # kg/mol
M_H2O = 18.01528e-3          # kg/mol
V_NORMAL = 22.41397e-3       # m3/mol of ideal gas at 0 degC, 101.325 kPa
G_ACC = 9.80665
STEP = 10.0                  # s


def caustic_density(w, t_c):
    """Density (kg/m3) of aqueous NaOH, mass fraction w (0.20-0.40), temperature t_c (degC)."""
    return 981.4 + 1255.5 * w - 335.0 * w * w - (0.36 + 0.75 * w) * (t_c - 20.0)


def caustic_fraction(rho, t_c):
    """Inverse of caustic_density for the mass fraction."""
    w = 0.32 + 0.0 * np.asarray(rho, dtype=float)
    for _ in range(30):
        f = caustic_density(w, t_c) - rho
        df = 1255.5 - 670.0 * w - 0.75 * (t_c - 20.0)
        w = w - f / df
    return w


def water_vapour_kpa(t_c):
    """Saturation pressure of water (kPa)."""
    return 0.1333224 * 10.0 ** (8.07131 - 1730.63 / (233.426 + t_c))


def busbar_resistance(bus, t_bus):
    """Resistance (ohm) of the two busbar poles at conductor temperature t_bus, joints included."""
    area_mm2 = bus["bars_per_pole"] * bus["bar_width_mm"] * bus["bar_thick_mm"]
    rho = bus["resistivity_20"] * (1.0 + bus["temp_coeff"] * (t_bus - 20.0))   # ohm mm2/m
    return rho * 2.0 * bus["route_m"] / area_mm2 + bus["joints"] * bus["joint_uohm_true"] * 1e-6


def rectifier_loss_kw(rect, i_ka):
    return rect["p0_kw"] + rect["a_v"] * i_ka + rect["b_mohm"] * i_ka * i_ka


def current_profile(d, t):
    """Set point (kA) at times t (s from the first record) and the ramp-limited actual current."""
    sp = np.full_like(t, d["i_set_ka"])
    for ev in d["setpoints"]:
        sp = np.where(t >= ev[0] * 3600.0, ev[1], sp)
    i = np.empty_like(sp)
    i[0] = sp[0]
    ramp = d["ramp_ka_per_min"] * STEP / 60.0
    for k in range(1, len(sp)):
        i[k] = i[k - 1] + max(-ramp, min(ramp, sp[k] - i[k - 1]))
    return sp, i


def simulate(d, rng):
    """Run the line for d['record_hours'] hours; returns a dict of 10-s series and the true results."""
    n_pre = int(round(d["pre_roll_h"] * 3600.0 / STEP))
    n_steps = n_pre + int(round(d["record_hours"] * 3600.0 / STEP))
    t = (np.arange(n_steps) - n_pre) * STEP
    hours = t / 3600.0
    sp, i_true = current_profile(d, t)
    i_true = i_true * (1.0 + 4e-5 * rng.standard_normal(n_steps))
    area = d["area_m2"]
    n_el = d["n_el"]
    j = i_true / area                                   # kA/m2

    # temperatures
    load = i_true / d["i_set_ka"]
    lag_t = np.empty(n_steps)
    lag_b = np.empty(n_steps)
    lag_t[0] = load[0]
    lag_b[0] = load[0] ** 2
    for k in range(1, n_steps):
        lag_t[k] = lag_t[k - 1] + (load[k] - lag_t[k - 1]) * STEP / 1500.0
        lag_b[k] = lag_b[k - 1] + (load[k] ** 2 - lag_b[k - 1]) * STEP / 2400.0
    phase = rng.uniform(0, 2 * math.pi)
    t_cell = (d["t_cell_set"] - 7.0 * (1.0 - lag_t) + 0.25 * np.sin(2 * math.pi * hours / 6.3 + phase)
              + 0.05 * rng.standard_normal(n_steps))
    t_amb = d["t_amb"] + 4.5 * np.sin(2 * math.pi * (hours + d["clock_hour0"] - 9.0) / 24.0)
    t_bus = t_amb + d["bus"]["rise_full_k"] * lag_b * (d["i_set_ka"] / d["rect"]["rated_ka"]) ** 2
    t_prod = d["t_prod"] + 0.3 * np.sin(2 * math.pi * hours / 4.1 + 2 * phase) + 0.05 * rng.standard_normal(n_steps)
    t_cool = d["t_cool"] + 0.3 * np.sin(2 * math.pi * (hours + d["clock_hour0"] - 10.0) / 24.0) \
        + 0.05 * rng.standard_normal(n_steps)
    p_h2 = d["p_h2_kpa"] + 0.3 * np.sin(2 * math.pi * hours / 2.7 + phase) + 0.05 * rng.standard_normal(n_steps)

    # catholyte loop under level (product valve) and concentration (DM water) control
    tank_area = math.pi * d["tank_diameter_m"] ** 2 / 4.0
    rho_cal = d["lt_rho_cal"]
    w = d["w_start"]
    m_hold = d["holdup_m3"] * caustic_density(w, d["t_cell_set"])
    m_tank = d["level_sp_mm"] / 1000.0 * tank_area * rho_cal
    m_loop = m_tank + m_hold
    trims = sorted(d["trims"])
    lots = sorted(d["lab_lots"])
    out = {k: np.empty(n_steps) for k in ("w", "m_tank", "m_hold", "gen", "prod", "dm", "ce", "u_cell",
                                          "level_ind", "w_an", "trim")}
    lic_i = 0.0
    aic_i = 0.0
    gen0 = d["ce_ref"] * n_el * d["i_set_ka"] * 1e3 / FARADAY
    prod_bias = gen0 * M_NAOH / d["w_set"]
    dm_bias = d["dm_bias_kg_s"]
    pump_off = d.get("pump_stop")
    dm_manual = d.get("dm_manual")
    noise_lvl = 0.0004 * rng.standard_normal(n_steps)
    noise_w = 0.0001 * rng.standard_normal(n_steps)
    dm_hold = None
    for k in range(n_steps):
        hk = hours[k]
        ce = d["ce_ref"] + d["dce_dj"] * (j[k] - d["j_ref"]) + d["dce_dw"] * (100.0 * w - 32.0)
        n_e = n_el * i_true[k] * 1e3 / FARADAY          # mol e-/s
        gen = ce * n_e                                   # mol NaOH/s
        y_w = water_vapour_kpa(t_cool[k]) / p_h2[k]
        vap = 0.5 * n_e * y_w / (1.0 - y_w)             # mol water/s leaving with the hydrogen
        water = (d["beta_w"] * gen - n_e - vap) * M_H2O  # kg/s from transport, cathode use, vapour
        level_ind = (m_tank / (tank_area * rho_cal)) * 1000.0 + noise_lvl[k] * 1000.0
        rho_p = caustic_density(w, t_prod[k])
        trim = 0.0
        for tt, val in trims:
            if hk >= tt:
                trim = val
        w_an = 100.0 * caustic_fraction(rho_p + d["dt311_bias"] + 100.0 * noise_w[k], t_prod[k]) + trim
        # product valve: PI on the indicated level
        err_l = (level_ind - d["level_sp_mm"]) / 1000.0
        lic_i = min(max(lic_i + 0.012 * err_l * STEP, -2.0 * prod_bias), 2.0 * prod_bias)
        prod = min(max(prod_bias + 12.0 * err_l + lic_i, 0.0), 3.0 * prod_bias)
        if pump_off and pump_off[0] <= hk < pump_off[1]:
            prod = 0.0
        # DM water: PI on the analyser
        err_w = w_an - d["w_set"] * 100.0
        if dm_manual and dm_manual[0] <= hk < dm_manual[1]:
            if dm_hold is None:
                dm_hold = dm
            dm = dm_hold
        else:
            dm_hold = None
            aic_i = min(max(aic_i + 0.0006 * err_w * STEP, -dm_bias), 2.0 * dm_bias)
            dm = min(max(dm_bias + 0.6 * err_w + aic_i, 0.0), 3.0 * dm_bias)
        out["w"][k] = w
        out["m_tank"][k] = m_tank
        out["m_hold"][k] = m_hold
        out["gen"][k] = gen
        out["prod"][k] = prod
        out["dm"][k] = dm
        out["ce"][k] = ce
        out["level_ind"][k] = level_ind
        out["w_an"][k] = w_an
        out["trim"][k] = trim
        # integrate the loop
        naoh = m_loop * w + (gen * M_NAOH - prod * w) * STEP
        m_loop = m_loop + (gen * M_NAOH + water + dm - prod) * STEP
        w = naoh / m_loop
        m_hold = d["holdup_m3"] * caustic_density(w, t_cell[k])
        m_tank = m_loop - m_hold

    # element voltages (V) and the electrolyser terminal voltage
    u_mean = (d["u0"] + d["k_j"] * j + d["c_t"] * (t_cell - 88.0) + d["c_w"] * (100.0 * out["w"] - 32.0)
              + d["u_drift_per_h"] * hours)
    spread = d["elem_spread"] * rng.standard_normal(n_el)
    spread -= spread.mean()
    u_el = u_mean[:, None] + spread[None, :] * (j / d["j_ref"])[:, None]
    u_el += 0.0004 * rng.standard_normal(u_el.shape)
    u_e = u_el.sum(axis=1)
    r_bus = busbar_resistance(d["bus"], t_bus)
    v_r = u_e + r_bus * i_true * 1e3
    loss = rectifier_loss_kw(d["rect"], i_true) * (1.0 + d["rect"]["site_dev"])
    p_ac = (v_r * i_true * 1e3 / 1e6) + loss / 1e3      # MW
    n_h2 = 0.5 * n_el * i_true * 1e3 / FARADAY           # mol/s dry
    y_w = water_vapour_kpa(t_cool) / p_h2
    q_act = n_h2 / (1.0 - y_w) * 8.314462618 * (t_cool + 273.15) / (p_h2 * 1e3)   # m3/s wet at the meter

    out.update(dict(t=t, hours=hours, sp=sp, i=i_true, j=j, t_cell=t_cell, t_bus=t_bus, t_prod=t_prod,
                    t_cool=t_cool, p_h2=p_h2, u_el=u_el, u_e=u_e, v_r=v_r, p_ac=p_ac, n_h2=n_h2,
                    q_act=q_act, r_bus=r_bus))
    keep = slice(n_pre, n_steps)
    for key, val in list(out.items()):
        if isinstance(val, np.ndarray) and val.shape[0] == n_steps:
            out[key] = val[keep]
    out["tank_area"] = tank_area
    hours = out["hours"]
    i_true = out["i"]
    u_e = out["u_e"]

    # true results over the run window
    t0, t1 = d["window_h"]
    sel = (hours >= t0 - 1e-9) & (hours < t1 - 1e-9)
    naoh_kg = np.sum(out["gen"][sel]) * M_NAOH * STEP
    charge_mol = np.sum(n_el * i_true[sel] * 1e3 / FARADAY) * STEP
    energy_j = np.sum(u_e[sel] * i_true[sel] * 1e3) * STEP
    out["truth"] = {
        "current_efficiency_pct": 100.0 * naoh_kg / M_NAOH / charge_mol,
        "dc_energy_kwh_per_t": energy_j / naoh_kg / 3600.0,
    }
    return out
