"""Physical model behind the synthetic chromatography records (generator side).

Flow path of one run, all at the run's volumetric flow Q:
  system pump (buffers A and B blended at %B) -> mixer (one well-mixed volume) -> tubing to the injection valve
  -> [during loading the sample pump feeds sample through the valve instead] -> tubing to the column
  -> column (N equal well-mixed stages, kinetic steric mass action) -> tubing and detector cells in the order the
  system file gives. Tubing and cells are plug flow.
Units inside: mL, s, mM (mmol/L), K.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

ALPHA_NACL = 0.0212  # 1/K, conductivity temperature coefficient of NaCl solutions about 25 C


def tube_volume(length_cm, id_mm):
    return np.pi * (id_mm / 20.0) ** 2 * length_cm  # mL


def column_volumes(col):
    Vc = np.pi * (col["diameter_cm"] / 2) ** 2 * col["bed_height_cm"]
    N = col["stages"]
    return dict(Vc=Vc, N=N, Vl=col["porosity_salt"] * Vc / N, Vlp=col["porosity_protein"] * Vc / N,
                Vs=(1 - col["porosity_salt"]) * Vc / N)


def _sparsity(N):
    S = lil_matrix((3 * N, 3 * N))
    for i in range(N):
        for a in range(3):
            for b in range(3):
                S[a * N + i, b * N + i] = 1
        if i:
            S[i, i - 1] = 1
            S[N + i, N + i - 1] = 1
    return S.tocsr()


def simulate_column(col, sma, Qfun, t_grid, salt_in, prot_in, t_end):
    """salt_in, prot_in: column inlet on t_grid; Qfun(t) the flow (mL/s). Returns outlet salt and protein."""
    v = column_volumes(col)
    N, Vl, Vlp, Vs = v["N"], v["Vl"], v["Vlp"], v["Vs"]
    Lam, kk = col["ionic_capacity_mM"], sma["k_kin_per_s"]
    nu, keq, sig = sma["nu"], sma["k_eq"], sma["sigma"]

    def rhs(t, y):
        Q = Qfun(t)
        cs, cp, qp = y[:N], y[N:2 * N], y[2 * N:]
        qbar = np.maximum(Lam - (nu + sig) * qp, 0.0)
        rate = kk * (cp * (qbar / Lam) ** nu - (qp / keq) * (np.maximum(cs, 1e-12) / Lam) ** nu)
        s0 = np.interp(t, t_grid, salt_in)
        p0 = np.interp(t, t_grid, prot_in)
        cs_up = np.concatenate([[s0], cs[:-1]])
        cp_up = np.concatenate([[p0], cp[:-1]])
        return np.concatenate([(Q * (cs_up - cs) + Vs * nu * rate) / Vl, (Q * (cp_up - cp) - Vs * rate) / Vlp, rate])

    y0 = np.concatenate([np.full(N, salt_in[0]), np.zeros(2 * N)])
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", jac_sparsity=_sparsity(N), rtol=1e-7, atol=1e-11,
                    max_step=1.0, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    Y = sol.sol(t_grid)
    return Y[N - 1], Y[2 * N - 1]


def well_mixed(t, u, V, Q, before):
    """Outlet of one well-mixed volume V that held `before` at t = 0, for the inlet u(t) (exact for u piecewise
    linear between grid points); Q[i] is the flow during the step that ends at t[i]."""
    out = np.empty_like(u)
    out[0] = before
    dt = t[1] - t[0]
    for i in range(1, len(t)):
        a = np.exp(-Q[i] * dt / V)
        g = (V / (Q[i] * dt)) * (1 - a)
        u0 = u[i - 1] if i > 1 else before
        out[i] = u[i] + (out[i - 1] - u0) * a - (u[i] - u0) * g
    return out


def delayed(vol_pumped, x, volume, before):
    """Plug flow through `volume` that held `before` at t = 0: the outlet now is what entered when the pumped
    volume was smaller by `volume` (valid for any flow history)."""
    return np.interp(vol_pumped - volume, vol_pumped, x, left=before)


def run_records(system, col, sma, buffers, run, protein, dt=0.1, method_Vc=None):
    """Simulate one run; return the time grid and the true salt (sodium) and protein along the flow path.
    method_Vc: column volume the method converts CV with (the system file's), which stays the same when the bed
    itself changes; run["B_delivery"] (default 1): fraction of the set %B the pump actually delivers."""
    Qe = run["flow_mL_min"] / 60.0
    Ql = run["load_flow_mL_min"] / 60.0
    Vc = column_volumes(col)["Vc"] if method_Vc is None else method_Vc
    t_load = run["load_volume_mL"] / Ql
    t0g = t_load + run["wash_CV"] * Vc / Qe
    t_grad = run["gradient_CV"] * Vc / Qe
    t_end = t0g + t_grad + run["hold_CV"] * Vc / Qe
    t = np.arange(0.0, t_end + dt / 2, dt)
    Q = np.where(t <= t_load + 1e-9, Ql, Qe)
    vol = Ql * np.minimum(t, t_load) + Qe * np.maximum(t - t_load, 0.0)

    def Qfun(tt):
        return Ql if tt < t_load else Qe

    pctB = np.where(t < t0g, run["start_pctB"],
                    np.where(t < t0g + t_grad, run["start_pctB"] + (run["end_pctB"] - run["start_pctB"]) * (t - t0g) / t_grad,
                             run["end_pctB"]))
    na = buffers.get("Na_from_buffer_mM", 0.0)
    A = buffers["A_NaCl_mM"] + na
    pump_salt = A + (buffers["B_NaCl_mM"] - buffers["A_NaCl_mM"]) * pctB * run.get("B_delivery", 1.0) / 100.0
    c_sample = run["sample_mg_mL"] / protein["molar_mass_g_mol"] * 1000.0
    s_sample = run["sample_NaCl_mM"] + na
    via_pump = run["sample_inlet"] == "pump_A"
    loading = t < t_load
    pump_prot = np.where(loading, c_sample, 0.0) if via_pump else np.zeros_like(t)
    if via_pump:
        pump_salt = np.where(loading, s_sample, pump_salt)
    Vm = system["mixer_mL"]
    vols = {s["name"]: tube_volume(s["length_cm"], s["id_mm"]) for s in system["tubing"]}
    mix_salt = well_mixed(t, pump_salt, Vm, Q, A)
    mix_prot = well_mixed(t, pump_prot, Vm, Q, 0.0)
    at_valve_salt = delayed(vol, mix_salt, vols["mixer_to_valve"], A)
    at_valve_prot = delayed(vol, mix_prot, vols["mixer_to_valve"], 0.0)
    if via_pump:
        valve_salt, valve_prot = at_valve_salt, at_valve_prot
    else:  # sample pump through the injection valve while loading, system flow otherwise
        valve_salt = np.where(loading, s_sample, at_valve_salt)
        valve_prot = np.where(loading, c_sample, 0.0)
    col_salt = delayed(vol, valve_salt, vols["valve_to_column"], A)
    col_prot = delayed(vol, valve_prot, vols["valve_to_column"], 0.0)
    out_salt, out_prot = simulate_column(col, sma, Qfun, t, col_salt, col_prot, t[-1])
    # after the column: the stated sequence of tubes and cells
    pos_salt, pos_prot, seen = out_salt, out_prot, {}
    for item in system["after_column"]:
        V = vols[item] if item in vols else system["detectors"][item]["cell_volume_mL"] / 2.0
        pos_salt, pos_prot = delayed(vol, pos_salt, V, A), delayed(vol, pos_prot, V, 0.0)
        if item in system["detectors"]:
            seen[item] = (pos_salt.copy(), pos_prot.copy())
            V2 = system["detectors"][item]["cell_volume_mL"] / 2.0
            pos_salt, pos_prot = delayed(vol, pos_salt, V2, A), delayed(vol, pos_prot, V2, 0.0)
    return dict(t=t, Q=Q, vol=vol, pctB=pctB, seen=seen, col_salt=col_salt, out_salt=out_salt, out_prot=out_prot,
                t_load=t_load, t0g=t0g, t_grad=t_grad, na=na)


def uv_reading(salt_mM, prot_mM, protein, uv, wavelength):
    """Absorbance at 280 or 300 nm in mAU, as the detector reports it (no noise, before the auto-zero)."""
    c_mg = prot_mM * protein["molar_mass_g_mol"] / 1000.0
    eps = protein[f"extinction_{wavelength}nm_mL_mg_cm"]
    base = uv_salt_slope_AU_per_M(uv) * (0.35 if wavelength == 300 else 1.0)
    A = eps * c_mg * uv["path_cm"] + base * salt_mM / 1000.0
    s = uv["stray_light_fraction"]
    A_obs = -np.log10((10.0 ** (-A) + s) / (1.0 + s))
    return 1000.0 * A_obs


def uv_salt_slope_AU_per_M(uv):
    return uv["_salt_baseline_AU_per_M_per_cm"] * uv["path_cm"]


def kappa25(salt_mM):
    """NaCl conductivity at 25 C, mS/cm (smooth fit to handbook values; the standards file tabulates it)."""
    c = np.asarray(salt_mM) / 1000.0
    return 118.98 * c - 41.45 * c ** 1.5 + 8.27 * c ** 2


def conductivity_reading(nacl_mM, T_C, cond, kappa_buffer):
    """Conductivity in mS/cm as the monitor reports it: the NaCl part from kappa25 plus the buffer's own
    conductivity at 25 C, the temperature dependence of NaCl solutions, then the monitor's own compensation."""
    k25 = kappa25(nacl_mM) + kappa_buffer
    kT = k25 * (1.0 + ALPHA_NACL * (T_C - 25.0))
    if cond["temperature_compensation_pct_per_K"] is None:
        return kT
    return kT / (1.0 + cond["temperature_compensation_pct_per_K"] / 100.0 * (T_C - 25.0))
