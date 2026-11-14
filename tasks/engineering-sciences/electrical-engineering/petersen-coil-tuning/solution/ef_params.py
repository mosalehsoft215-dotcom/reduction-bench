"""Earth-fault parameters of a compensated 20 kV network from one event folder.

    python3 ef_params.py <event folder> <results.json>

The folder holds ef_event.cfg / ef_event.dat (the disturbance record), event_log.txt and station.json.

Method (the same for every event):
(a) Record. Scaling and ratios from the cfg; every channel on its own time axis, (input number - 1) x 25 us late,
    because one converter samples eight inputs in turn at 5 kHz (the cfg skew fields are empty); samples at the end of
    the converter range are left out.
(b) Frequency fitted on U_L1L2; phasors at it with the 5th and 7th harmonics.
(c) Coil. During the fault and before the auxiliary resistor closes, U_NE / I_ASC is the coil inductance alone (the
    measuring winding sits on the core); the decaying offset of the fault inception is fitted with the phasors. If the
    resistor closes before the feeder trips, the step of I_ASC / U_NE at the contactor gives its referred value.
    Between the network neutral and the inductance: main winding at the oil temperature and one third of the earthing
    transformer's zero-sequence impedance.
(d) Ring-down after the faulted feeder opens: the network without that feeder rings with the coil branch (with the
    resistor if it is in), and with the fixed coils of the other feeders; capacitance and leakage fitted in the time
    domain to U_NE and I_ASC.
(e) The faulted feeder comes back from the primary test: its relay reading corrected for the core-balance CT, scaled
    from the test voltage and frequency, its own fixed coil of the test day removed, and every change since the test
    that touches it applied (cables moved in or out, cables added).
(f) Answers at 20 kV, 50 Hz, resistor out, all feeders and all of today's fixed coils in: capacitive current, coil
    current, detuning, active current, and the natural displacement from the pre-fault U_NE through the series
    impedance, as two independent phase-capacitance differences put on a balanced supply.
"""
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

U_N = 20000.0 / np.sqrt(3.0)
W_N = 2 * np.pi * 50.0
A120 = np.exp(2j * np.pi / 3)
HARM = (1, 5, 7)


def cable_c(km, uf_per_km):
    return 3 * uf_per_km * 1e-6 * km


def read_record(folder, inputs, slot_us, drop_clipped=True):
    cfg = (folder / "ef_event.cfg").read_text().splitlines()
    na, nd = (int(re.sub(r"\D", "", x)) for x in cfg[1].split(",")[1:3])
    analog = [cfg[2 + i].split(",") for i in range(na)]
    digital = [cfg[2 + na + i].split(",") for i in range(nd)]
    rate, n = cfg[4 + na + nd].split(",")
    rate, n = float(rate), int(n)
    stamp = lambda s: datetime.strptime(s.strip(), "%d/%m/%Y,%H:%M:%S.%f")
    first, trig = stamp(cfg[5 + na + nd]), stamp(cfg[6 + na + nd])
    rows = np.loadtxt(folder / "ef_event.dat", delimiter=",", dtype=np.int64)
    t0 = np.arange(n) / rate - (trig - first).total_seconds()
    out, counts = {}, {}
    for i, a in enumerate(analog):
        gain = float(a[5]) * (float(a[10]) / float(a[11]) if a[12].strip().upper() == "S" else 1.0)
        value = gain * rows[:, 2 + i] + float(a[6])
        if drop_clipped:
            value = np.where(np.abs(rows[:, 2 + i]) >= int(a[9]) - 2, np.nan, value)
        out[a[1]] = (t0 + (inputs[a[1]] - 1) * slot_us * 1e-6, value)
        counts[a[1]] = gain
    dig = {d[1]: (t0, rows[:, 2 + na + i]) for i, d in enumerate(digital)}
    return out, dig, counts, trig


def basis(t, w):
    cols = [np.ones_like(t)]
    for h in HARM:
        cols += [np.sqrt(2.0) * np.cos(h * w * t), -np.sqrt(2.0) * np.sin(h * w * t)]
    return np.column_stack(cols)


def phasor(sig, lo, hi, w):
    t, y = sig
    m = (t >= lo) & (t <= hi) & np.isfinite(y)
    c = np.linalg.lstsq(basis(t[m], w), y[m], rcond=None)[0]
    return c[1] + 1j * c[2]


def phasor_decaying(sig, lo, hi, w, tau):
    t, y = sig
    m = (t >= lo) & (t < hi) & np.isfinite(y)
    t, y = t[m], y[m]
    cols = [np.ones_like(t), np.exp(-(t - lo) / tau)]
    for h in HARM:
        cols += [np.cos(h * w * t), -np.sin(h * w * t)]
    c = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)[0]
    return (c[2] + 1j * c[3]) / np.sqrt(2.0)


def phasor_offset(sig, lo, hi, w):
    """phasor with an offset that decays with a time constant fitted to the data (the inception transient decays with
    the loop of the coil and the fault, whose resistance is not known)"""
    from scipy.optimize import minimize_scalar as ms

    def sse(log_tau):
        t, y = sig
        m = (t >= lo) & (t < hi) & np.isfinite(y)
        tt, yy = t[m], y[m]
        cols = [np.ones_like(tt), np.exp(-(tt - lo) / np.exp(log_tau))]
        for h in HARM:
            cols += [np.cos(h * w * tt), -np.sin(h * w * tt)]
        a = np.column_stack(cols)
        return np.sum((yy - a @ np.linalg.lstsq(a, yy, rcond=None)[0]) ** 2)

    grid = np.log(np.geomspace(0.005, 2.0, 40))
    g0 = grid[np.argmin([sse(g) for g in grid])]
    best = ms(sse, bounds=(g0 - 0.2, g0 + 0.2), method="bounded").x
    return phasor_decaying(sig, lo, hi, w, np.exp(best)), np.exp(best)


def system_frequency(sig):
    t, y = sig
    m = np.isfinite(y)
    t, y = t[m], y[m]

    def misfit(f):
        b = basis(t, 2 * np.pi * f)
        return np.sum((y - b @ np.linalg.lstsq(b, y, rcond=None)[0]) ** 2)

    grid = np.arange(49.0, 51.0001, 0.02)
    f0 = grid[np.argmin([misfit(f) for f in grid])]
    return minimize_scalar(misfit, bounds=(f0 - 0.02, f0 + 0.02), method="bounded", options=dict(xatol=1e-9)).x


def ring_columns(c, g, branch, w, t_u, t_i, t0):
    rs, ls, lp, rr, lf = branch
    yf = 0.0 if lf is None else 1.0 / lf
    if rr is None:
        a = np.array([[-g / c, -1.0 / c, -1.0 / c], [1.0 / (ls + lp), -rs / (ls + lp), 0.0], [yf, 0.0, 0.0]])
        out_u = np.array([lp / (ls + lp), -lp * rs / (ls + lp), 0.0])
        out_i = np.array([0.0, 1.0, 0.0])
    else:
        a = np.array([[-g / c, -1.0 / c, 0.0, -1.0 / c], [1.0 / ls, -(rs + rr) / ls, rr / ls, 0.0],
                      [0.0, rr / lp, -rr / lp, 0.0], [yf, 0.0, 0.0, 0.0]])
        out_u = np.array([0.0, rr, -rr, 0.0])
        out_i = np.array([0.0, 1.0, 0.0, 0.0])
    lam, vec = np.linalg.eig(a)
    cols_u, cols_i = [], []
    for k in range(len(lam)):
        if lam[k].real < -2000.0 or lam[k].imag < -1e-9 or abs(lam[k]) < 1e-9:
            continue
        z_u = (out_u @ vec[:, k]) * np.exp(lam[k] * (t_u - t0))
        z_i = (out_i @ vec[:, k]) * np.exp(lam[k] * (t_i - t0))
        cols_u += [z_u.real] + ([z_u.imag] if lam[k].imag > 1e-9 else [])
        cols_i += [z_i.real] + ([z_i.imag] if lam[k].imag > 1e-9 else [])
    src = np.zeros(len(a))
    src[0] = 1.0 / c
    for h in HARM:
        x = np.linalg.solve(1j * h * w * np.eye(len(a)) - a, src)
        z_u, z_i = (out_u @ x) * np.exp(1j * h * w * t_u), (out_i @ x) * np.exp(1j * h * w * t_i)
        cols_u += [z_u.real, z_u.imag]
        cols_i += [z_i.real, z_i.imag]
    return np.column_stack(cols_u), np.column_stack(cols_i)


def ring_fit(rec, counts, branch, w, lo, hi):
    (t_u, y_u), (t_i, y_i) = rec["U_NE"], rec["I_ASC"]
    m_u = (t_u >= lo) & (t_u <= hi) & np.isfinite(y_u)
    m_i = (t_i >= lo) & (t_i <= hi) & np.isfinite(y_i)
    t_u, y_u, t_i, y_i = t_u[m_u], y_u[m_u] / counts["U_NE"], t_i[m_i], y_i[m_i] / counts["I_ASC"]

    def residual(q):
        b_u, b_i = ring_columns(q[0] * 1e-6, q[1] * 1e-3, branch, w, t_u, t_i, lo)
        b = np.vstack([b_u / counts["U_NE"], b_i / counts["I_ASC"]])
        y = np.concatenate([y_u, y_i])
        return y - b @ np.linalg.lstsq(b, y, rcond=None)[0]

    best = None
    for c0 in np.arange(10.0, 90.0, 2.5):
        for g0 in (0.1, 0.5, 1.5, 3.0):
            cost = np.sum(residual((c0, g0)) ** 2)
            if best is None or cost < best[0]:
                best = (cost, c0, g0)
    fit = least_squares(residual, best[1:], x_scale=(1.0, 0.1), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    # the forced 50 Hz part of the fit is the natural unbalance injection of the remaining network
    b_u, b_i = ring_columns(fit.x[0] * 1e-6, fit.x[1] * 1e-3, branch, w, t_u, t_i, lo)
    b = np.vstack([b_u / counts["U_NE"], b_i / counts["I_ASC"]])
    coef = np.linalg.lstsq(b, np.concatenate([y_u, y_i]), rcond=None)[0]
    k = b.shape[1] - 2 * len(HARM)
    inject = (coef[k] - 1j * coef[k + 1]) / np.sqrt(2.0)          # rms phasor of J, c du/dt = ... + J
    return fit.x[0] * 1e-6, fit.x[1] * 1e-3, inject


def event_times(folder, tripped, trig):
    times = {}
    for line in (folder / "event_log.txt").read_text().splitlines():
        m = re.match(r"(\d\d:\d\d:\d\d\.\d{3})\s+(\S+)\s+(.+?)\s{2,}(\S.*)$", line)
        if not m:
            continue
        hms, bay, sig, state = m.groups()
        t = datetime.strptime(f"{trig:%Y-%m-%d} {hms}", "%Y-%m-%d %H:%M:%S.%f")
        times[(bay, sig.strip(), state.strip())] = (t - trig).total_seconds()
    return times[(tripped, "Q0", "OPEN")]


def solve(folder, resistor="measured", skew="inputs", clipped="drop", add_back=True, changes=True, cbct=True,
          unbalance_test=True, winding_check=True, board_lag=True):
    folder = Path(folder)
    st = json.loads((folder / "station.json").read_text())
    coil, et = st["arc_suppression_coil"], st["earthing_transformer"]["z0_ohm_per_phase_50Hz"]
    tripped = st["event"]["faulted_feeder"]
    inputs = st["recorder"]["inputs_used"]
    if skew == "adjacent":
        inputs = {k: i + 1 for i, k in enumerate(("U_L1L2", "U_L2L3", "U_NE", "I_ASC"))}
    slot = 1e6 / (st["recorder"]["sampling_Hz"] * 8) if skew else 0.0
    rec, dig, counts, trig = read_record(folder, inputs, slot, clipped == "drop")
    t_trip = event_times(folder, tripped, trig)
    f = system_frequency(rec["U_L1L2"])
    w = 2 * np.pi * f
    t_end = rec["U_NE"][0][-1]

    t_d, res_on = dig["ASC1 RES ON"]
    t_res = float(t_d[np.argmax(res_on > 0)]) if res_on.any() else np.inf
    end1 = min(t_res, t_trip) - 0.008
    coil_r = coil["main_winding_ohm_at_20C"] * (234.5 + coil["oil_temperature_C_on_event_day"]) / (234.5 + 20.0)
    rs, ls = coil_r + et[0] / 3.0, et[1] / W_N / 3.0
    if board_lag:
        # The measuring winding senses the induced voltage of a loss-free core, so before the resistor closes U_NE
        # leads I_ASC by exactly 90 degrees. Any other angle is a lag of the current input board: shift its samples.
        ratio = phasor_offset(rec["I_ASC"], 0.03, end1, w)[0] / phasor_offset(rec["U_NE"], 0.03, end1, w)[0]
        lag = -(np.angle(ratio) + np.pi / 2) / w
        t_i, x_i = rec["I_ASC"]
        rec["I_ASC"] = (t_i - lag, x_i)
    if winding_check and coil["measuring_winding"].get("remark"):
        # the winding was rewound, so its ratio is no longer the nameplate's. During a bolted fault the network
        # neutral sits at minus the faulted phase voltage, which the line VTs give (V connection: V1 = (2 U12 + U23)
        # / 3, ...); the neutral is the winding's EMF plus the drop of the coil current across the coil copper and
        # the earthing transformer, so the winding's scale s solves s U_NE + Z I_ASC = -V_k
        uf = phasor_offset(rec["U_NE"], 0.03, end1, w)[0]
        jf = phasor_offset(rec["I_ASC"], 0.03, end1, w)[0]
        u12, u23 = phasor(rec["U_L1L2"], 0.03, end1, w), phasor(rec["U_L2L3"], 0.03, end1, w)
        v1 = (2 * u12 + u23) / 3.0
        phases = [v1, v1 - u12, v1 - u12 - u23]
        drop = (rs + 1j * w * ls) * jf
        k = int(np.argmin([abs(np.angle(-(uf + drop) / vp)) for vp in phases]))
        scale = (-phases[k] - drop) / uf
        if abs(np.angle(scale)) < np.radians(1.0) and abs(abs(scale) - 1.0) < 0.05:
            t_u, x_u = rec["U_NE"]
            rec["U_NE"] = (t_u, x_u * abs(scale))
    i1, tau_i = phasor_offset(rec["I_ASC"], 0.03, end1, w)
    u1, tau_u = phasor_offset(rec["U_NE"], 0.03, end1, w)
    lp = -1.0 / (w * (i1 / u1).imag)
    y1 = i1 / u1
    rr = None
    if t_res < t_trip - 0.05 and resistor != "ignored":
        if resistor == "measured":
            u2, i2 = phasor(rec["U_NE"], t_res + 0.01, t_trip - 0.008, w), phasor(rec["I_ASC"], t_res + 0.01, t_trip - 0.008, w)
            rr = 1.0 / (i2 / u2 - y1).real
        else:
            ratio = coil["auxiliary_winding"]["ratio_V"]
            rr = coil["auxiliary_winding"]["resistor_ohm_nameplate"] * (ratio[0] / ratio[1]) ** 2

    fixed = {c["feeder"]: c["current_A_at_20kV_50Hz_with_its_earthing_transformer"] for c in st["fixed_coils_today"]}
    i_fix_ring = sum(a for fd, a in fixed.items() if fd != tripped)
    lf_ring = U_N / (W_N * i_fix_ring) if i_fix_ring else None
    c, g, inject_rest = ring_fit(rec, counts, (rs, ls, lp, rr, lf_ring), w, t_trip + 0.006, t_end)
    c_rest, g_rest = c, g
    # pre-fault line voltages and displacement, used for the unbalance of the tripped feeder and for U0
    u0, i0 = phasor(rec["U_NE"], -1.0, -0.012, w), phasor(rec["I_ASC"], -1.0, -0.012, w)
    v0 = u0 + (rs + 1j * w * ls) * i0
    uab, ubc = phasor(rec["U_L1L2"], -1.0, -0.012, w), phasor(rec["U_L2L3"], -1.0, -0.012, w)
    m_line = np.array([[(uab + ubc).real, ubc.real], [(uab + ubc).imag, ubc.imag]])

    test = st["primary_test"]
    i_fix_all = sum(fixed.values())
    y_fix = (lambda wx: 1.0 / (1j * wx * U_N / (W_N * i_fix_all))) if i_fix_all else (lambda wx: 0.0)
    if add_back:
        row = next(r for r in test["readings"] if r["feeder"] == tripped)
        mag, ang = row["residual_A"], row["angle_deg"]
        if cbct:
            cb = st["feeder_relays"]["core_balance_ct"]
            mag, ang = mag / (1 + cb["ratio_error_percent"] / 100.0), ang - cb["phase_displacement_deg_secondary_leading"]
        w_t = 2 * np.pi * test["frequency_Hz"]
        u_t = test["busbar_kV"] * 1000.0 / np.sqrt(3.0)
        reading0 = mag * np.exp(1j * np.radians(ang))
        # the fixed coil that sat on this feeder on the test day read with it
        fixed_test = dict(fixed)
        for ch in st["changes_since_test"] if changes else []:
            if ch["type"] == "fixed_coil_added":
                fixed_test.pop(ch["feeder"], None)
            elif ch["type"] == "fixed_coil_removed":
                fixed_test[ch["feeder"]] = ch["current_A_at_20kV_50Hz"]

        def feeder_today(reading):
            """capacitance and leakage of the tripped feeder today from its test-day admittance"""
            y_t = reading / u_t
            if tripped in fixed_test:
                y_t += 1j * fixed_test[tripped] / U_N * W_N / w_t
            c_add, g_add = y_t.imag / w_t, y_t.real
            for ch in st["changes_since_test"] if changes else []:
                if ch["type"] == "cable_moved" and ch["src"] == tripped:
                    c_add -= cable_c(ch["km"], ch["uf"])
                elif ch["type"] == "cable_moved" and ch["dst"] == tripped:
                    c_add += cable_c(ch["km"], ch["uf"])
                elif ch["type"] == "cable_added" and ch["feeder"] == tripped:
                    c_add += cable_c(ch["km"], ch["uf"])
            return c_add, g_add

        c_add, g_add = feeder_today(reading0)
        if unbalance_test:
            # the relay also read the feeder's own natural unbalance current j w_t sum_p C_p E_p (solid L1 fault, U0 =
            # -E1 as reference; cable changes are balanced, so the deviations are today's). Its two phase-capacitance
            # differences come from the record: the unbalance injection of the whole network before the fault minus
            # that of the network left after the trip (the forced part of the ring-down fit)
            l12, l23 = -u_t * (1 - A120 ** 2), -u_t * (A120 ** 2 - A120)
            for _ in range(2):
                y_full = (g_rest + g_add) + 1j * w * (c_rest + c_add) + 1.0 / (rs + 1j * w * (ls + lp)) + y_fix(w)
                s_feeder = -y_full * v0 + inject_rest
                dd1, dd2 = np.linalg.solve(m_line, [(s_feeder / (1j * w)).real, (s_feeder / (1j * w)).imag])
                c_add, g_add = feeder_today(reading0 - 1j * w_t * (dd1 * (l12 + l23) + dd2 * l23))
        c, g = c_rest + c_add, g_rest + g_add

    def admittance(wx):
        branch = 1.0 / (rs + 1j * wx * (ls + lp))
        return g + 1j * wx * c + branch + y_fix(wx), branch

    y50, branch = admittance(W_N)
    i_c, i_l = W_N * c * U_N, -branch.imag * U_N
    source = -admittance(w)[0] * v0
    d1, d2 = np.linalg.solve(m_line, [(source / (1j * w)).real, (source / (1j * w)).imag])
    e12, e23 = 20000.0 * np.exp(1j * np.pi / 6), 20000.0 * np.exp(-1j * np.pi / 2)
    unbalance = abs(1j * W_N * (d1 * (e12 + e23) + d2 * e23))
    out = dict(capacitive_current_A=i_c, coil_current_A=i_l, detuning_percent=100.0 * (i_c - i_l - i_fix_all) / i_c,
               active_current_A=y50.real * U_N, neutral_displacement_V=unbalance / abs(y50))
    return {k: float(v) for k, v in out.items()}, dict(frequency_Hz=f, coil_H=lp, resistor_ohm=rr, ring_C_uF=c * 1e6)


if __name__ == "__main__":
    res, info = solve(sys.argv[1])
    Path(sys.argv[2]).parent.mkdir(parents=True, exist_ok=True)
    Path(sys.argv[2]).write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res), json.dumps(info, default=str))
