"""Twin solver: every event analysed again by separately written code (numpy only; nothing imported from solution/
or authoring/provenance/).

    python twin_solver.py <events folder> <truth.json> <output json>

Where it differs from solution/ef_params.py:
  - its own line-by-line COMTRADE reader; each input's slot in the converter cycle from station.json; clipped samples
    removed before anything is fitted;
  - frequency from the phase advance of U_L1L2 between the two halves of the record, phasors by Hann-weighted
    demodulation fits rather than least squares with harmonics;
  - the coil inductance from late fault phasors only (no inception-offset model);
  - ring-down without a circuit simulation: one damped oscillation at complex frequency s plus the forced terms and a
    slow exponential; the ratio of the two complex amplitudes at s is the coil branch as the measuring winding sees
    it, which gives the auxiliary resistor (when it is in) from the ring-down itself rather than from the contactor
    step; the ringing network then follows in closed form from the admittances at the neutral summing to zero at s;
  - unbalance through the positive- and negative-sequence parts of the supply rather than phase differences.
"""
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

U_NOM = 20000.0
U_PH = U_NOM / math.sqrt(3.0)
F_NOM = 50.0


def load(folder, st):
    lines = (Path(folder) / "ef_event.cfg").read_text().splitlines()
    n_an, n_dig = int(lines[1].split(",")[1].rstrip("Aa")), int(lines[1].split(",")[2].rstrip("Dd"))
    rate, n = lines[4 + n_an + n_dig].split(",")
    when = lambda s: datetime.strptime(s, "%d/%m/%Y,%H:%M:%S.%f")
    first, trig = when(lines[5 + n_an + n_dig]), when(lines[6 + n_an + n_dig])
    raw = np.array([[int(v) for v in line.split(",")] for line in (Path(folder) / "ef_event.dat").read_text().split()])
    base = (raw[:, 0] - 1) / float(rate) - (trig - first).total_seconds()
    slot = 1.0 / (float(rate) * 8)
    chans, lsb = {}, {}
    for k in range(n_an):
        f = lines[2 + k].split(",")
        gain = float(f[5]) * (float(f[10]) / float(f[11]) if f[12].strip() in ("S", "s") else 1.0)
        keep = np.abs(raw[:, 2 + k]) < 32760
        lag = (st["recorder"]["inputs_used"][f[1]] - 1) * slot
        chans[f[1]] = ((base + lag)[keep], (gain * raw[:, 2 + k] + float(f[6]))[keep])
        lsb[f[1]] = gain
    res_on = raw[:, 2 + n_an + 1]
    t_res = float(base[np.argmax(res_on > 0)]) if res_on.any() else math.inf
    return chans, lsb, trig, t_res


def opened(folder, trig, feeder):
    for line in (Path(folder) / "event_log.txt").read_text().splitlines():
        if re.match(rf"\d\d:\d\d:\d\d\.\d+\s+{feeder}\s+Q0\s+OPEN", line):
            hms = line.split()[0]
            return (datetime.strptime(f"{trig:%Y-%m-%d} {hms}", "%Y-%m-%d %H:%M:%S.%f") - trig).total_seconds()
    raise ValueError("no breaker opening in the event log")


def demod(sig, f, lo, hi):
    """rms phasor at frequency f between lo and hi, phase referred to t = 0"""
    t, y = sig
    m = (t >= lo) & (t <= hi)
    # Hann-weighted least squares (cosine, sine, offset): unlike a plain weighted sum it stays unbiased when samples
    # are missing, as they are where the coil-current input clipped
    win = np.sqrt(0.5 - 0.5 * np.cos(2 * np.pi * (t[m] - lo) / (hi - lo)))
    arg = 2 * np.pi * f * t[m]
    cols = np.column_stack([np.cos(arg), np.sin(arg), np.ones_like(arg)])
    c = np.linalg.lstsq(cols * win[:, None], y[m] * win, rcond=None)[0]
    return (c[0] - 1j * c[1]) / math.sqrt(2.0)


def frequency(sig):
    t = sig[0]
    mid = 0.5 * (t[0] + t[-1])
    f = F_NOM
    for _ in range(6):
        a, b = demod(sig, f, t[0], mid), demod(sig, f, mid, t[-1])
        f += np.angle(b / a) / (2 * np.pi * (0.5 * (t[-1] - t[0])))
    return f


def oscillation(chans, lsb, f, lo):
    """complex frequency of the free oscillation and its complex amplitudes in U_NE and I_ASC"""
    w = 2 * np.pi * f
    data = []
    for name in ("U_NE", "I_ASC"):
        t, y = chans[name]
        m = t >= lo
        data.append((t[m], y[m] / lsb[name]))

    def columns(q, t):
        tau = t - lo
        cols = [np.exp(q[0] * tau) * np.cos(q[1] * tau), -np.exp(q[0] * tau) * np.sin(q[1] * tau), np.exp(-q[2] * tau)]
        for h in (1, 5, 7):
            cols += [np.cos(h * w * t), np.sin(h * w * t)]
        return np.column_stack(cols)

    def residual(q):
        out, amp = [], []
        for t, y in data:
            b = columns(q, t)
            c = np.linalg.lstsq(b, y, rcond=None)[0]
            out.append(y - b @ c)
            amp.append(complex(c[0], c[1]))
        return np.concatenate(out), amp

    best = None
    for fr in np.arange(40.0, 70.01, 0.5):
        for sig in (-8.0, -16.0, -30.0):
            r = residual(np.array([sig, 2 * np.pi * fr, 0.5]))[0]
            if best is None or r @ r < best[0]:
                best = (r @ r, np.array([sig, 2 * np.pi * fr, 0.5]))
    q, damp = best[1], 1e-3
    r = residual(q)[0]
    for _ in range(60):
        jac = np.empty((len(r), 3))
        for k in range(3):
            step = np.zeros(3)
            step[k] = 1e-6 * max(abs(q[k]), 1.0)
            jac[:, k] = (residual(q + step)[0] - residual(q - step)[0]) / (2 * step[k])
        h = jac.T @ jac
        dq = np.linalg.solve(h + damp * np.diag(np.diag(h)), -jac.T @ r)
        r_new = residual(q + dq)[0]
        if r_new @ r_new < r @ r:
            done = np.all(np.abs(dq) < 1e-10 * np.maximum(np.abs(q), 1.0))
            q, r, damp = q + dq, r_new, damp / 5.0
            if done:
                break
        else:
            damp *= 10.0
            if damp > 1e12:
                break
    r, amp = residual(q)
    n_u = len(data[0][0])
    rms = (float(np.sqrt(np.mean(r[:n_u] ** 2))), float(np.sqrt(np.mean(r[n_u:] ** 2))))
    forced = []
    for (tt, y), name in zip(data, ("U_NE", "I_ASC")):
        c = np.linalg.lstsq(columns(q, tt), y, rcond=None)[0]
        forced.append((c[3] - 1j * c[4]) / math.sqrt(2.0) * lsb[name])      # rms phasor of the 50 Hz part
    return complex(q[0], q[1]), amp[0] * lsb["U_NE"], amp[1] * lsb["I_ASC"], q[2], rms, forced


def analyse(folder):
    st = json.loads((Path(folder) / "station.json").read_text())
    chans, lsb, trig, t_res = load(folder, st)
    tripped = st["event"]["faulted_feeder"]
    t_open = opened(folder, trig, tripped)
    coil = st["arc_suppression_coil"]
    z_et = complex(*st["earthing_transformer"]["z0_ohm_per_phase_50Hz"])
    r_main = coil["main_winding_ohm_at_20C"] * (234.5 + coil["oil_temperature_C_on_event_day"]) / (234.5 + 20.0)
    fixed = {c["feeder"]: c["current_A_at_20kV_50Hz_with_its_earthing_transformer"] for c in st["fixed_coils_today"]}
    i_fix_all = sum(fixed.values())
    i_fix_ring = sum(a for fd, a in fixed.items() if fd != tripped)
    f = frequency(chans["U_L1L2"])
    w = 2 * np.pi * f
    w50 = 2 * np.pi * F_NOM
    # Current input board: the winding EMF of a loss-free core leads the coil current by a quarter period, so the
    # mean phase gap from 0.25 s after inception to the resistor or the trip, less 90 degrees, is the board's lag.
    hi0 = min(t_res, t_open) - 0.008
    gap = np.angle(demod(chans["U_NE"], f, 0.25, hi0) / demod(chans["I_ASC"], f, 0.25, hi0))
    board = (gap - math.pi / 2) / w
    tt, yy = chans["I_ASC"]
    chans["I_ASC"] = (tt - board, yy)
    if coil["measuring_winding"].get("remark"):
        # a rewound winding's ratio is checked against the bolted-fault condition: neutral = minus the faulted phase
        # voltage, where the neutral is the winding EMF plus the coil current's drop over copper and transformer
        hi = min(t_res, t_open) - 0.008
        uf, jf = demod(chans["U_NE"], f, 0.25, hi), demod(chans["I_ASC"], f, 0.25, hi)
        a, b = demod(chans["U_L1L2"], f, 0.25, hi), demod(chans["U_L2L3"], f, 0.25, hi)
        va = (2 * a + b) / 3.0
        drop = (r_main + z_et.real / 3.0 + 1j * z_et.imag / 3.0 * w / w50) * jf
        cand = [(-(v + drop) / uf) for v in (va, va - a, va - a - b)]
        best = min(cand, key=lambda z: abs(np.angle(z)))
        if abs(np.angle(best)) < math.radians(1.0) and abs(abs(best) - 1.0) < 0.05:
            tt, yy = chans["U_NE"]
            chans["U_NE"] = (tt, yy * abs(best))
            lsb["U_NE"] *= abs(best)

    s, a_u, a_i, slow, rms, (u_post, i_post) = oscillation(chans, lsb, f, t_open + 0.006)
    end = min(t_res, t_open) - 0.008
    u1, i1 = demod(chans["U_NE"], f, max(0.35, end - 0.6), end), demod(chans["I_ASC"], f, max(0.35, end - 0.6), end)
    l_fault = abs(u1) / (w * abs(i1))
    resistor_in = t_res < t_open
    y_coil = a_i / a_u
    r_res = 1.0 / (y_coil - 1.0 / (s * l_fault)).real if resistor_in else math.inf

    r_series = r_main + z_et.real / 3.0
    l_series = z_et.imag / w50 / 3.0
    z_par = s * l_fault if not resistor_in else s * l_fault * r_res / (s * l_fault + r_res)
    z_s = r_series + s * l_series + z_par
    rest = -1.0 / z_s - (i_fix_ring / (U_PH / w50)) / s if i_fix_ring else -1.0 / z_s
    cap = rest.imag / s.imag
    leak = rest.real - s.real * cap

    # the network left after the trip at the system frequency, and its forced 50 Hz state after the trip
    z_par_w = 1j * w * l_fault if not resistor_in else 1j * w * l_fault * r_res / (1j * w * l_fault + r_res)
    y_rest_w = leak + 1j * w * cap + 1.0 / (r_series + 1j * w * l_series + z_par_w)
    if i_fix_ring:
        y_rest_w += i_fix_ring / U_PH * w50 / (1j * w)
    v0_post = u_post + (r_series + 1j * w * l_series) * i_post

    test = st["primary_test"]
    row = next(r for r in test["readings"] if r["feeder"] == tripped)
    cb = st["feeder_relays"]["core_balance_ct"]
    w_t = 2 * np.pi * test["frequency_Hz"]
    y_t = row["residual_A"] / (1 + cb["ratio_error_percent"] / 100.0) * np.exp(
        1j * math.radians(row["angle_deg"] - cb["phase_displacement_deg_secondary_leading"])) / (test["busbar_kV"] * 1000.0 / math.sqrt(3.0))
    on_test = dict(fixed)
    for ch in st["changes_since_test"]:
        if ch["type"] == "fixed_coil_added":
            on_test.pop(ch["feeder"], None)
        elif ch["type"] == "fixed_coil_removed":
            on_test[ch["feeder"]] = ch["current_A_at_20kV_50Hz"]
    if tripped in on_test:
        y_t += 1j * on_test[tripped] / U_PH * w50 / w_t
    def tripped_now(y_read):
        c7, g7 = y_read.imag / w_t, y_read.real
        for ch in st["changes_since_test"]:
            c_cable = 3 * ch.get("uf", 0.0) * 1e-6 * ch.get("km", 0.0)
            if ch["type"] == "cable_moved":
                c7 += c_cable * ((ch["dst"] == tripped) - (ch["src"] == tripped))
            elif ch["type"] == "cable_added" and ch["feeder"] == tripped:
                c7 += c_cable
        return c7, g7

    c7, g7 = tripped_now(y_t)
    cap_all, leak_all = cap + c7, leak + g7

    def total(wx):
        coil_b = 1.0 / (r_series + 1j * wx * (l_series + l_fault))
        y_fx = (i_fix_all / U_PH * w50 / (1j * wx)) if i_fix_all else 0.0
        return leak_all + 1j * wx * cap_all + coil_b + y_fx, coil_b

    y50, coil50 = total(w50)
    i_c = w50 * cap_all * U_PH
    i_l = -coil50.imag * U_PH
    t_first = chans["U_NE"][0][0]
    u0, i0 = demod(chans["U_NE"], f, t_first, -0.012), demod(chans["I_ASC"], f, t_first, -0.012)
    v0 = u0 + (r_series + 1j * w * l_series) * i0
    uab, ubc = demod(chans["U_L1L2"], f, t_first, -0.012), demod(chans["U_L2L3"], f, t_first, -0.012)
    a = np.exp(2j * np.pi / 3)
    uca = -uab - ubc
    e_pos = (uab + a * ubc + a * a * uca) / 3.0 / (math.sqrt(3.0) * np.exp(1j * np.pi / 6))
    e_neg = (uab + a * a * ubc + a * uca) / 3.0 / (math.sqrt(3.0) * np.exp(-1j * np.pi / 6))
    mat = np.array([[(e_pos + e_neg).real, (-e_pos + e_neg).imag], [(e_pos + e_neg).imag, (e_pos - e_neg).real]])

    def unbalance(drive):
        """K = C1 + a^2 C2 + a C3 from drive = sum_p C_p E_p = K e_pos + conj(K) e_neg"""
        d_re, d_im = np.linalg.solve(mat, [drive.real, drive.imag])
        return complex(d_re, d_im)

    # the tripped feeder's relay read its own unbalance current too: j w_t K_F E1 with E1 = -U at the solid L1 test
    # fault; K_F = K of the whole network before the fault minus K of the network left after the trip
    k_rest = unbalance(-y_rest_w * v0_post / (1j * w))
    for _ in range(2):
        k_f = unbalance(-total(w)[0] * v0 / (1j * w)) - k_rest
        c7, g7 = tripped_now(y_t + 1j * w_t * k_f)
        cap_all, leak_all = cap + c7, leak + g7
    y50, coil50 = total(w50)
    i_c = w50 * cap_all * U_PH
    i_l = -coil50.imag * U_PH
    k_all = unbalance(-total(w)[0] * v0 / (1j * w))
    d_re, d_im = k_all.real, k_all.imag
    u_ne0 = abs(w50 * U_PH * complex(d_re, d_im) / y50)
    res = dict(capacitive_current_A=i_c, coil_current_A=i_l, detuning_percent=100.0 * (i_c - i_l - i_fix_all) / i_c,
               active_current_A=y50.real * U_PH, neutral_displacement_V=u_ne0)
    info = dict(frequency_hz=f, oscillation_hz=s.imag / (2 * np.pi), decay_per_s=-s.real, coil_H=l_fault,
                resistor_ohm=r_res if resistor_in else None, fit_rms_counts=rms)
    return {k: float(v) for k, v in res.items()}, info


if __name__ == "__main__":
    BOUND = dict(capacitive_current_A=("rel", 0.1), coil_current_A=("rel", 0.05), detuning_percent=("abs", 0.04),
                 active_current_A=("rel", 0.8), neutral_displacement_V=("rel", 0.6))   # the verifier's limits
    root, truth = Path(sys.argv[1]), json.loads(Path(sys.argv[2]).read_text())
    report = {}
    for ev in sorted(truth):
        result, details = analyse(root / ev)
        table = {}
        for key, (kind, tol) in BOUND.items():
            err = (result[key] / truth[ev][key] - 1.0) * 100.0 if kind == "rel" else result[key] - truth[ev][key]
            table[key] = dict(value=result[key], error=err, error_over_bound=abs(err) / tol)
        report[ev] = dict(result=table, worst=max(v["error_over_bound"] for v in table.values()), details=details)
        print(ev, f"worst {report[ev]['worst']:.3f}", json.dumps(details, default=float), flush=True)
    Path(sys.argv[3]).write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="ascii", newline="\n")
