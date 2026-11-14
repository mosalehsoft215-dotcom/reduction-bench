"""Second, independent reduction of a vault session (same call as the reference).

- Weight lifts (coil terminated by the digitizer input Z alone), peeled one at a time: each lift is fitted on its own
  window (from 1 s before its logged time to 0.5 s before the next event) for its onset, w0, damping, amplitude and a
  straight baseline; a window whose rms misfit exceeds 1 % of the lift's peak (an earthquake) is not peeled and not
  used, the others are subtracted from the rest of the record before the next lift is fitted.  Of the clean windows,
  those with an amplitude under half their median (a lift the solenoid missed) are dropped; w0, the lifts' damping
  h_o and the amplitude are the medians of the rest, and amplitude / (m_w g) is the lift scale S_d G Z / (R_c + Z).
- Calibration-coil steps (coil terminated by R_p = R_s Z / (R_s + Z)): the twelve signed windows (1.5 s before to
  10 s after each logged time), each detrended by a straight line through its first 1.5 s and last 2 s; windows whose
  peak-to-peak size is under half (a missed command) or over twice (an earthquake) the median are dropped, the rest
  are combined by their sample-wise median and fitted with the step response for its onset, damping and amplitude,
  w0 held at the lifts' value.
- q = stacked amplitude / (I x lift scale) = K R_p (R_c + Z) / ((R_c + R_p) Z).  With G = K l_sig / l_cal (one gap)
  and h - h_o = G^2 (Z - R_p) / (2 m w0 (R_c + R_p)(R_c + Z)), R_c is the smallest positive root of the cubic
  2 m w0 (h - h_o) R_p^2 (R_c + Z)^3 = (l_sig / l_cal)^2 q^2 Z^2 (Z - R_p)(R_c + R_p) (numpy.roots); then K, G and
  sensitivity = lift scale (R_c + Z) / Z x R_p / (R_c + R_p) x |H(5 Hz)|.  The listed S_d, the sheet K and the
  vault thermometer are not used.
"""
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402


def vel(t, m, w0, h):
    wd = w0 * np.sqrt(1 - h * h)
    x = np.clip(t, 0, None)
    return np.where(t > 0, np.exp(-h * w0 * x) * np.sin(wd * x) / (m * wd), 0.0)


def run(folder):
    folder = Path(folder)
    sen = json.loads((folder / "sensor.json").read_text())
    ses = json.loads((folder / "session.json").read_text())
    ev = np.genfromtxt(folder / "events.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    rate, delay, sd = ses["sample_rate_hz"], ses["digitizer_fir_delay_s"], ses["digitizer_counts_per_volt"]
    m, rs = sen["moving_mass_kg"], sen["shunt_resistance_ohm"]

    def events(rec, up):
        sel = ev["record"] == rec
        return ev["time_s"][sel].astype(float), np.where(ev["event"][sel] == up, 1.0, -1.0)

    # weight lifts, peeled
    y = np.genfromtxt(folder / "open_circuit_record.csv", skip_header=1).astype(float)
    t = np.arange(y.size) / rate
    times, signs = events("open_circuit", "weight lifted")
    spec = np.abs(np.fft.rfft((y - y.mean()) * np.hanning(y.size)))
    freq = np.fft.rfftfreq(y.size, 1 / rate)
    band = (freq > 0.4) & (freq < 12.0)
    w_guess = 2 * np.pi * freq[band][np.argmax(spec[band])]
    rows = []
    for k, (t_ev, s) in enumerate(zip(times, signs)):
        end = times[k + 1] - 0.5 if k + 1 < len(times) else t[-1]
        sel = (t >= t_ev - 1.0) & (t <= end)
        x, yy = t[sel] - t_ev, y[sel]

        def model(p, x=x):
            d, w0, hm, a, c0, c1 = p
            return a * s * vel(x - d, m, w0, hm) + c0 + c1 * x
        scale = np.ptp(yy) * m * w_guess / 2
        best = None
        for hm0 in (0.02, 0.08):
            p0 = [delay, w_guess, hm0, scale, float(np.median(yy[x < 0])), 0.0]
            sol = least_squares(lambda p: model(p) - yy, p0, x_scale="jac",
                                bounds=([0.0, 0.5 * w_guess, 0.001, 0.0, -np.inf, -np.inf],
                                        [0.5, 1.5 * w_guess, 0.5, np.inf, np.inf, np.inf]))
            if best is None or sol.cost < best.cost:
                best = sol
        d, w0, hm, a, _, _ = best.x
        misfit = float(np.sqrt(np.mean(best.fun ** 2)))
        clean = misfit <= 0.01 * a / (m * w0)                 # a quake leaves a misfit of percents of the lift
        rows.append((w0, hm, a, clean))
        if clean:                                             # peel this lift off the rest of the record
            y = y - a * s * vel(t - t_ev - d, m, w0, hm)
    rows = np.array(rows)
    # a window a quake disturbed (misfit over 1 % of the lift's peak) or a lift the solenoid missed (amplitude under
    # half the median of the clean windows) is not used
    good = rows[:, 3] > 0
    good &= rows[:, 2] >= 0.5 * np.median(rows[good, 2])
    w0, hm, amp = (float(np.median(c)) for c in rows[good, :3].T)
    scale = amp / (ses["weight_lift_mass_g"] * 1e-3 * ses["local_gravity_m_per_s2"])

    # calibration-coil steps, stacked
    yc = np.genfromtxt(folder / "calibration_record.csv", skip_header=1).astype(float)
    tc = np.arange(yc.size) / rate
    pieces = []
    for t_ev, s in zip(*events("calibration", "current on")):
        sel = (tc >= t_ev - 1.5) & (tc < t_ev + 10.0)
        x, yy = tc[sel] - t_ev, yc[sel]
        ref = (x < 0) | (x > 8.0)
        pieces.append(s * (yy - np.polyval(np.polyfit(x[ref], yy[ref], 1), x)))
    n = min(p.size for p in pieces)
    stack = np.array([p[:n] for p in pieces])
    # a command the relay missed leaves a flat window and a quake a wild one: windows whose size (peak-to-peak) is
    # under half or over twice the median are left out, and the rest are combined by their sample-wise median
    size = np.ptp(stack, axis=1)
    keep = (size > 0.5 * np.median(size)) & (size < 2.0 * np.median(size))
    xs, ys = x[:n], np.median(stack[keep], axis=0)
    best = None
    for h0 in (0.4, 0.7, 0.9):
        sol = least_squares(lambda p: p[2] * vel(xs - p[0], m, w0, p[1]) - ys,
                            [delay, h0, np.ptp(ys) * m * w0], x_scale="jac",
                            bounds=([0.0, 0.05, -np.inf], [0.5, 0.98, np.inf]))
        if best is None or sol.cost < best.cost:
            best = sol
    h = best.x[1]

    # the coil is terminated by Z alone during the lifts and by R_p = R_s || Z during the steps:
    # scale = S_d G Z / (R_c + Z), q = step amplitude / (I scale) = K R_p (R_c + Z) / ((R_c + R_p) Z), G = rho K and
    # h - h_o = G^2 (Z - R_p) / (2 m w0 (R_c + R_p)(R_c + Z)) give the cubic
    # D R_p^2 (R_c + Z)^3 = rho^2 q^2 Z^2 (Z - R_p) (R_c + R_p), D = 2 m w0 (h - h_o); its smallest positive root is R_c
    z = ses["digitizer_input_impedance_ohm"]
    rp = rs * z / (rs + z)
    q = best.x[2] / (ses["calibrator_current_ma"] * 1e-3 * scale)
    rho = sen["signal_coil_wire_in_gap_m"] / sen["calibration_coil_wire_in_gap_m"]
    dd = 2 * m * w0 * (h - hm)
    # in x = R_c / Z: dd R_p^2 Z^2 (x + 1)^3 - rho^2 q^2 Z^2 (Z - R_p) (x + R_p / Z) = 0
    c3 = dd * rp * rp * z * z
    c1 = rho * rho * q * q * z * z * (z - rp)
    roots = np.roots([c3, 3 * c3, 3 * c3 - c1, c3 - c1 * rp / z])
    rc = z * float(min(r.real for r in roots if abs(r.imag) <= 1e-9 * max(abs(r), 1e-12) and r.real > 0))
    k = q * z * (rc + rp) / (rp * (rc + z))
    r = 5.0 * 2 * np.pi / w0
    resp = r * r / np.sqrt((1 - r * r) ** 2 + (2 * h * r) ** 2)
    sens = scale * (rc + z) / z * rp / (rc + rp) * resp
    return {"generator_constant_v_s_per_m": float(k * rho), "sensitivity_5hz_counts_per_m_per_s": float(sens)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(run(sys.argv[1])))
