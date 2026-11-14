"""Calibrates one optical-trap session: stiffness per watt, the reference recording's stiffness and detector
sensitivity, and the radius of the bead in that recording.

    python3 trap_reduce.py <session folder> <answer.json>

1. Every recording: int16 counts -> volts; one-sided periodogram of the whole record; bins averaged in blocks of 40
   (the drive line and its neighbours left out) and fitted with the spectrum of an instantaneously sampled
   Ornstein-Uhlenbeck process (aliasing included), times the detector's response.
2. Detector: recordings taken after the log says a silicon QPD was fitted carry its parasitic filter
   alpha^2 + (1 - alpha^2) / (1 + (f / f3)^2); alpha and f3 are fitted per recording, then fixed at their medians.
3. Drive recordings: the line power above the fitted background gives the sensitivity beta, beta^2 times the fitted
   diffusion constant (in V^2/s) gives D, and gamma = k T / D.
4. Segments from the log: a detector change or a new bead starts a new sensitivity; a new bead a new radius; a beam
   realignment a new objective transmission; an x channel recorded without sum normalisation makes the
   sensitivity inversely proportional to the laser power (beta * setpoint constant).
5. Radius: gamma = 6 pi eta(T) R / F(R / h) over the drive recordings of the bead, h = R + xi * stage travel. Where
   a drive recording sits close enough to the coverslip, the focal-shift factor xi is fitted too; otherwise xi = 0.8.
6. Stiffness kappa = 2 pi gamma fc, with gamma measured (drive) or from R and Faxen; per watt, a fit through the origin
   against the objective power of the recordings taken before any realignment.
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

KB = 1.380649e-23
BLOCK = 40


def viscosity(t_c):
    return 2.414e-5 * 10.0 ** (247.8 / (t_c + 273.15 - 140.0))


def faxen(x):
    return 1.0 - 9.0 / 16.0 * x + x ** 3 / 8.0 - 45.0 / 256.0 * x ** 4 - x ** 5 / 16.0


def periodogram(volts, fs):
    n = volts.size
    spec = np.fft.rfft(volts - volts.mean())
    power = 2.0 / (fs * n) * np.abs(spec) ** 2
    freq = np.fft.rfftfreq(n, 1.0 / fs)
    return freq[1:-1], power[1:-1]


def aliased(f, fc, dv, fs):
    dt = 1.0 / fs
    c = math.exp(-2 * math.pi * fc * dt)
    var = dv / (2 * math.pi * fc)
    return 2 * dt * var * (1 - c * c) / (1 + c * c - 2 * c * np.cos(2 * math.pi * f * dt))


def lorentzian(f, fc, dv, fs):
    return dv / (2 * math.pi ** 2 * (fc ** 2 + f ** 2))


def response(f, alpha, f3):
    return alpha ** 2 + (1 - alpha ** 2) / (1 + (f / f3) ** 2)


def binned(f, p, fs, drive_hz, duration):
    keep = (f >= 60.0) & (f <= min(0.92 * fs / 2, 20000.0))
    if drive_hz:
        k = np.abs(f - drive_hz) <= 1.5 / duration
        keep &= ~k
    f, p = f[keep], p[keep]
    n = f.size // BLOCK * BLOCK
    return f[:n].reshape(-1, BLOCK).mean(1), p[:n].reshape(-1, BLOCK).mean(1)


def levenberg(resid, x0, steps=60):
    x = np.array(x0, float)
    lam = 1e-3
    r = resid(x)
    for _ in range(steps):
        jac = np.empty((r.size, x.size))
        for j in range(x.size):
            dx = np.zeros_like(x)
            dx[j] = 1e-6 * max(1.0, abs(x[j]))
            jac[:, j] = (resid(x + dx) - r) / dx[j]
        a = jac.T @ jac
        g = jac.T @ r
        step = np.linalg.solve(a + lam * np.diag(np.diag(a) + 1e-12), -g)
        trial = resid(x + step)
        if np.sum(trial ** 2) < np.sum(r ** 2):
            x, r, lam = x + step, trial, lam * 0.3
            if np.max(np.abs(step)) < 1e-9:
                break
        else:
            lam *= 10.0
    return x


def fit_spectrum(fb, pb, fs, filt=None, free_filter=False, alias=True):
    shape = aliased if alias else lorentzian
    low = pb[: max(3, pb.size // 50)].mean()
    fc0 = 1000.0
    dv0 = low * 2 * math.pi ** 2 * fc0 ** 2

    def model(x):
        fc, dv = math.exp(x[0]), math.exp(x[1])
        out = shape(fb, fc, dv, fs)
        if free_filter:
            out = out * response(fb, min(max(x[2], 0.0), 1.0), math.exp(x[3]))
        elif filt is not None:
            out = out * response(fb, *filt)
        return out

    x = np.array([math.log(fc0), math.log(dv0)] + ([0.4, math.log(8000.0)] if free_filter else []))
    for _ in range(4):                      # weights from the model, not from the noisy data (no 1/n bias)
        weight = model(x)
        x = levenberg(lambda v: math.sqrt(BLOCK) * (pb - model(v)) / weight, x)
    fc, dv = math.exp(x[0]), math.exp(x[1])
    return (fc, dv, (x[2], math.exp(x[3]))) if free_filter else (fc, dv)


def clock(text):
    h, m, *s = (int(v) for v in text.split(":"))
    return h * 60 + m + (s[0] / 60 if s else 0.0)


def log_events(folder):
    out = []
    for line in (folder / "session_log.txt").read_text().splitlines():
        m = re.match(r"(\d\d:\d\d)\s+(.*)", line)
        if not m:
            continue
        t, text = clock(m.group(1)), m.group(2).lower()
        if "silicon" in text and "qpd" in text:
            out.append((t, "silicon"))
        elif "realigned" in text:
            out.append((t, "realign"))
        elif "new bead" in text:
            out.append((t, "bead"))
        elif "without sum normalisation" in text:
            out.append((t, "raw"))
    return out


def reduce(folder, detector_filter=True, alias=True, focal="fit", faxen_on=True, segments=True, ignore=()):
    folder = Path(folder)
    sheet = json.loads((folder / "session_sheet.json").read_text())
    vpc = sheet["volts_per_count"]
    temps = {}
    for line in (folder / "bath_thermometer.csv").read_text().splitlines()[1:]:
        if line.strip():
            t, val = line.split(",")
            temps[t] = float(val)
    meter = np.loadtxt(folder / "power_meter.csv", delimiter=",", skiprows=1, usecols=(1, 2))
    transmission = float(np.mean(meter[:, 1] / meter[:, 0]))
    events = [e for e in log_events(folder) if e[1] not in ignore] if segments else []
    recs = []
    for row in sheet["recordings"]:
        t0 = clock(row["start"])
        n_si = sum(1 for t, kind in events if kind == "silicon" and t < t0)
        n_bead = sum(1 for t, kind in events if kind == "bead" and t < t0)
        n_power = sum(1 for t, kind in events if kind == "realign" and t < t0)
        raw = any(kind == "raw" and t < t0 for t, kind in events)
        volts = np.load(folder / row["file"]).astype(float) * vpc
        f, p = periodogram(volts, row["sample_rate_hz"])
        fb, pb = binned(f, p, row["sample_rate_hz"], row["drive_frequency_hz"], row["duration_s"])
        recs.append(dict(row=row, f=f, p=p, fb=fb, pb=pb, silicon=n_si > 0 and detector_filter, bead=n_bead,
                         det=n_si, power=n_power, temp=temps[row["start"]], raw=raw))
    si = [r for r in recs if r["silicon"]]
    filt = None
    if si:
        fits = [fit_spectrum(r["fb"], r["pb"], r["row"]["sample_rate_hz"], free_filter=True, alias=alias)[2]
                for r in si]
        filt = (float(np.median([a for a, _ in fits])), float(np.median([b for _, b in fits])))
    for r in recs:
        fs = r["row"]["sample_rate_hz"]
        r["fc"], r["dv"] = fit_spectrum(r["fb"], r["pb"], fs, filt if r["silicon"] else None, alias=alias)
        fd, amp = r["row"]["drive_frequency_hz"], r["row"]["drive_amplitude_nm"]
        if fd:
            k = int(np.argmin(np.abs(r["f"] - fd)))
            background = (aliased if alias else lorentzian)(np.array([fd]), r["fc"], r["dv"], fs)[0]
            if r["silicon"]:
                background *= response(fd, *filt)
            w_ex = (r["p"][k] - background) / r["row"]["duration_s"]
            w_th = (amp * 1e-9) ** 2 / (2 * (1 + r["fc"] ** 2 / fd ** 2))
            r["beta"] = math.sqrt(w_th / w_ex)                                      # m / V
            r["gamma_meas"] = KB * (r["temp"] + 273.15) / (r["beta"] ** 2 * r["dv"])
    # sensitivity per (detector, bead) segment, radius and focal factor per bead
    beta_seg, radius, xi_of = {}, {}, {}
    for key in {(r["det"], r["bead"], r["raw"]) for r in recs}:
        group = [r for r in recs if (r["det"], r["bead"], r["raw"]) == key and "beta" in r]
        if key[2]:   # un-normalised x: volts grow with the light, so beta * setpoint is the constant
            beta_seg[key] = float(np.mean([r["beta"] * r["row"]["laser_setpoint_w"] for r in group])) if group else None
        else:
            beta_seg[key] = float(np.mean([r["beta"] for r in group])) if group else None
    for b in {r["bead"] for r in recs}:
        drv = [r for r in recs if r["bead"] == b and "gamma_meas" in r]
        travel = np.array([r["row"]["stage_travel_um"] * 1e-6 for r in drv])
        eta = np.array([viscosity(r["temp"]) for r in drv])
        gam = np.array([r["gamma_meas"] for r in drv])

        def resid(v, fixed_xi=None):
            rad = math.exp(v[0])
            xi = fixed_xi if fixed_xi is not None else v[1]
            h = rad + xi * travel
            corr = faxen(rad / h) if faxen_on else 1.0
            return np.log(6 * math.pi * eta * rad / corr) - np.log(gam)

        far = travel >= 5.0e-6
        near = (~far).sum() >= 2 and far.sum() >= 1
        if focal == "fit" and near and faxen_on:
            # radius from the far recordings, where Faxen is small and the focal factor hardly matters
            v = levenberg(lambda v: resid(v, 0.8)[far], [math.log(5e-7)])
            rad = math.exp(v[0])
            # then the focal factor from the recordings close to the coverslip, at that radius
            w = levenberg(lambda w: resid([math.log(rad), w[0]])[~far], [0.85])
            radius[b], xi_of[b] = rad, float(w[0])
        else:
            xi = 1.0 if focal == "travel" else 0.8
            v = levenberg(lambda v: resid(v, xi), [math.log(5e-7)])
            radius[b], xi_of[b] = math.exp(v[0]), xi
    for r in recs:
        if "gamma_meas" in r:
            r["gamma"] = r["gamma_meas"]
        else:
            rad, xi = radius[r["bead"]], xi_of[r["bead"]]
            h = rad + xi * r["row"]["stage_travel_um"] * 1e-6
            corr = faxen(rad / h) if faxen_on else 1.0
            r["gamma"] = 6 * math.pi * viscosity(r["temp"]) * rad / corr
        r["kappa"] = 2 * math.pi * r["gamma"] * r["fc"]                            # N / m
    first = [r for r in recs if r["power"] == 0]
    p_obj = np.array([transmission * r["row"]["laser_setpoint_w"] for r in first])
    kap = np.array([r["kappa"] for r in first])
    per_watt = float(np.sum(p_obj * kap) / np.sum(p_obj ** 2))
    ref = next(r for r in recs if r["row"]["id"] == sheet["reference_recording"])
    return {
        "stiffness_per_watt_pn_nm_w": per_watt * 1e3,
        "reference_stiffness_pn_nm": ref["kappa"] * 1e3,
        "reference_sensitivity_nm_v": beta_seg[(ref["det"], ref["bead"], ref["raw"])] * 1e9 /
        (ref["row"]["laser_setpoint_w"] if ref["raw"] else 1.0),
        "bead_radius_nm": radius[ref["bead"]] * 1e9,
        "detector_filter": filt,
        "focal_factor": xi_of[ref["bead"]],
    }


if __name__ == "__main__":
    out = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1))
