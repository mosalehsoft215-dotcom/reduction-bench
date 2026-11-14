"""Independent second reduction of a spin-echo measurement: linear demodulation of the echo groups.

    python3 demod_solver.py <measurement folder> <answer file>

- Counts are dead-time corrected, the closed-beam rate of their group is taken off, and they are divided by the
  incident beam, row by row, with NumPy arrays. With the source steady (logbook) the counting time measures the
  beam; a run during which the logbook records a source change keeps its own monitor.
- Each echo group is fitted linearly with three terms: a constant, cos and sin of 2 pi I / P, both carrying the
  wavelength-spread envelope around the current echo-point estimate. For a scatterer and Fourier time, the echo
  point is the weighted circular mean of the phases of all its detector groups. It is iterated four times, with the
  envelope re-centred on it each time, and the amplitude is the length of the (cos, sin) pair after the last pass.
- Solvent subtraction and the standard's normalisation follow the stated physics. The solvent path comes from a
  logbook clause that names the solvent's cell, when there is one.
- The Zimm fit starts from a weighted line through ln(-ln S) against ln tau and is refined with
  scipy.optimize.curve_fit. The incoherent fraction uses the count-weighted mean flipping ratio at the reference group.
"""
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

K_B = 1.380649e-23


def main(folder):
    sh = json.loads((folder / "measurement_sheet.json").read_text())
    full = np.genfromtxt(folder / "counts.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    closed = full[full["scatterer"] == "closed"]
    table = full[full["scatterer"] != "closed"]
    lam = sh["wavelength_A"]
    # room background per group from the closed-beam counts
    bg_rate = np.zeros(len(sh["detector_groups"]))
    for g in range(len(bg_rate)):
        sel = closed["group"] == g
        if sel.any():
            bg_rate[g] = closed["counts"][sel].sum() / closed["time_s"][sel].sum()
    # beam normalisation: with the source steady the counting time is the beam measure; a run during which the
    # logbook records a source change keeps its monitor, scaled to the others through their monitor rate
    log = (folder / "logbook.txt").read_text().lower()
    troubled = {w for w in ("sample", "solvent", "standard") if f"during the {w} run" in log}
    steady = np.isin(table["scatterer"], [w for w in ("sample", "solvent", "standard") if w not in troubled])
    r0 = table["monitor"][steady].sum() / table["time_s"][steady].sum()
    beam = np.where(steady, table["time_s"] * r0, table["monitor"]).astype(float)
    per = sh["phase_coil_period_A_at_8A"] * 8.0 / lam
    spread = sh["wavelength_spread_sigma_rel"]
    raw = table["counts"].astype(float)
    corr = raw / (1 - raw * sh["detector_dead_time_s"] / table["time_s"]) - bg_rate[table["group"]] * table["time_s"]
    rate = corr / beam
    var = np.maximum(raw, 1) / beam ** 2
    who, grp, tix, kind, ph = table["scatterer"], table["group"], table["tau_index"], table["kind"], table["phase_A"]
    n_g, n_t = len(sh["detector_groups"]), len(sh["main_coil_currents_A"])

    def pick(w, g, k, kd):
        return (who == w) & (grp == g) & (tix == k) & (kind == kd)

    amplitude = {}
    for w in ("sample", "solvent", "standard"):
        for k in range(n_t):
            centre = 0.0
            for _ in range(4):
                sums, results = 0j, []
                for g in range(n_g):
                    sel = pick(w, g, k, "echo")
                    I, y, v = ph[sel], rate[sel], var[sel]
                    arg = 2 * math.pi * I / per
                    env = np.exp(-0.5 * (2 * math.pi * (I - centre) / per * spread) ** 2)
                    X = np.column_stack([np.ones_like(I), np.cos(arg) * env, np.sin(arg) * env])
                    wt = 1 / v
                    A = (X * wt[:, None]).T @ X
                    coef = np.linalg.solve(A, (X * wt[:, None]).T @ y)
                    cov = np.linalg.inv(A)
                    results.append((coef, cov))
                    z = complex(coef[1], coef[2])
                    sums += z * abs(z) / max(cov[1, 1] + cov[2, 2], 1e-30)
                centre = math.atan2(sums.imag, sums.real) / (2 * math.pi) * per
            for g, (coef, cov) in enumerate(results):
                proj = coef[1] * math.cos(2 * math.pi * centre / per) + coef[2] * math.sin(2 * math.pi * centre / per)
                amplitude[(w, g, k)] = (proj, cov[1, 1] + cov[2, 2])

    def single(w, g, k, kd):
        sel = pick(w, g, k, kd)
        return float(rate[sel][0])

    path = 1.0
    for clause in re.split(r"[;.\n]\s*", (folder / "logbook.txt").read_text()):
        found = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*mm", clause)
        if found and "solvent" in clause:                 # the clause that says where the solvent went
            path = float(found.group(1)) / sh["sample_cell_path_mm"]
    scale = (1 - sh["polymer_volume_fraction"]) * sh["transmission"]["sample"] / sh["transmission"]["solvent"] / path
    q = np.array([4 * math.pi / lam * math.sin(math.radians(d["two_theta_deg"]) / 2) for d in sh["detector_groups"]])
    tau = 0.1863 * np.array(sh["main_coil_currents_A"]) * sh["main_coil_T_m_per_A"] * lam ** 3
    qq, tt, ss, ww, flips = [], [], [], [], {g: [] for g in range(n_g)}
    for g in range(n_g):
        for k in range(n_t):
            up = {w: single(w, g, k, "up") for w in ("sample", "solvent", "standard")}
            dn = {w: single(w, g, k, "down") for w in ("sample", "solvent", "standard")}
            a_s, v_s = amplitude[("sample", g, k)]
            a_b, v_b = amplitude[("solvent", g, k)]
            a_r, v_r = amplitude[("standard", g, k)]
            chains = 0.5 * ((up["sample"] - dn["sample"]) - scale * (up["solvent"] - dn["solvent"]))
            S = (a_s - scale * a_b) / chains / (a_r / (0.5 * (up["standard"] - dn["standard"])))
            qq.append(q[g]); tt.append(tau[k]); ss.append(S)
            ww.append(1 / ((v_s + scale ** 2 * v_b) / chains ** 2 * (0.5 * (up["standard"] - dn["standard"]) / a_r) ** 2 + 1e-12))
            total = (up["sample"] + dn["sample"]) - scale * (up["solvent"] + dn["solvent"])
            p0 = (up["standard"] - dn["standard"]) / (up["standard"] + dn["standard"])
            flips[g].append((2 * chains / total / p0, total))
    qq, tt, ss, ww = map(np.array, (qq, tt, ss, ww))
    T = sh["temperature_K"]
    xdata = np.vstack([qq, tt])

    def zimm(x, eta, beta):
        gam = K_B * T * (x[0] * 1e10) ** 3 / (6 * math.pi * eta * 1e-3) * 1e-9
        return np.exp(-(gam * x[1]) ** beta)

    ok = (ss > 0.05) & (ss < 0.95)
    gfac = K_B * T * (qq * 1e10) ** 3 / (6 * math.pi * 1e-3) * 1e-9
    slope, icept = np.polyfit(np.log(tt[ok] * gfac[ok]), np.log(-np.log(ss[ok])), 1)
    eta0 = math.exp(-icept / slope)
    popt, _ = curve_fit(zimm, xdata, ss, p0=[eta0, slope], sigma=1 / np.sqrt(ww), absolute_sigma=True,
                        bounds=([0.05, 0.2], [50.0, 1.5]))
    g_ref = int(np.argmin(abs(q - sh["reference_q_inv_A"])))
    ratio = sum(r * w for r, w in flips[g_ref]) / sum(w for _, w in flips[g_ref])
    return {"zimm_viscosity_mpa_s": float(popt[0]), "stretch_exponent": float(popt[1]),
            "incoherent_fraction_at_reference": 0.75 * (1 - ratio)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(main(Path(sys.argv[1])), indent=1) + "\n")
