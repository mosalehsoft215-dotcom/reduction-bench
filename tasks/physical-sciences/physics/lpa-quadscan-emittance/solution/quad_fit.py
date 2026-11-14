"""Reference analysis: chromatic quadrupole-scan fit of the laser-plasma electron beam.

1. Spectrometer: background (median of the outer pixels) is removed from every line-out, counts per
   pixel are turned into counts per MeV with the pixel widths of the calibration, each shot is
   normalised to unit area and the shots are averaged. This gives the mean kinetic energy and the
   measured rms relative spread of the whole run. The spectrometer has no focusing, so the measured
   spread contains the beam's own horizontal spot on its screen (R11 x + R12 x' from the quadrupole
   entrance, plus the 0.10 mm screen resolution), converted to energy with the local dispersion of the
   calibration. That blur depends on the beam the scan measures, so the fit of step 3 and the
   deconvolution are iterated: fit, compute the blur from the fitted beam, remove it in quadrature
   from the measured spread, refit.
2. Scan: every line-out gets a straight background from its outer pixels, is centred on its centroid
   and normalised; the 20 shots of a setting are averaged.
3. Fit: the 21 averaged profiles are fitted together with quad_model.profile, a mixture over the
   Gaussian momentum distribution of the run (each momentum focused by its own k = k0 / (1 + delta)),
   blurred by the 12 um screen resolution. Free: geometric emittance, beta and alpha at the
   quadrupole entrance, and an amplitude and offset per setting. Weights from the shot-to-shot scatter.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from quad_model import C_GEV, ME, gamma_rel, momentum_gev, profile

G_PER_A, RES, PIX = 0.4900, 12.0e-6, 5.0e-6                       # notes.txt
R11_SPEC, R12_SPEC, SPEC_RES_MM, SPEC_MM_PER_PX = 0.92, 1.85, 0.10, 0.10
KEYS = ("emittance_n_mm_mrad", "beta_m", "alpha", "energy_MeV", "energy_spread_rms_percent")


def load_csv(path):
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    return a


def spectrum(d):
    raw = load_csv(d / "spectrometer_lineouts.csv")[:, 1:]
    cal = load_csv(d / "spectrometer_calibration.csv")[:, 1]
    edges = np.r_[cal[0] - 0.5 * (cal[1] - cal[0]), 0.5 * (cal[1:] + cal[:-1]), cal[-1] + 0.5 * (cal[-1] - cal[-2])]
    width = np.diff(edges)
    dens = []
    for y in raw:
        b = np.median(np.r_[y[:40], y[-40:]])
        f = (y - b) / width
        dens.append(f / np.sum(f * width))
    f = np.mean(dens, axis=0)
    w = f * width
    mean = np.sum(w * cal) / np.sum(w)
    rel = np.sqrt(np.sum(w * (cal - mean) ** 2) / np.sum(w)) / mean     # measured, blur included
    e_tot = mean + ME
    kin = (e_tot ** 2 - ME ** 2) / (e_tot * mean)                    # dT/T = kin * dp/p
    j = np.interp(mean, cal, np.arange(cal.size))
    disp = np.interp(j, np.arange(cal.size - 1) + 0.5, np.diff(np.log(cal))) / SPEC_MM_PER_PX   # dlnT/ds, per mm
    return mean, rel, kin, disp


def blur_rel(eps, beta, alpha, disp):
    """relative kinetic-energy blur of the spectrometer: the beam's horizontal spot on its screen"""
    s2 = eps * (R11_SPEC ** 2 * beta - 2 * R11_SPEC * R12_SPEC * alpha + R12_SPEC ** 2 * (1 + alpha ** 2) / beta) * 1e6
    return np.sqrt(s2 + SPEC_RES_MM ** 2) * disp


def scan_profiles(d):
    a = load_csv(d / "quadscan_lineouts.csv")
    sett, cur, y = a[:, 0].astype(int), a[:, 1], a[:, 3:]
    n = y.shape[1]
    px = np.arange(n)
    xc = (px - n / 2 + 0.5) * PIX
    outer = np.r_[0:150, n - 150:n]
    profs, errs, currents = [], [], []
    for s in np.unique(sett):
        rows = []
        for yy in y[sett == s]:
            c = np.polyfit(px[outer], yy[outer], 1)
            z = yy - np.polyval(c, px)
            core = np.abs(px - np.argmax(np.convolve(z, np.ones(9) / 9, "same"))) < 300
            x0 = np.sum(z[core] * xc[core]) / np.sum(z[core])
            zs = np.interp(xc + x0, xc, z, left=0, right=0)              # centred
            rows.append(zs / (np.sum(zs) * PIX))
        rows = np.array(rows)
        profs.append(rows.mean(0))
        errs.append(rows.std(0, ddof=1) / np.sqrt(len(rows)) + 2e-3 * rows.mean(0).max())
        currents.append(cur[sett == s][0])
    return xc, np.array(currents), np.array(profs), np.array(errs)


def fit(xc, currents, profs, errs, p_gev, spread_p, chromatic=True, res=RES, starts=((2.0e-9, 2.0, -40.0),)):
    k0 = C_GEV * G_PER_A * currents / p_gev
    nset = len(currents)
    sp = spread_p if chromatic else 0.0

    def resid(q):
        eps, beta, alpha = q[0] * 1e-9, q[1], q[2]
        out = []
        for i in range(nset):
            m = q[3 + 2 * i] * profile(xc, k0[i], eps, beta, alpha, sp, res, n=21 if sp > 0 else 1) + q[4 + 2 * i]
            out.append((profs[i] - m) / errs[i])
        return np.concatenate(out)

    best = None
    for e0, b0, a0 in starts:
        q0 = [e0 * 1e9, b0, a0] + [1.0, 0.0] * nset
        lo = [1e-3, 0.01, -500.0] + [0.5, -50.0] * nset
        hi = [1e3, 100.0, 500.0] + [2.0, 50.0] * nset
        f = least_squares(resid, q0, bounds=(lo, hi), x_scale="jac", max_nfev=300)
        if best is None or f.cost < best.cost:
            best = f
    return best.x, 2 * best.cost / (best.fun.size - best.x.size)


def solve(d):
    """full analysis; returns the fitted beam and the run's spectrum figures"""
    t_mean, measured_t, kin, disp = spectrum(d)
    p = momentum_gev(t_mean)
    xc, cur, profs, errs = scan_profiles(d)
    spread_t = measured_t
    starts = ((2.0e-9, 2.0, -40.0), (4.0e-9, 1.0, -20.0))
    for _ in range(4):                          # fit <-> spectrometer-blur deconvolution, to self-consistency
        q, chi2 = fit(xc, cur, profs, errs, p, spread_t / kin, starts=starts)
        blur = blur_rel(q[0] * 1e-9, q[1], q[2], disp)
        spread_t = np.sqrt(max(measured_t ** 2 - blur ** 2, 1e-8))
        starts = ((q[0] * 1e-9, q[1], q[2]),)
    q, chi2 = fit(xc, cur, profs, errs, p, spread_t / kin, starts=starts)
    return dict(q=q, chi2=chi2, t_mean=t_mean, measured_t=measured_t, spread_t=spread_t, kin=kin, blur=blur,
                p=p, xc=xc, cur=cur, profs=profs, errs=errs)


def main(data_dir, out_path):
    r = solve(Path(data_dir))
    q, chi2, t_mean, spread_t, measured_t, blur = r["q"], r["chi2"], r["t_mean"], r["spread_t"], r["measured_t"], r["blur"]
    epsn = q[0] * 1e-9 * gamma_rel(t_mean) * np.sqrt(1 - 1 / gamma_rel(t_mean) ** 2)
    out = dict(zip(KEYS, (round(epsn * 1e6, 5), round(float(q[1]), 5), round(float(q[2]), 4), round(float(t_mean), 3),
                          round(100 * float(spread_t), 4))))
    out.update(measured_spectrometer_spread_percent=round(100 * float(measured_t), 4),
               spectrometer_blur_percent=round(100 * float(blur), 4), reduced_chi2=round(float(chi2), 4))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
