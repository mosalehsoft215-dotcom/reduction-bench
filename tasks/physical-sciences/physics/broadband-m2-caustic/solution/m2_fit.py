"""Reference analysis: chromatic caustic fit of the broadband laser beam.

1. Spectrometer: background (median of the outer pixels) is removed from every line-out, each pulse is normalised
   to unit area and the pulses are averaged. The spectrometer has no slit: a wavelength lands at its calibrated
   position plus the beam's own spot along the dispersion axis (R11 x + R12 x' from the ETL plane, rms
   sqrt(eps(lambda) (R11^2 beta - 2 R11 R12 alpha + R12^2 gamma))), widened by the camera blur. The record is fitted
   with that forward model: a Gaussian spectrum (centroid, rms relative width) placed with the calibration, every
   wavelength spread by its own spot, integrated over the pixels. The spot depends on the beam the caustic measures,
   so steps 1 and 3 are iterated: fit the caustic, compute the spot from the fitted beam, refit the record, refit.
2. Caustic: every line-out gets a straight background from its outer pixels, is centred on its centroid and
   normalised; the 20 pulses of a setting are averaged.
3. Fit: the 21 averaged profiles are fitted together with beam_model.profile, a mixture over the Gaussian spectrum
   of the run in which every wavelength has its own ETL power, P633 (n(lambda) - 1) / (n(632.8 nm) - 1), and its own
   emittance M2 lambda / 4 pi, weighted by the camera quantum efficiency (camera_qe.csv: the camera counts
   photo-electrons, the spectrometer counts power), blurred by the 6 um camera resolution. Free: M2, beta and alpha
   at the ETL, and an amplitude and offset per setting. Weights from the pulse-to-pulse scatter.
4. The waist: Rayleigh range beta / (1 + alpha^2), distance upstream of the ETL -alpha beta / (1 + alpha^2), and the
   second-moment waist radius 2 sqrt(eps z_R) at the centroid wavelength.
"""
import json
import sys
from functools import partial
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from beam_model import profile, qe_si, waist_from_twiss

RES, PIX = 6.0e-6, 5.5e-6                                                 # notes.txt
R11_SPEC, R12_SPEC, SPEC_RES_MM, SPEC_MM_PER_PX = 0.10, 0.250, 0.08, 0.10
KEYS = ("m2", "waist_radius_um", "waist_distance_mm", "centroid_wavelength_nm", "bandwidth_rms_percent")


def load_csv(path):
    return np.loadtxt(path, delimiter=",", skiprows=1)


def spectrum_data(d):
    """averaged spectrometer record: fraction of the counts in each pixel (background removed, each pulse normalised),
    its standard error, and the pixel centres and edges along the dispersion axis (mm) with the calibration"""
    raw = load_csv(d / "spectrometer_lineouts.csv")[:, 1:]
    cal = load_csv(d / "spectrometer_calibration.csv")[:, 1]
    fr = []
    for y in raw:
        z = y - np.median(np.r_[y[:30], y[-30:]])
        fr.append(z / z.sum())
    fr = np.array(fr)
    sc = (np.arange(cal.size) + 0.5) * SPEC_MM_PER_PX
    return cal, sc, fr.mean(0), fr.std(0, ddof=1) / np.sqrt(len(fr)) + 1e-5


def spot_term(m2, beta, alpha):
    """(rms spot)^2 / lambda on the spectrometer camera, mm^2 per um: eps(lambda) (R11^2 beta - 2 R11 R12 alpha + ...)"""
    eps1 = m2 * 1e-6 / (4 * np.pi)
    return eps1 * (R11_SPEC ** 2 * beta - 2 * R11_SPEC * R12_SPEC * alpha + R12_SPEC ** 2 * (1 + alpha ** 2) / beta) * 1e6


def spectrometer_model(mu_nm, rel, spot, cal, sc, n=41):
    """expected fraction per pixel: a Gaussian spectrum (centroid mu, relative rms rel) placed with the on-axis
    calibration and spread by each wavelength's own spot (sqrt(spot * lambda + camera blur^2)) along the dispersion"""
    from scipy.special import ndtr
    xq, wq = np.polynomial.hermite_e.hermegauss(n)
    wq = wq / wq.sum()
    lam = mu_nm * (1 + xq * rel)
    ok = lam > 0
    pos = np.interp(np.log(lam[ok]), np.log(cal), sc)
    b = np.sqrt(spot * lam[ok] / 1000 + SPEC_RES_MM ** 2)
    lo, hi = sc - 0.5 * SPEC_MM_PER_PX, sc + 0.5 * SPEC_MM_PER_PX
    return sum(w * (ndtr((hi - c) / bb) - ndtr((lo - c) / bb)) for c, bb, w in zip(pos, b, wq[ok]))


def fit_spectrum(data, spot):
    cal, sc, f, e = data
    m0 = np.sum(f * cal) / np.sum(f)
    r0 = np.sqrt(np.sum(f * (cal - m0) ** 2) / np.sum(f)) / m0
    r = least_squares(lambda q: (q[2] * spectrometer_model(q[0], q[1], spot, cal, sc) + q[3] - f) / e,
                      [m0, r0, 1.0, 0.0], x_scale=[1.0, 0.001, 0.01, 1e-5])
    return float(r.x[0]), float(r.x[1]), m0, r0


def scan_profiles(d):
    a = load_csv(d / "caustic_lineouts.csv")
    sett, cur, y = a[:, 0].astype(int), a[:, 1], a[:, 3:]
    n = y.shape[1]
    px = np.arange(n)
    xc = (px - n / 2 + 0.5) * PIX
    outer = np.r_[0:120, n - 120:n]
    profs, errs, currents = [], [], []
    for s in np.unique(sett):
        rows = []
        for yy in y[sett == s]:
            c = np.polyfit(px[outer], yy[outer], 1)
            z = yy - np.polyval(c, px)
            core = np.abs(px - np.argmax(np.convolve(z, np.ones(15) / 15, "same"))) < 800
            x0 = np.sum(z[core] * xc[core]) / np.sum(z[core])
            zs = np.interp(xc + x0, xc, z, left=0, right=0)              # centred
            rows.append(zs / (np.sum(zs) * PIX))
        rows = np.array(rows)
        profs.append(rows.mean(0))
        errs.append(rows.std(0, ddof=1) / np.sqrt(len(rows)) + 2e-3 * rows.mean(0).max())
        currents.append(cur[sett == s][0])
    return xc, np.array(currents), np.array(profs), np.array(errs)


def fit(xc, currents, profs, errs, lam_um, spread, chromatic=True, eps_scales=True, calibrated=True, res=RES,
        starts=((1.5, 8.0, -150.0),), qe=True):
    nset = len(currents)

    def resid(q):
        m2, beta, alpha = q[0], q[1], q[2]
        out = []
        for i in range(nset):
            m = q[3 + 2 * i] * profile(xc, currents[i], lam_um, m2, beta, alpha, spread, res, n=21 if spread > 0 else 1,
                                       chromatic=chromatic, eps_scales=eps_scales, calibrated=calibrated,
                                       qe=qe if callable(qe) else (qe_si if qe else None)) + q[4 + 2 * i]
            out.append((profs[i] - m) / errs[i])
        return np.concatenate(out)

    best = None
    for m0, b0, a0 in starts:
        q0 = [m0, b0, a0] + [1.0, 0.0] * nset
        lo = [0.3, 0.05, -5000.0] + [0.5, -50.0] * nset
        hi = [20.0, 500.0, 5000.0] + [2.0, 50.0] * nset
        f = least_squares(resid, q0, bounds=(lo, hi), x_scale="jac", max_nfev=300)
        if best is None or f.cost < best.cost:
            best = f
    return best.x, 2 * best.cost / (best.fun.size - best.x.size)


def solve(d, deconvolve=True, spectrometer="forward", mono=False, **kw):
    """full analysis; spectrometer "forward" (the record fitted with its own spot) or "moments" (the plain mean and rms
    of the record, spot ignored) for the alternative analyses"""
    data = spectrum_data(d)
    xc, cur, profs, errs = scan_profiles(d)
    if kw.get("qe", True) is True:
        kw["qe"] = partial(qe_si, table=load_csv(d / "camera_qe.csv"))
    starts = ((1.5, 8.0, -150.0), (1.2, 15.0, -250.0))
    lam_nm, spread, m0, r0 = fit_spectrum(data, 0.0)
    if spectrometer == "moments" or not deconvolve:
        lam_nm, spread = m0, r0
    for _ in range(4 if (deconvolve and spectrometer == "forward") else 1):
        q, chi2 = fit(xc, cur, profs, errs, lam_nm / 1000, 0.0 if mono else spread, starts=starts, **kw)
        if not deconvolve or spectrometer != "forward":
            break
        lam_nm, spread, _, _ = fit_spectrum(data, spot_term(q[0], q[1], q[2]))
        starts = ((q[0], q[1], q[2]),)
    if deconvolve and spectrometer == "forward":
        q, chi2 = fit(xc, cur, profs, errs, lam_nm / 1000, 0.0 if mono else spread, starts=starts, **kw)
    return dict(q=q, chi2=chi2, lam_nm=lam_nm, measured=r0, measured_centroid=m0, spread=spread)


def results(r):
    m2, beta, alpha = (float(v) for v in r["q"][:3])
    zr, dw = waist_from_twiss(beta, alpha)
    eps = m2 * r["lam_nm"] * 1e-9 / (4 * np.pi)
    return dict(zip(KEYS, (round(m2, 5), round(2e6 * float(np.sqrt(eps * zr)), 5), round(1e3 * float(dw), 5),
                           round(float(r["lam_nm"]), 4), round(100 * float(r["spread"]), 5))))


def main(data_dir, out_path):
    r = solve(Path(data_dir))
    out = results(r)
    out.update(record_rms_bandwidth_percent=round(100 * float(r["measured"]), 5),
               record_mean_wavelength_nm=round(float(r["measured_centroid"]), 4), reduced_chi2=round(float(r["chi2"]), 4))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
