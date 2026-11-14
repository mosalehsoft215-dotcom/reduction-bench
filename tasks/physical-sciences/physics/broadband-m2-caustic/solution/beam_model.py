"""Second-moment optics of the caustic line, shared by the generator and the reference (copy in solution/).

The scanned element is an electrically tunable lens (ETL), thin, followed by a drift DRIFT to the camera.
Its optical power is calibrated at 632.8 nm, P633(I) = P_OFFSET + P_SLOPE * I; at another wavelength the power
scales with the lens fluid's refractivity, P(lambda) = P633 (n(lambda) - 1) / (n(0.6328) - 1), Cauchy fluid
n(lambda) = N_A + N_B / lambda^2 (lambda in um).

Each wavelength of the pulse propagates as its own beam of second-moment emittance eps = M2 lambda / (4 pi), with
the same Twiss parameters (beta, alpha) at the ETL for every wavelength (the beam comes from a stable focus whose
Rayleigh range does not depend on the wavelength). At the camera its rms size is
    sigma^2 = eps (M11^2 beta - 2 M11 M12 alpha + M12^2 (1 + alpha^2) / beta),   M11 = 1 - DRIFT P, M12 = DRIFT.
"""
import numpy as np

DRIFT = 1.500                       # m, ETL to camera
P_OFFSET, P_SLOPE = 5.00, 0.0400    # dioptres, dioptres per mA, at 632.8 nm
N_A, N_B = 1.5300, 0.01100          # lens fluid: n = N_A + N_B / lambda_um^2
LAM_CAL = 0.6328                    # um


def n_fluid(lam_um):
    return N_A + N_B / np.asarray(lam_um, float) ** 2


def power(current_ma, lam_um, calibrated=True):
    """optical power (dioptres) of the ETL at the wavelength; calibrated=False takes the 632.8 nm value as is"""
    p633 = P_OFFSET + P_SLOPE * np.asarray(current_ma, float)
    if not calibrated:
        return p633 * np.ones_like(np.asarray(lam_um, float))
    return p633 * (n_fluid(lam_um) - 1) / (n_fluid(LAM_CAL) - 1)


def size2(p, eps, beta, alpha, drift=DRIFT):
    m11, m12 = 1 - drift * p, drift
    return eps * (m11 ** 2 * beta - 2 * m11 * m12 * alpha + m12 ** 2 * (1 + alpha ** 2) / beta)


def hermite(n=31):
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


QE_TABLE = ((650, 0.62), (700, 0.55), (750, 0.46), (800, 0.37), (850, 0.28), (900, 0.20), (950, 0.13),
            (1000, 0.075), (1050, 0.04))        # caustic camera quantum efficiency (camera_qe.csv)


def qe_si(lam_um, table=QE_TABLE):
    t = np.array(table, float)
    return np.interp(np.asarray(lam_um, float) * 1e3, t[:, 0], t[:, 1], left=t[0, 1], right=t[-1, 1])


def profile(xc, current_ma, lam_um, m2, beta, alpha, spread, res, n=31, chromatic=True, eps_scales=True,
            calibrated=True, qe=qe_si):
    """expected line-out shape (unit area, per metre) at camera positions xc: a mixture over the Gaussian spectrum
    (relative rms `spread` of the wavelength about lam_um) of Gaussians, each wavelength with its own lens power and
    emittance, blurred by the camera resolution. chromatic=False: one power for all wavelengths; eps_scales=False:
    the same emittance for all wavelengths; qe: the camera's quantum efficiency weighting the wavelengths it records
    (None: every wavelength counts by its power)."""
    xq, wq = hermite(n) if spread > 0 else (np.zeros(1), np.ones(1))
    if qe is not None:
        wq = wq * qe(lam_um * (1 + xq * spread))
        wq = wq / wq.sum()
    out = np.zeros_like(np.asarray(xc, float))
    for d, w in zip(xq * spread, wq):
        lam = lam_um * (1 + d)
        p = power(current_ma, lam if chromatic else lam_um, calibrated)
        eps = m2 * (lam if eps_scales else lam_um) * 1e-6 / (4 * np.pi)
        s2 = size2(p, eps, beta, alpha) + res ** 2
        out += w * np.exp(-xc ** 2 / (2 * s2)) / np.sqrt(2 * np.pi * s2)
    return out


def twiss_from_waist(zr, d):
    """(beta, alpha) at the ETL of a beam whose waist (Rayleigh range zr) lies a distance d upstream"""
    return zr + d * d / zr, -d / zr


def waist_from_twiss(beta, alpha):
    """(Rayleigh range, distance of the waist upstream of the ETL) from (beta, alpha) at the ETL"""
    zr = beta / (1 + alpha * alpha)
    return zr, -alpha * zr
