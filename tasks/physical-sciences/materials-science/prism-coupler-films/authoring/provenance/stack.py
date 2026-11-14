"""Plane-wave reflectance of the prism / air gap / film / substrate stack at the prism base (exact Airy sums).

N is the effective index n_p sin(angle in the prism at the base).  Each layer j has kz_j = k sqrt(n_j^2 - N^2)
(Im >= 0); TE uses q_j = kz_j, TM q_j = kz_j / n_j^2.  A small extinction kappa of the film stands for the guided
light that leaves the pressed spot of the film (and its scattering), which is what turns a mode into a dip.
"""
import numpy as np

LAM_UM = 0.6328
PRISM_N = 1.9648


def kz(n, N, k):
    z = np.sqrt((n * n - N * N).astype(complex) if np.ndim(N) else complex(n * n - N * N))
    return k * np.where(np.imag(z) < 0, -z, z)


def reflectance(N, film_n, film_kappa, film_um, gap_um, sub_n, pol, film_ne=None, film_ny=None):
    """film_n is the in-plane (ordinary) index; film_ne the index along the film normal (optic axis), seen only by
    TM light: kz = k (n_o / n_e) sqrt(n_e^2 - N^2), admittance kz / n_o^2."""
    k = 2 * np.pi / LAM_UM
    nf = film_n + 1j * film_kappa
    ne = (film_n if film_ne is None else film_ne) + 1j * film_kappa
    layers = [PRISM_N, 1.0, nf, sub_n]
    z = [kz(n, N, k) for n in layers]
    q = [zj / (n * n) if pol == "TM" else zj for zj, n in zip(z, layers)]
    if pol == "TE" and film_ny is not None:     # biaxial film: TE light is polarised along y
        ny = film_ny + 1j * film_kappa
        w = np.sqrt(ny * ny - N * N + 0j)
        z[2] = k * np.where(np.imag(w) < 0, -w, w)
        q[2] = z[2]
    if pol == "TM":
        w = np.sqrt(ne * ne - N * N + 0j)
        w = np.where(np.imag(w) < 0, -w, w)
        z[2] = k * nf / ne * w
        q[2] = z[2] / (nf * nf)
    r = lambda i, j: (q[i] - q[j]) / (q[i] + q[j])
    phase_f = np.exp(2j * z[2] * film_um)
    lower = (r(1, 2) + r(2, 3) * phase_f) / (1 + r(1, 2) * r(2, 3) * phase_f)
    phase_g = np.exp(2j * z[1] * gap_um)
    total = (r(0, 1) + lower * phase_g) / (1 + r(0, 1) * lower * phase_g)
    return np.abs(total) ** 2
