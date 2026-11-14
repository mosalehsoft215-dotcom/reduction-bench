"""Linear optics of the focus-scan line, shared by the generator and the reference.

The electrons come from a small source (rms size `size`, rms divergence `div`) that lies a distance
`dist` upstream of the entrance of the active plasma lens. The lens (length LENS, field linear in the
radius, gradient G) is followed by a drift DRIFT to the profile screen. An electron with relative
momentum deviation delta sees k = 0.299792458 G[T/m] / p0[GeV/c] / (1 + delta). With (M11, M12) of
lens plus drift, the rms size of a momentum slice at the screen is
    sigma^2 = M11^2 size^2 + (M11 dist + M12)^2 div^2,
and the source is taken to be the same for every slice. A beam whose centroid enters the lens at
(x, x') from the lens axis is steered: the centroid of a slice lands at M11 x + M12 x', which depends
on the slice's momentum through k. At the screen every profile is blurred by the optics (a Gaussian) and
by the depth of the crystal (a box: each electron is seen as a line, see lens_fit.py).
"""
import numpy as np
from scipy.special import ndtr

ME = 0.51099895                       # electron rest energy, MeV
C_GEV = 0.299792458                   # k = C_GEV * G[T/m] / p[GeV/c]
LENS, DRIFT = 0.200, 1.500            # m


def momentum_gev(t_mev):
    return np.sqrt(t_mev ** 2 + 2 * t_mev * ME) / 1000.0


def m11_m12(k):
    """thick focusing lens (k > 0) followed by the drift; returns (M11, M12)"""
    k = np.asarray(k, float)
    sk = np.sqrt(k)
    c, s, sp = np.cos(sk * LENS), np.sin(sk * LENS) / sk, -sk * np.sin(sk * LENS)
    return c + DRIFT * sp, s + DRIFT * c


def size2(k, size, div, dist):
    m11, m12 = m11_m12(k)
    return m11 ** 2 * size ** 2 + (m11 * dist + m12) ** 2 * div ** 2


def hermite(n=31):
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


def skew_alpha(skew):
    """shape parameter of the skew-normal distribution with the given skewness (|skew| < 0.99)"""
    from scipy.optimize import brentq
    if skew == 0:
        return 0.0
    f = lambda d: (4 - np.pi) / 2 * (d * np.sqrt(2 / np.pi)) ** 3 / (1 - 2 * d * d / np.pi) ** 1.5 - abs(skew)
    d = brentq(f, 1e-9, 1 - 1e-12)
    return np.sign(skew) * d / np.sqrt(1 - d * d)


def skew_nodes(skew, n=161, span=7.0):
    """quadrature nodes and weights of the standardised (zero mean, unit variance) skew-normal distribution"""
    if skew == 0:
        return hermite(41)
    a = skew_alpha(skew)
    b = a / np.sqrt(1 + a * a) * np.sqrt(2 / np.pi)
    om = 1 / np.sqrt(1 - b * b)
    xi = -om * b
    x = np.linspace(-span, span, n)
    t = (x - xi) / om
    w = np.exp(-t * t / 2) * ndtr(a * t)
    w = w / w.sum()
    m = np.sum(w * x)
    sd = np.sqrt(np.sum(w * (x - m) ** 2))
    return (x - m) / sd, w


def centre(k, entry):
    """screen position of the centroid ray of a beam that enters the lens at entry = (x, x') measured
    from the lens axis (up to a constant that does not depend on k)"""
    m11, m12 = m11_m12(k)
    return m11 * entry[0] + m12 * entry[1]


def smeared(x, c, s2, wq, box):
    """sum over momentum slices (weights wq) of a Gaussian of variance s2 centred on c, convolved with a
    box of full width `box` (the line of light of one electron in the crystal, see lens_fit.py)"""
    s = np.sqrt(s2)[:, None]
    d = x[None, :] - c[:, None]
    if box <= 0:
        return np.sum(wq[:, None] * np.exp(-d ** 2 / (2 * s ** 2)) / (np.sqrt(2 * np.pi) * s), axis=0)
    return np.sum(wq[:, None] * (ndtr((d + box / 2) / s) - ndtr((d - box / 2) / s)), axis=0) / box


def profile(xc, k0, size, div, dist, spread, res, n=31, entry=(0.0, 0.0), box=0.0, nodes=None, kp=0.0):
    """expected profile of one shot (unit area, per metre) at positions xc measured from the shot's
    centroid, for nominal k0: a mixture over a momentum distribution of rms `spread` (Gaussian unless
    `nodes` give its standardised shape) of Gaussians blurred by the optics (`res`) and by the depth of
    the crystal (`box`); every momentum slice sits on its own centroid ray. kp is an extra focusing that the
    beam's own wake in the discharge plasma adds; it acts about the beam's own centroid, so it changes the
    sizes but does not steer the centroid rays"""
    xq, wq = (nodes if nodes is not None else hermite(n)) if spread > 0 else (np.zeros(1), np.ones(1))
    k = k0 / (1 + xq * spread)
    s2 = size2(k + kp / (1 + xq * spread), size, div, dist) + res ** 2
    c = centre(k, entry)
    return smeared(np.asarray(xc, float), c - np.sum(wq * c), s2, wq, box)


def run_profile(xc, k0, size, div, dist, within, jitter, res, entry=(0.0, 0.0), n=21, nj=7, box=0.0):
    """average of the centred shots of one lens current: the shots differ in mean momentum (Gaussian, rms
    `jitter`), each is the mixture of `profile` over its own momentum distribution (rms `within`) and is
    shifted to its own centroid before the average"""
    xj, wj = hermite(nj) if jitter > 0 else (np.zeros(1), np.ones(1))
    xq, wq = hermite(n) if within > 0 else (np.zeros(1), np.ones(1))
    x = np.asarray(xc, float)
    out = np.zeros_like(x)
    for dj, wshot in zip(xj * jitter, wj):
        k = k0 / (1 + dj + xq * within)
        s2 = size2(k, size, div, dist) + res ** 2
        c = centre(k, entry)
        out += wshot * smeared(x, c - np.sum(wq * c), s2, wq, box)
    return out
