"""Forward model of the NATA / acrylamide quenching experiment (generator and reference; identical copy in solution/).

Quenching (Smoluchowski): the acrylamide molecules are placed at random (Poisson, number density n) around each
NATA at the moment of excitation and diffuse relative to it with the mutual diffusion coefficient D; one that
reaches the encounter distance R quenches the excited NATA at once (perfectly absorbing sink). For independent
quenchers the survival probability of an excited NATA is exp(-n * (volume swept)), which gives
    no quencher inside R at excitation (otherwise quenched at once):   exp(-(4/3) pi R^3 n)
    survival of the others:          I(t) = exp(-t/tau0 - 4 pi R D n t - 8 R^2 n sqrt(pi D t))
with n = [Q] * 6.02214076e20 cm^-3 per M, R in cm and D in cm^2/s.

Lifetimes: time-correlated single-photon counting at 10 MHz. In each excitation period the photons that would
reach the detector are Poisson; the electronics time only the first one of a period. If lam_i is the mean number of
photons per period in channel i, the expected histogram is
    N_exc * (1 - exp(-lam_i)) * exp(-sum_{j<i} lam_j)                       (first-photon / pile-up distortion).
lam = A * (IRF (x) I) / int I + s * IRF + dark. The convolution treats the IRF as piecewise constant within bins of
width d and integrates the decay exactly over each bin (the double integral is the second difference of the second
antiderivative of I, evaluated by Gauss-Legendre quadrature in sqrt(t), where I is smooth); `shift` delays the
emission against the IRF.

Steady state: the excitation beam crosses the 4 mm path of the cuvette; the emission is seen from depths
x0..x1 (cm) along that path. The intensity is proportional to the NATA concentration times
    g(alpha) = mean over x0..x1 of 10**(-alpha x),   alpha = eps_Q c_Q + eps_F c_F   (absorbance per cm at 295 nm)
times the relative quantum yield exp(-(4/3) pi R^3 n) * int_0^inf I(t) dt / tau0.
"""
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erfc, erfcx, ndtr

NCH = 4096
DT = 50.0 / NCH       # ns per channel as set on the TAC (50 ns range over 4096 channels); the calibrated width differs
T = (np.arange(NCH) + 0.5) * DT
N_PER_M = 6.02214076e20   # molecules per cm^3 in a 1 M solution
_GX, _GW = np.polynomial.legendre.leggauss(10)
_GX, _GW = (_GX + 1) / 2, _GW / 2


def coefficients(q, tau0, r_a, d):
    """decay rate k (ns^-1), transient coefficient beta (ns^-1/2) and static exponent at quencher concentration q (M);
    r_a in Angstrom, d in cm^2/s"""
    n = np.asarray(q, float) * N_PER_M
    r = r_a * 1e-8
    k = 1.0 / tau0 + 4 * np.pi * r * d * n * 1e-9
    beta = 8 * r ** 2 * n * np.sqrt(np.pi * d * 1e-9)
    return k, beta, 4.0 / 3.0 * np.pi * r ** 3 * n


def survival(t, k, beta):
    t = np.maximum(t, 0.0)
    return np.exp(-k * t - beta * np.sqrt(t))


def mean_time(k, beta):
    """int_0^inf exp(-k t - beta sqrt(t)) dt (ns)"""
    x = beta / (2 * np.sqrt(k))
    return (1 - np.sqrt(np.pi) * x * erfcx(x)) / k


def relative_yield(q, tau0, r_a, d):
    """F / F0 of the NATA at quencher concentration q, before inner filter and dilution"""
    k, beta, stat = coefficients(q, tau0, r_a, d)
    return np.exp(-stat) * mean_time(k, beta) / tau0


def _bin_integrals(a, b, k, beta):
    """int_a^b f and int_a^b (b - u) f(u) du for f = exp(-k u - beta sqrt(u)), 0 <= a <= b, quadrature in sqrt(u)"""
    sa, sb = np.sqrt(a), np.sqrt(b)
    w = sa[:, None] + (sb - sa)[:, None] * _GX
    u = w * w
    f = np.exp(-k * u - beta * w) * (2 * w * (sb - sa)[:, None] * _GW)
    return f.sum(1), (f * (b[:, None] - u)).sum(1)


def kernel(k, beta, shift, d, n):
    """(m_lo, K): K[p] is the decay integrated over bin j + m_lo + p per unit IRF mass spread uniformly over bin j,
    the emission delayed by `shift` (ns); bins of width d"""
    m_lo = int(np.floor(shift / d)) - 1
    m = np.arange(m_lo - 1, n)
    t = m * d - shift
    i0, i1 = _bin_integrals(np.clip(t, 0, None), np.clip(t + d, 0, None), k, beta)
    return m_lo, i0[:-1] + np.diff(i1) / d


def convolve_decay(h, k, beta, shift, d=DT):
    """IRF masses h (bins of width d) convolved with the decay, integrated over each bin"""
    n = len(h)
    m_lo, kern = kernel(k, beta, shift, d, n)
    c = fftconvolve(h, kern)
    out = np.zeros(n)
    idx = np.arange(n) - m_lo
    ok = (idx >= 0) & (idx < len(c))
    out[ok] = c[idx[ok]]
    return out


def emg_cdf(x, sigma, tau):
    """CDF of a Gaussian (sigma) convolved with an exponential tail (tau); x from the Gaussian centre"""
    x = np.asarray(x, float)
    z = (sigma / tau - x / sigma) / np.sqrt(2)
    tail = np.empty_like(x)
    pos = z > 0
    tail[pos] = 0.5 * np.exp(-x[pos] ** 2 / (2 * sigma ** 2)) * erfcx(z[pos])
    tail[~pos] = 0.5 * np.exp(sigma ** 2 / (2 * tau ** 2) - x[~pos] / tau) * erfc(z[~pos])
    return ndtr(x / sigma) - tail


def irf_masses(edges, shift=0.0, sigma=0.0175, tau_irf=0.028, t0=1.5):
    """probability of the instrument response in each bin between `edges` (ns)"""
    return np.diff(emg_cdf(np.asarray(edges) - t0 - shift, sigma, tau_irf))


def first_photon(lam, n_exc):
    """expected TCSPC histogram when only the first photon of each period is timed"""
    before = np.concatenate([[0.0], np.cumsum(lam)[:-1]])
    return n_exc * (1 - np.exp(-lam)) * np.exp(-before)


def coates(counts, n_exc):
    """Coates correction: mean photons per period in each channel from a first-photon histogram"""
    before = np.concatenate([[0.0], np.cumsum(counts)[:-1]])
    return -np.log1p(-counts / (n_exc - before))


def window_mean(alpha, x0=0.10, x1=0.30):
    """mean of 10**(-alpha x) over depths x0..x1 (cm)"""
    alpha = np.asarray(alpha, float)
    k = alpha * np.log(10)
    safe = np.where(k > 1e-12, k, 1.0)
    return np.where(k > 1e-12, (np.exp(-k * x0) - np.exp(-k * x1)) / (safe * (x1 - x0)), 1.0)
