"""An independently written analysis (no code shared with solution/):

  step              reference (solution/kidfit.py)                this script
  resonance fits    all sweeps together, line and coupling shared  every sweep on its own (delay removed
                                                                   from the off-resonance phase slope first),
                                                                   1/Qi errors propagated from Ql, |Qc|, phi
  conductivity      Mattis-Bardeen integrals, adaptive quadrature  the same integrals after the substitutions
                    with algebraic end-point weights               E = D cosh t and E = D - hf/2 + (hf/2) cos t,
                                                                   200-point Gauss-Legendre
  film              full surface impedance coth(k d)               linear response, gamma = [1 + 2u/sinh(2u)]/2
  final fit         fr and 1/Qi jointly                            the same, on a log-scaled parameter set
usage: python independent_correct.py DATA_DIR [variant]
       variant: full; gao (closed low-temperature forms of the conductivities); thin (gamma = 1);
       nophi (symmetric resonance); hight (no TLS, T >= 230 mK only); bcs (gap from Tc)
"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
from scipy.special import i0, k0, psi

KB, H = 86.17333, 4.135667696e-9
TC = 1.312
U = 62.0 / 57.0
XG, WG = np.polynomial.legendre.leggauss(200)


def fit_one(f, z, symmetric=False):
    n = len(f) // 8
    e = np.r_[np.arange(n), np.arange(len(f) - n, len(f))]
    tau = -np.polyfit(f[e], np.unwrap(np.angle(z))[e], 1)[0] / (2 * np.pi)
    zc = z * np.exp(2j * np.pi * f * tau)
    env = zc[e].mean()
    zn = zc / env
    k = np.argmax(np.abs(1 - zn))
    fr0 = f[k]
    width = np.ptp(f[np.abs(1 - zn) > np.abs(1 - zn[k]) / np.sqrt(2)])
    ql0 = fr0 / max(width, f[1] - f[0])

    def model(p):
        a, al, t, dfr, lql, lqc, phi = p
        if symmetric:
            phi = 0.0
        ql, qc = np.exp(lql), np.exp(lqc)
        return a * np.exp(1j * al - 2j * np.pi * f * (tau + t * 1e-9)) * (1 - ql / qc * np.exp(1j * phi) / (1 + 2j * ql * (f / (fr0 + dfr * 1e3) - 1)))

    def res(p):
        d = model(p) - z
        return np.r_[d.real, d.imag]
    p0 = [abs(env), np.angle(env), 0.0, 0.0, np.log(ql0), np.log(ql0 / max(abs(1 - zn[k]), 1e-3)), 0.0]
    r = least_squares(res, p0, method="lm", x_scale="jac", xtol=1e-15, ftol=1e-15, max_nfev=20000)
    a, al, t, dfr, lql, lqc, phi = r.x
    if symmetric:
        phi = 0.0
    cov = np.linalg.pinv(r.jac.T @ r.jac) * np.sum(r.fun ** 2) / (r.fun.size - 7)
    ql, qc = np.exp(lql), np.exp(lqc)
    qi_inv = 1 / ql - np.cos(phi) / qc
    g = np.zeros(7)
    g[4], g[5], g[6] = -1 / ql, np.cos(phi) / qc, 0.0 if symmetric else np.sin(phi) / qc
    return fr0 + dfr * 1e3, np.sqrt(cov[3, 3]) * 1e3, qi_inv, np.sqrt(g @ cov @ g)


def sigma(T, D, hf):
    kT = KB * T
    fe = lambda E: 0.5 * (1 - np.tanh(E / (2 * kT)))
    tm = np.arccosh(1 + 60 * kT / D)
    t = 0.5 * tm * (XG + 1)
    E = D * np.cosh(t)
    s1 = 2 / hf * 0.5 * tm * np.sum(WG * (E * E + D * D + hf * E) / np.sqrt((E + hf) ** 2 - D * D) * (fe(E) - fe(E + hf)))
    th = 0.5 * np.pi * (XG + 1)
    E = D - hf / 2 + hf / 2 * np.cos(th)
    s2 = 1 / hf * 0.5 * np.pi * np.sum(WG * (E * E + D * D + hf * E) / (np.sqrt(D + E) * np.sqrt(E + hf + D)) * (1 - 2 * fe(E + hf)))
    return s1, s2


def response(T, D, f0, fd, alpha, qo, tls=True, closed=False, gamma=None):
    hf = H * f0
    kT = KB * T
    if gamma is None:
        gamma = 0.5 * (1 + 2 * U / np.sinh(2 * U))
    if closed:
        xi = hf / (2 * kT)
        ds2 = -2 * np.exp(-D / kT) * np.exp(-xi) * i0(xi)
        s1_s2 = (4 / np.pi) * np.exp(-D / kT) * np.sinh(xi) * k0(xi) / (1 + ds2)
    else:
        s20 = sigma(0.01, D, hf)[1]
        s1, s2 = sigma(T, D, hf)
        ds2, s1_s2 = s2 / s20 - 1, s1 / s2
    x = hf / (2 * np.pi * kT)
    tshift = fd / np.pi * (np.real(psi(0.5 + x / 1j)) - np.log(x)) if tls else 0.0
    tloss = fd * np.tanh(hf / (2 * kT)) if tls else 0.0
    return f0 * (1 + 0.5 * alpha * gamma * ds2 + tshift), 1 / qo + alpha * gamma * s1_s2 + tloss


if __name__ == "__main__":
    d = sys.argv[1]
    var = sys.argv[2] if len(sys.argv) > 2 else "full"
    A = np.loadtxt(os.path.join(d, "sweeps.csv"), delimiter=",", skiprows=1)
    Ts = np.unique(A[:, 0])
    out = np.array([fit_one(A[A[:, 0] == t, 1], A[A[:, 0] == t, 2] + 1j * A[A[:, 0] == t, 3], symmetric=(var == "nophi")) for t in Ts])
    fr, sfr, qi, sqi = out.T
    T = Ts / 1000.0
    use = T >= 0.23 if var == "hight" else np.ones(T.size, bool)
    tls = var != "hight"
    fixD = 1.764 * KB * TC if var == "bcs" else None
    kw = dict(tls=tls, closed=(var == "gao"), gamma=1.0 if var == "thin" else None)

    def res(p):
        f0, la, D, fd, lqo = p
        D = fixD if fixD is not None else D
        m = [response(t, D, f0 * 1e9, fd * 1e-5, np.exp(la), np.exp(lqo) * 1e5, **kw) for t in T[use]]
        mf, mq = np.array(m).T
        return np.r_[(mf - fr[use]) / sfr[use], (mq - qi[use]) / sqi[use]]
    best = None
    for D0 in (190.0, 215.0):
        p0 = [fr[0] / 1e9, np.log(0.07), D0, 1.0 if tls else 0.0, np.log(3.0)]
        lo = [p0[0] - 0.01, np.log(1e-4), 120, 0.0, np.log(0.01)]
        hi = [p0[0] + 0.01, 0.0, 400, 100.0 if tls else 1e-9, np.log(1000)]
        r = least_squares(res, p0, bounds=(lo, hi), x_scale=[1e-6, 0.05, 1.0, 0.1, 0.1], xtol=1e-14, ftol=1e-14)
        best = r if best is None or r.cost < best.cost else best
    r = best
    D = fixD if fixD is not None else r.x[2]
    alpha = np.exp(r.x[1])
    gt = os.path.join(d, "..", "ground_truth.json")
    tr = json.load(open(gt)) if os.path.exists(gt) else dict(Delta0_ueV=np.nan, alpha_k=np.nan)
    print("%s %-5s Delta0 %.3f ueV (%+.2f %%)  alpha %.5f (%+.2f %%)  chi2 %.1f / %d"
          % (d, var, D, 100 * (D / tr["Delta0_ueV"] - 1), alpha, 100 * (alpha / tr["alpha_k"] - 1), np.sum(r.fun ** 2), r.fun.size), flush=True)
