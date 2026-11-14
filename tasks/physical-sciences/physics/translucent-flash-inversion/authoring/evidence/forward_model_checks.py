"""Forward-model checks for the v3 generator (no fitting).

1. The closed-form through-thickness response against an independent finite-volume solve in z.
2. In-plane convergence: the generator basis (M=24, P=24) and the reference basis (M=16, P=16) against M=32, P=32,
   at the generator inputs, in units of the per-row noise sigma.
Writes forward_model_checks.json.
"""
import json
from pathlib import Path
import numpy as np
from scipy.special import jv, jnp_zeros, ive
from numpy.polynomial.legendre import leggauss

R = 5.0
PROG = {1: (780., .70, .55, 4.0), 2: (940., .92, 1.45, 71.0), 3: (1064., .78, 2.20, 149.0)}
MASKS = {'M1': (.75, 3.20, 8., 72.), 'M2': (1.15, 4.10, 132., 66.), 'M3': (1.45, 4.45, 244., 61.)}
THK = [.52, .68, .89, 1.16, 1.46]; DRV = [.18, .55, 1.00, 1.55]; FREQS = np.geomspace(.045, 42., 14)


def gl(a, b, n):
    x, w = leggauss(n); return .5 * (b - a) * x + .5 * (b + a), .5 * (b - a) * w


def basis(M, P, nq=600):
    ms, ks = [0], [0.]
    for m in range(M + 1):
        for root in jnp_zeros(m, P):
            if m == 0: ms.append(0); ks.append(root / R)
            else: ms += [m, -m]; ks += [root / R] * 2
    ms = np.array(ms); ks = np.array(ks); am = np.abs(ms)
    norm = np.where(ks == 0, np.pi * R * R, np.pi * R * R * (1 - (am / np.where(ks == 0, 1, ks * R)) ** 2) * jv(am, ks * R) ** 2)
    avg = {}
    for mk, (ri, ro, c, wd) in MASKS.items():
        r, w = gl(ri, ro, nq); rad = (jv(am[:, None], ks[:, None] * r[None, :]) * r * w).sum(1)
        h = np.deg2rad(wd / 2); ang = np.where(am == 0, 2 * h, 2 * np.sin(ms * h) / np.where(ms == 0, 1, ms)) * np.exp(1j * ms * np.deg2rad(c))
        avg[mk] = rad * ang
    src = {}
    r, w = gl(0, R, 2 * nq)
    for p, (wl, wb, d, ang) in PROG.items():
        g = np.exp(-2 * (r - d) ** 2 / wb ** 2); x = 4 * r * d / wb ** 2
        Z = 2 * np.pi * np.sum(g * ive(0, x) * r * w)
        rad = np.array([np.sum(g * ive(m, x) * jv(m, k * r) * r * w) for m, k in zip(am, ks)])
        src[p] = 2 * np.pi * rad / (Z * norm)
    return dict(m=ms, k=ks, avg=avg, src=src)


def slab(mu, Rr, L):
    e = np.exp(-mu * L); T = (1 - Rr) ** 2 * e / (1 - Rr ** 2 * e * e); A = (1 - Rr) * (1 - e) / (1 - Rr * e); return T, A


def axial(L, mu, Rf, kap, Ri, alpha, sv, b):
    """detector-band-weighted temperature for unit incident flash power, exact in z"""
    g = np.sqrt(sv / alpha)
    eL = np.exp(-mu * L); emL = np.exp(mu * L); D = 1 - Rf ** 2 * eL ** 2
    a1 = mu * (1 - Rf) / D; a2 = a1 * Rf * eL ** 2
    den = alpha * (g * g - mu * mu); A1 = a1 / den; A2 = a2 / den
    th0 = A1 + A2; dth0 = mu * (A2 - A1); thL = A1 * eL + A2 * emL; dthL = mu * (A2 * emL - A1 * eL)
    E = np.exp(-g * L); Pp = -(alpha * g + b); Qq = (alpha * g - b) * E
    r1 = b * th0 - alpha * dth0; r2 = b * thL + alpha * dthL; det = Pp * Pp - Qq * Qq
    C1 = (Pp * r1 - Qq * r2) / det; C2 = (Pp * r2 - Qq * r1) / det
    ek = np.exp(-kap * L); K = kap * (1 - Ri) / (1 - Ri ** 2 * ek ** 2)
    T1 = (eL - ek) / (kap - mu) + Ri * ek * (1 - np.exp(-(kap + mu) * L)) / (kap + mu)
    T2 = (emL - ek) / (kap + mu) + Ri * ek * (1 - np.exp((mu - kap) * L)) / (kap - mu)
    T3 = (E - ek) / (kap - g) + Ri * ek * (1 - np.exp(-(kap + g) * L)) / (kap + g)
    T4 = (1 - np.exp(-(kap + g) * L)) / (kap + g) + Ri * ek * (ek - E) / (g - kap)
    return K * (A1 * T1 + A2 * T2 + C1 * T3 + C2 * T4)


def fd_axial(L, mu, Rf, kap, Ri, alpha, s, b, N=4000):
    """finite-volume check of axial(): unit control volumes, coating capacity/surface loss on the end faces"""
    z = (np.arange(N) + .5) * L / N; h = L / N
    eL = np.exp(-mu * L); S = mu * (1 - Rf) / (1 - Rf ** 2 * eL ** 2) * (np.exp(-mu * z) + Rf * np.exp(-mu * (2 * L - z)))
    ek = np.exp(-kap * L); w = kap * (1 - Ri) / (1 - Ri ** 2 * ek ** 2) * (np.exp(-kap * (L - z)) + Ri * np.exp(-kap * (L + z)))
    main = np.full(N, -2 * alpha / h ** 2 - s, complex); off = np.full(N - 1, alpha / h ** 2, complex)
    # face: flux alpha*(T1-Tface)/(h/2) = b*Tface -> Tface = 2a/(2a+bh) T1 ; flux into cell 1 = -b*Tface
    gface = b * 2 * alpha / (2 * alpha + b * h)
    main[0] = -alpha / h ** 2 - s - gface / h; main[-1] = main[0]
    import scipy.linalg as sl
    ab = np.zeros((3, N), complex); ab[0, 1:] = off; ab[1] = main; ab[2, :-1] = off
    T = sl.solve_banded((1, 1), ab, -S)
    return np.sum(w * T) * h


TRUE = dict(alpha=.71, pen=.79, loss=.105, tau=.031, alpha_r=1.31, surf=.016, coat=.060, aslope=.24, lheat=.035,
            phi={'A': 9., 'B': -17., 'C': 21.}, pexp=.68, Rf=.112, kap=6.35, Ri=.09, leak=.0048)


def thermal_table(B, P, lumped=False, unitnorm=False, qmode='A'):
    """complex thermal response for every (L, program, drive, freq) x (session, mask): dict[(s,mk)] -> array (5,3,4,14)"""
    mu0 = 1 / P['pen']; out = {}
    Ls = np.array(THK)[:, None, None, None, None]; om = 2 * np.pi * FREQS[None, None, None, :, None]
    dls = np.array(DRV)[None, None, :, None, None]
    res = np.zeros((5, 3, 4, 14, len(B['m'])), complex)
    for ip, p in enumerate(PROG):
        mu = mu0 * (940. / PROG[p][0]) ** P['pexp']
        T, A = slab(mu, P['Rf'], Ls)
        q = dls * (A if qmode == 'A' else 1 - T)
        al = P['alpha'] * (1 + P['aslope'] * q); lv = P['loss'] + P['lheat'] * q * q
        kk = P['alpha_r'] * B['k'][None, None, None, None, :] ** 2
        if lumped:
            sv = lv + 2 * P['surf'] / Ls + 1j * om * (1 + 2 * P['coat'] / Ls) + kk; b = 0 * om
        else:
            sv = lv + 1j * om + kk; b = P['surf'] + 1j * om * P['coat']
        th = axial(Ls, mu, P['Rf'], P['kap'], P['Ri'], al, sv, b)
        if unitnorm:
            ek = np.exp(-P['kap'] * Ls); Eir = (1 - P['Ri']) * (1 - ek) / (1 - P['Ri'] * ek); th = th / (A * Eir)
        res[:, ip] = th[:, 0]
    for s in 'ABC':
        for mk in MASKS:
            v = np.stack([B['src'][p] * np.exp(-1j * B['m'] * np.deg2rad(PROG[p][3] + P['phi'][s])) * B['avg'][mk] for p in PROG])
            out[(s, mk)] = np.einsum('lpdfj,pj->lpdf', res, v)
    return out


if __name__ == '__main__':
    rng = np.random.default_rng(1); worst = 0
    for _ in range(24):
        L = rng.choice(THK); mu = rng.uniform(1, 1.5); al = rng.uniform(.6, .9); om = 2*np.pi*rng.choice(FREQS)
        s = .12+1j*om+1.31*rng.uniform(0, 12)**2; b = .016+1j*om*.06
        a = axial(L, mu, .112, 6.35, .09, al, s, b); f = fd_axial(L, mu, .112, 6.35, .09, al, s, b)
        worst = max(worst, abs(a-f)/abs(a))
    dl = np.array(DRV)[None, None, :, None]
    ref = thermal_table(basis(32, 32, 900), TRUE)
    def dev(o):
        return max(float(np.max(np.abs(dl*(o[k]-ref[k]))/(7e-5+3.3e-4*np.abs(dl*ref[k])))) for k in ref)
    out = {'axial_closed_form_vs_finite_volume_max_relative_difference': worst,
           'in_plane_max_deviation_in_sigma_vs_M32_P32': {'generator_M24_P24': dev(thermal_table(basis(24, 24, 900), TRUE)),
                                                           'reference_M16_P16': dev(thermal_table(basis(16, 16, 500), TRUE)),
                                                           'coarse_M8_P8': dev(thermal_table(basis(8, 8, 500), TRUE))}}
    (Path(__file__).with_name('forward_model_checks.json')).write_bytes(json.dumps(out, indent=2).encode()); print(out)
