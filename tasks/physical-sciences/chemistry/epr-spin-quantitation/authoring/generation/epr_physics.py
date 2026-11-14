"""Forward model used to write every spectrum of the task (the truth).

Rigid-limit powder spectra of S = 1/2 centres, computed exactly:
  * for each orientation of the field in the molecular frame the full spin Hamiltonian
    H(B) = muB B n.g.S + S.A.I - gn muN B n.I is diagonalised and every electron-spin transition
    (allowed and forbidden) is located at the fixed microwave frequency by Newton steps on B;
  * the transition probability is |<a|(b.g).S muB - gn muN b.I|b>|^2 averaged over the two
    directions b of the linearly polarised B1 perpendicular to B0, divided by |d(Ea-Eb)/dB|
    (the field-swept Jacobian).  This is the physics behind the Aasa-Vanngard g-average.
  * each line is a Gaussian in field whose width depends on orientation and on its hyperfine
    displacement (strain); natural copper is 69.15 % 63Cu and 30.85 % 65Cu.
Detection: first harmonic of 100 kHz field modulation (exact Bessel kernel), resonator profile
along the tube axis, Bloch saturation of narrow spin packets (area factor 1/sqrt(1+s)), Curie law.
"""
import numpy as np
from scipy.special import erf, j1

MUB = 1.39962449361      # MHz / G
MUN = 7.6225932291e-4    # MHz / G
H_OVER_K = 4.799243073e-5  # K / MHz  (h/k)


def _ops(I):
    half = np.array([[0.0, 0.5], [0.5, 0.0]]), np.array([[0.0, -0.5j], [0.5j, 0.0]]), np.diag([0.5, -0.5])
    m = np.arange(I, -I - 1, -1)
    n = len(m)
    up = np.zeros((n, n))
    for k in range(1, n):
        up[k - 1, k] = np.sqrt(I * (I + 1) - m[k] * (m[k] + 1))
    nuc = (up + up.T) / 2, (up - up.T) / 2j, np.diag(m).astype(float)
    S = np.array([np.kron(o, np.eye(n)) for o in half])
    Iv = np.array([np.kron(np.eye(2), o) for o in nuc])
    return S, Iv


def orientation_grid(rhombic, n_theta=None, n_phi=None):
    """Gauss-Legendre in cos(theta) on [0,1]; uniform phi on [0, pi/2] (D2h powder symmetry)."""
    if not rhombic:
        n_theta = n_theta or 320
        x, w = np.polynomial.legendre.leggauss(n_theta)
        ct, wt = (x + 1) / 2, w / 2
        st = np.sqrt(1 - ct ** 2)
        return np.stack([st, np.zeros_like(st), ct], 1), wt
    n_theta = n_theta or 96
    n_phi = n_phi or 48
    x, w = np.polynomial.legendre.leggauss(n_theta)
    ct, wt = (x + 1) / 2, w / 2
    ph = (np.arange(n_phi) + 0.5) * (np.pi / 2) / n_phi
    st = np.sqrt(1 - ct ** 2)
    v = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(ct, np.ones(n_phi))], -1).reshape(-1, 3)
    return v, np.repeat(wt, n_phi) / n_phi


def resonances(g, A, I, gn, nu, n_vec):
    """All electron-spin transitions at frequency nu (MHz) for field directions n_vec.
    Returns fields (No, P), weights W/|dE/dB| (No, P) in units where isotropic g, A=0 gives g/2 per pair."""
    S, Iv = _ops(I)
    d = S.shape[1]
    half = d // 2
    g = np.asarray(g, float)
    A = np.asarray(A, float)
    F = sum(A[k] * S[k] @ Iv[k] for k in range(3))
    Z = np.einsum('nk,kij->nij', n_vec * g * MUB, S) - gn * MUN * np.einsum('nk,kij->nij', n_vec, Iv)
    geff = np.sqrt(((n_vec * g) ** 2).sum(1))
    pairs = [(a, b) for a in range(half, d) for b in range(half)]
    No, P = len(n_vec), len(pairs)
    B = np.repeat((nu / (geff * MUB))[:, None], P, 1)
    ia = np.array([p[0] for p in pairs])
    ib = np.array([p[1] for p in pairs])
    for _ in range(40):
        H = Z[:, None] * B[..., None, None] + F
        E, V = np.linalg.eigh(H.reshape(-1, d, d))
        E = E.reshape(No, P, d)
        V = V.reshape(No, P, d, d)
        dZ = np.einsum('npia,nij,npja->npa', V.conj(), Z, V).real
        dE = E[:, np.arange(P), ia] - E[:, np.arange(P), ib]
        slope = dZ[:, np.arange(P), ia] - dZ[:, np.arange(P), ib]
        step = (dE - nu) / slope
        B = B - step
        if np.max(np.abs(step)) < 1e-8:
            break
    t1 = np.cross(n_vec, [0.0, 0.0, 1.0])
    flat = np.linalg.norm(t1, axis=1) < 1e-9
    t1[flat] = [1.0, 0.0, 0.0]
    t1 /= np.linalg.norm(t1, axis=1)[:, None]
    t2 = np.cross(n_vec, t1)
    W = np.zeros((No, P))
    Vt = V.transpose(0, 1, 3, 2)            # (orientation, pair, eigenstate, component)
    va = Vt[:, np.arange(P), ia, :]
    vb = Vt[:, np.arange(P), ib, :]
    for t in (t1, t2):
        M = np.einsum('nk,kij->nij', t * g * MUB, S) - gn * MUN * np.einsum('nk,kij->nij', t, Iv)
        amp = np.einsum('npi,nij,npj->np', va.conj(), M, vb)
        W += 0.5 * np.abs(amp) ** 2
    weight = W / np.abs(slope) / MUB * 4.0 / (2 * I + 1)
    return B, weight


def absorption(centre, nu, field):
    """Absorption per spin on the field grid (G); its integral over field is the centre's
    intensity factor (equal to the Aasa-Vanngard g-average to 0.02 %)."""
    rhombic = len(set(np.round(centre['g'][:2], 6))) > 1 or len(set(np.round(centre['A'][:2], 6))) > 1
    n_vec, w_or = orientation_grid(rhombic)
    edges = np.concatenate([[field[0] - (field[1] - field[0]) / 2], (field[1:] + field[:-1]) / 2,
                            [field[-1] + (field[-1] - field[-2]) / 2]])
    out = np.zeros(len(field))
    for frac, ascale, gn in centre['isotopes']:
        parts = [resonances(centre['g'], np.asarray(centre['A']) * ascale, centre['I'], gn, nu, n_vec[k:k + 512])
                 for k in range(0, len(n_vec), 512)]
        B = np.concatenate([q[0] for q in parts])
        wt = np.concatenate([q[1] for q in parts])
        wt = wt * w_or[:, None] * frac
        centre_field = (B * wt).sum(1) / wt.sum(1)
        shift = B - centre_field[:, None]
        strain = np.sqrt((centre['strain_G'] ** 2 * n_vec ** 2).sum(1))
        sig = np.sqrt(centre['sigma0_G'] ** 2 + strain[:, None] ** 2 + (centre['kappa'] * shift) ** 2)
        b = B.ravel()
        s = sig.ravel()
        a = wt.ravel()
        keep = a > 1e-12 * a.max()
        b, s, a = b[keep], s[keep], a[keep]
        for k in range(0, len(b), 1500):
            bb, ss, aa = b[k:k + 1500, None], s[k:k + 1500, None], a[k:k + 1500, None]
            cdf = 0.5 * (1 + erf((edges[None, :] - bb) / (np.sqrt(2) * ss)))
            out += (aa * np.diff(cdf, axis=1)).sum(0)
    return out / np.diff(edges)


def first_harmonic(absorp, step, bm_pp):
    """Exact first-harmonic field-modulation signal (peak-to-peak amplitude bm_pp, G):
    S(B) = (2/pi) int_0^pi a(B + bm/2 cos t) cos t dt; tends to (bm/2) da/dB for small bm."""
    n = len(absorp)
    pad = 1 << int(np.ceil(np.log2(4 * n)))
    k = 2 * np.pi * np.fft.rfftfreq(pad, d=step)
    spec = np.fft.rfft(absorp, pad)
    kernel = 2j * j1(k * bm_pp / 2)
    return np.fft.irfft(spec * kernel, pad)[:n]


def slice_grid(segments, active_mm, n=64):
    """Midpoints and widths of thin slices of the frozen column inside the active length."""
    zs, dz = [], []
    for lo, hi in segments:
        lo2, hi2 = max(lo, -active_mm / 2), min(hi, active_mm / 2)
        if hi2 <= lo2:
            continue
        m = max(4, int(np.ceil(n * (hi2 - lo2) / active_mm)))
        e = np.linspace(lo2, hi2, m + 1)
        zs.append((e[1:] + e[:-1]) / 2)
        dz.append(np.diff(e))
    return np.concatenate(zs), np.concatenate(dz)
