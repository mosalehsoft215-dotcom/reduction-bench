"""Powder lineshapes for the spin-count program (rigid limit, S = 1/2, coaxial rhombic g and A)."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import j1

MU_B = 1.39962449361      # MHz/G
MU_N = 7.6225932291e-4    # MHz/G
CU_ISOTOPES = ((0.6915, 1.0, 1.4824), (0.3085, 1.0711, 1.5878))   # abundance, A scale, nuclear g
N14_ISOTOPES = ((1.0, 1.0, 0.40376),)


def _spin_matrices(I):
    dim_n = int(round(2 * I + 1))
    mz = I - np.arange(dim_n)
    raise_n = np.diag(np.sqrt(I * (I + 1) - mz[1:] * (mz[1:] + 1)), 1)
    nuc = [(raise_n + raise_n.T) / 2, (raise_n - raise_n.T) / 2j, np.diag(mz)]
    pauli = [np.array([[0, 1], [1, 0]]) / 2, np.array([[0, -1j], [1j, 0]]) / 2, np.diag([0.5, -0.5])]
    S = np.stack([np.kron(p, np.eye(dim_n)) for p in pauli])
    In = np.stack([np.kron(np.eye(2), q) for q in nuc])
    return S, In


def g_average(g):
    """Aasa-Vanngard average of a g tensor (principal values)."""
    g = np.asarray(g, float)
    return 2.0 / 3.0 * np.sqrt(np.mean(g ** 2)) + np.mean(g) / 3.0


def _clamped_operator(coarse, fine):
    """Matrix that maps values on the coarse grid to the clamped cubic spline on the fine grid
    (zero slope at both ends: a powder pattern is even about theta = 0, pi/2 and phi = 0, pi/2)."""
    return CubicSpline(coarse, np.eye(len(coarse)), bc_type='clamped', axis=0)(fine)


class Centre:
    """Exact resonance fields and B1-averaged transition rates on a coarse octant grid, mapped by
    clamped bicubic splines onto a fine grid (Gauss-Legendre in cos theta, midpoints in phi).
    Parameters: gx gy gz Ax Ay Az (MHz) w0 w_perp w_par (G) kappa; widths
    sqrt(w0^2 + w_perp^2 (nx^2 + ny^2) + w_par^2 nz^2 + (kappa * hyperfine displacement)^2)."""
    names = ['gx', 'gy', 'gz', 'Ax', 'Ay', 'Az', 'w0', 'w_perp', 'w_par', 'kappa']

    def __init__(self, nu, I, isotopes, coarse=(33, 7), fine=(256, 24)):
        self.nu, self.isotopes = nu, isotopes
        self.S, self.I = _spin_matrices(I)
        d = self.S.shape[1]
        h = d // 2
        self.up = np.repeat(np.arange(h, d), h)
        self.lo = np.tile(np.arange(h), h)
        self.P = h * h
        tc = np.linspace(0, np.pi / 2, coarse[0])
        pc = np.linspace(0, np.pi / 2, coarse[1])
        T, F = np.meshgrid(tc, pc, indexing='ij')
        self.shape_c = T.shape
        self.n_c = np.stack([np.sin(T) * np.cos(F), np.sin(T) * np.sin(F), np.cos(T)], -1).reshape(-1, 3)
        self.e1 = np.stack([np.cos(T) * np.cos(F), np.cos(T) * np.sin(F), -np.sin(T)], -1).reshape(-1, 3)
        self.e2 = np.stack([-np.sin(F), np.cos(F), np.zeros_like(F)], -1).reshape(-1, 3)
        x, w = np.polynomial.legendre.leggauss(fine[0])
        ct = (x + 1) / 2
        th = np.arccos(ct)
        ph = (np.arange(fine[1]) + 0.5) * (np.pi / 2) / fine[1]
        self.Mt = _clamped_operator(tc, th)
        self.Mp = _clamped_operator(pc, ph)
        self.wf = np.outer(w / 2, np.full(fine[1], 1.0 / fine[1]))
        st = np.sqrt(1 - ct ** 2)
        self.nx2 = np.outer(st ** 2, np.cos(ph) ** 2)
        self.ny2 = np.outer(st ** 2, np.sin(ph) ** 2)
        self.nz2 = np.outer(ct ** 2, np.ones(fine[1]))

    def _exact(self, g, A, gn):
        n = self.n_c
        S2 = self.S.reshape(3, -1)
        I2 = self.I.reshape(3, -1)
        d = self.S.shape[1]
        F = sum(A[k] * self.S[k] @ self.I[k] for k in range(3))
        Z = ((n * g * MU_B) @ S2 - gn * MU_N * (n @ I2)).reshape(-1, d, d)
        geff = np.sqrt(((n * g) ** 2).sum(1))
        K = np.sqrt(((A * g * n) ** 2).sum(1)) / geff
        P, up, lo = self.P, self.up, self.lo
        h = d // 2
        m_sum = (up - h - (h - 1) / 2.0) + ((h - 1) / 2.0 - lo)      # m_I(upper) + m_I(lower), first order
        B = (self.nu - np.outer(K, m_sum) / 2.0) / (geff * MU_B)[:, None]
        cols = np.arange(P)
        for _ in range(2):          # from the first-order start, two Newton steps reach 1e-6 G
            E, V = np.linalg.eigh(Z[:, None] * B[..., None, None] + F)
            ZV = Z[:, None] @ V
            dE = np.real(np.sum(V.conj() * ZV, axis=-2))
            gap = E[:, cols, up] - E[:, cols, lo]
            slope = dE[:, cols, up] - dE[:, cols, lo]
            B = B - (gap - self.nu) / slope
        Vt = np.swapaxes(V, 2, 3)              # (orientation, pair, eigenstate, component)
        vu = Vt[:, cols, up, :]
        vl = Vt[:, cols, lo, :]
        rate = np.zeros_like(B)
        for e in (self.e1, self.e2):
            M = ((e * g * MU_B) @ S2 - gn * MU_N * (e @ I2)).reshape(-1, d, d)
            Mv = (M[:, None] @ vl[..., None])[..., 0]
            rate += 0.5 * np.abs(np.sum(vu.conj() * Mv, axis=-1)) ** 2
        rate = rate / np.abs(slope)
        return B.reshape(*self.shape_c, P), rate.reshape(*self.shape_c, P)

    def _to_fine(self, a):
        return np.matmul(self.Mp, np.tensordot(self.Mt, a, (1, 0)))

    def lines(self, p):
        gx, gy, gz, Ax, Ay, Az, w0, wp, wz, kap = p
        g = np.array([gx, gy, gz])
        strain2 = wp ** 2 * (self.nx2 + self.ny2) + wz ** 2 * self.nz2
        Bs, Ws, Ss = [], [], []
        for frac, scale, gn in self.isotopes:
            Bc, Rc = self._exact(g, np.array([Ax, Ay, Az]) * scale, gn)
            B = self._to_fine(Bc)
            R = np.clip(self._to_fine(Rc), 0, None) * self.wf[..., None] * frac
            mean_b = (B * R).sum(-1, keepdims=True) / R.sum(-1, keepdims=True)
            width = np.sqrt(w0 ** 2 + strain2[..., None] + (kap * (B - mean_b)) ** 2)
            Bs.append(B.ravel())
            Ws.append(R.ravel())
            Ss.append(width.ravel())
        return np.concatenate(Bs), np.concatenate(Ws), np.concatenate(Ss)

    def g_avg(self, p):
        return g_average(p[:3])


class Copper(Centre):
    lower = [1.98, 1.98, 2.10, 0.0, 0.0, 50.0, 2.0, 0.0, 0.0, 0.0]
    upper = [2.15, 2.15, 2.45, 160.0, 160.0, 750.0, 40.0, 40.0, 60.0, 0.4]

    def __init__(self, nu, coarse=(25, 5), fine=(256, 24)):
        super().__init__(nu, 1.5, CU_ISOTOPES, coarse, fine)


class Nitroxide(Centre):
    lower = [2.0040, 2.0030, 1.9990, 5.0, 5.0, 80.0, 0.8, 0.0, 0.0, 0.0]
    upper = [2.0130, 2.0100, 2.0060, 30.0, 30.0, 120.0, 8.0, 6.0, 6.0, 0.2]

    def __init__(self, nu, coarse=(21, 5), fine=(160, 24)):
        super().__init__(nu, 1.0, N14_ISOTOPES, coarse, fine)


def render(field, B, w, s, n_groups=36):
    """Sum of unit-area Gaussians (centre B, weight w, sigma s) on an even grid, by width groups:
    each group is histogrammed with linear sharing and convolved with its Gaussian by FFT."""
    step = field[1] - field[0]
    n = len(field)
    n_pad = 1 << int(np.ceil(np.log2(2 * n)))
    k = 2 * np.pi * np.fft.rfftfreq(n_pad, d=step)
    pos = (B - field[0]) / step
    ok = (pos > 0) & (pos < n - 2)
    pos, w, s = pos[ok], w[ok], s[ok]
    if s.max() / s.min() < 1.02:
        edges = np.array([s.min() - 1, s.max() + 1])
    else:
        edges = np.quantile(s, np.linspace(0, 1, n_groups + 1))
        edges[-1] += 1e-9
    out = np.zeros(n_pad)
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (s >= a) & (s < b)
        if not np.any(sel):
            continue
        i0 = np.floor(pos[sel]).astype(int)
        fr = pos[sel] - i0
        h = np.bincount(i0, w[sel] * (1 - fr), n_pad) + np.bincount(i0 + 1, w[sel] * fr, n_pad)
        sig = np.sqrt(np.average(s[sel] ** 2, weights=w[sel] + 1e-300))
        out += np.fft.irfft(np.fft.rfft(h) * np.exp(-0.5 * (k * sig) ** 2), n_pad)
    return out[:n] / step


def absorption(model, p, field):
    """Unit-area absorption on an evenly spaced field grid."""
    out = render(field, *model.lines(p))
    area = out.sum() * (field[1] - field[0])
    return out / area if area > 0 else out


def modulation_kernel(n_pad, step, bm, m_values, weights):
    k = 2 * np.pi * np.fft.rfftfreq(n_pad, d=step)
    ker = np.zeros(len(k), complex)
    for mv, wv in zip(m_values, weights):
        ker += wv * 2j * j1(k * bm * mv / 2)
    return ker / np.sum(weights)


def first_harmonic(absorp, kernel, n_pad):
    return np.fft.irfft(np.fft.rfft(absorp, n_pad) * kernel, n_pad)[:len(absorp)]
