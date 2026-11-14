"""Free-vibration spectrum of an anisotropic rectangular parallelepiped.

Rayleigh-Ritz method of Visscher et al. (J. Acoust. Soc. Am. 90, 2154, 1991)
with a product basis of normalised Legendre polynomials instead of plain
monomials. Both span the same space (total degree <= N), so the eigenvalues are
identical; the Legendre basis just keeps the mass matrix at identity and the
problem well conditioned at high N.

The elastic tensor must have at least orthorhombic symmetry with its axes along
the sample edges (true for a tetragonal 4/mmm crystal cut on {110} or {100}
faces). The problem then splits into 8 parity blocks (Ohno 1976), solved
independently.

Units: lengths in m, stiffness in Pa, density in kg/m^3, frequencies in Hz.
"""

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.linalg import eigh

VOIGT_INDEX = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def voigt_to_tensor(cv):
    c = np.empty((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    c[i, j, k, l] = cv[VOIGT_INDEX[i, j], VOIGT_INDEX[k, l]]
    return c


def tensor_to_voigt(c):
    cv = np.empty((6, 6))
    for p, (i, j) in enumerate(VOIGT_PAIRS):
        for q, (k, l) in enumerate(VOIGT_PAIRS):
            cv[p, q] = c[i, j, k, l]
    return cv


def rotate_voigt(cv, rot):
    """Rows of `rot` are the new axes written in the old frame."""
    c = voigt_to_tensor(cv)
    c = np.einsum("ia,jb,kc,ld,abcd->ijkl", rot, rot, rot, rot, c)
    return tensor_to_voigt(c)


def tetragonal_voigt(c11, c33, c44, c66, c12, c13):
    cv = np.zeros((6, 6))
    cv[0, 0] = cv[1, 1] = c11
    cv[2, 2] = c33
    cv[3, 3] = cv[4, 4] = c44
    cv[5, 5] = c66
    cv[0, 1] = cv[1, 0] = c12
    cv[0, 2] = cv[2, 0] = cv[1, 2] = cv[2, 1] = c13
    return cv


def _legendre_tables(order):
    nodes, weights = npleg.leggauss(order + 2)
    p = np.empty((order + 1, nodes.size))
    dp = np.empty_like(p)
    for n in range(order + 1):
        coef = np.zeros(n + 1)
        coef[n] = 1.0
        norm = np.sqrt((2 * n + 1) / 2.0)
        p[n] = norm * npleg.legval(nodes, coef)
        dp[n] = norm * npleg.legval(nodes, npleg.legder(coef))
    i00 = (p * weights) @ p.T
    i10 = (dp * weights) @ p.T      # derivative on the first function
    i11 = (dp * weights) @ dp.T
    edge = np.sqrt((2 * np.arange(order + 1) + 1) / 2.0)   # p_n(+1)
    return {(0, 0): i00, (1, 0): i10, (0, 1): i10.T, (1, 1): i11}, edge


class Parallelepiped:
    """Rayleigh-Ritz problem for one sample shape and polynomial order N."""

    def __init__(self, dims_m, order):
        self.dims = np.asarray(dims_m, dtype=float)
        self.half = self.dims / 2.0
        self.order = order
        self.tables, self.edge = _legendre_tables(order)
        lmn = [(l, m, n) for l in range(order + 1)
               for m in range(order + 1 - l)
               for n in range(order + 1 - l - m)]
        self.lmn = np.array(lmn)

        # parity blocks: displacement component i expanded in (l, m, n)
        self.blocks = []
        for gx in range(2):
            for gy in range(2):
                for gz in range(2):
                    comps, lams = [], []
                    for i in range(3):
                        par = (self.lmn + np.eye(3, dtype=int)[i]) % 2
                        sel = np.where((par[:, 0] == gx) & (par[:, 1] == gy)
                                       & (par[:, 2] == gz))[0]
                        comps.extend([i] * sel.size)
                        lams.extend(sel.tolist())
                    self.blocks.append(((gx, gy, gz), np.array(comps),
                                        np.array(lams)))

    def _stiffness(self, c4, comps, lams):
        """Gamma[(i,a),(k,b)] = sum_jl C_ijkl * int d_j(phi_a) d_l(phi_b) dV."""
        idx = [self.lmn[lams, axis] for axis in range(3)]
        gam = np.zeros((comps.size, comps.size))
        for j in range(3):
            for l in range(3):
                d = np.ones_like(gam)
                for axis in range(3):
                    t = self.tables[(int(axis == j), int(axis == l))]
                    d *= t[np.ix_(idx[axis], idx[axis])]
                d /= self.half[j] * self.half[l]
                gam += c4[:, j, :, l][np.ix_(comps, comps)] * d
        return gam

    def solve(self, cv_pa, density, n_modes=None, vectors=False):
        """Sorted elastic modes (the six rigid-body modes are dropped)."""
        c4 = voigt_to_tensor(cv_pa)
        freqs, groups, shapes = [], [], []
        for key, comps, lams in self.blocks:
            gam = self._stiffness(c4, comps, lams)
            if vectors:
                w2, vec = eigh(gam)
            else:
                w2, vec = eigh(gam, eigvals_only=True), None
            w2 = w2 / density
            for idx, val in enumerate(w2):
                freqs.append(np.sqrt(max(val, 0.0)) / (2 * np.pi))
                groups.append(key)
                if vectors:
                    shapes.append((comps, lams, vec[:, idx] / np.sqrt(density)))
        freqs = np.array(freqs)
        order = np.argsort(freqs)
        cutoff = 1e-3 * freqs[order[6]]
        keep = [k for k in order if freqs[k] > cutoff]
        if len(order) - len(keep) != 6:
            raise RuntimeError("expected exactly six rigid-body modes")
        if n_modes is not None:
            keep = keep[:n_modes]
        out = {"freq": freqs[keep], "group": [groups[k] for k in keep]}
        if vectors:
            out["shape"] = [shapes[k] for k in keep]
        return out

    def displacement(self, shape, corner_signs):
        """Displacement of a mass-normalised mode at the corner (+-a, +-b, +-c)."""
        comps, lams, coef = shape
        val = np.ones(lams.size)
        for axis in range(3):
            n = self.lmn[lams, axis]
            sign = 1.0 if corner_signs[axis] > 0 else -1.0
            val *= self.edge[n] * sign ** n
        val /= np.sqrt(np.prod(self.half))
        u = np.zeros(3)
        np.add.at(u, comps, coef * val)
        return u
