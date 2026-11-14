"""Free vibrations of an anisotropic, optionally piezoelectric, parallelepiped.

Rayleigh-Ritz method of Visscher et al. (J. Acoust. Soc. Am. 90, 2154, 1991),
extended to piezoelectric solids as in Ohno (Phys. Chem. Minerals 17, 371,
1990): the electric potential is expanded in the same basis as the
displacements and eliminated, which adds a stiffening term
Omega Phi^-1 Omega^T to the elastic matrix. The field outside the sample is
neglected, so on bare faces the natural boundary condition is D.n = 0.
Optionally the two z faces carry electrodes. Grounded electrodes (shorted to
each other): the potential functions carry a factor (1 - (z/c)^2). Floating
electrodes (connected to nothing): each face is still an equipotential but
the two may differ, so one more potential function, z/c, joins the basis; its
coefficient is the voltage between the electrodes, and stationarity with
respect to it sets the net charge on each electrode to zero. That function is
orthogonal to the others in the dielectric energy (they vanish on both faces),
so it adds a rank-one term to the fully symmetric block only.

Electrodes of finite mass (silver paint, say) add a surface density on the z
faces to the kinetic energy; the mass matrix is then rho I + sigma S, with S
the face integrals of the basis. Their stiffness is neglected.

Basis: products of normalised Legendre polynomials of total degree <= N. They
span the same space as Visscher's monomials, so the eigenvalues are the same;
the mass matrix is simply the identity. With bare faces the constant potential
carries no field and is left out of the potential basis.

Symmetry: the elastic tensor must be at least orthorhombic with axes along the
edges. A piezoelectric tensor of point group 4mm (or 6mm, mm2) with its polar
axis along z is even under the x and y mirrors and odd under the z mirror.
The problem then still splits into 8 blocks (Ohno 1976); in block (gx, gy, gz)
the potential functions have parity (gx, gy, 1 - gz).

Units: lengths in m, stiffness in Pa, e in C/m^2, permittivity in F/m,
density in kg/m^3, frequencies in Hz.
"""

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.linalg import cho_factor, cho_solve, eigh

VOIGT_INDEX = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
EPS0 = 8.8541878128e-12


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


def piezo_to_tensor(ev):
    """3x6 Voigt piezoelectric matrix (e_k,alpha) -> e_kij."""
    e = np.empty((3, 3, 3))
    for k in range(3):
        for i in range(3):
            for j in range(3):
                e[k, i, j] = ev[k, VOIGT_INDEX[i, j]]
    return e


def rotate_piezo(ev, rot):
    e = np.einsum("ka,ib,jc,abc->kij", rot, rot, rot, piezo_to_tensor(ev))
    out = np.empty((3, 6))
    for k in range(3):
        for a, (i, j) in enumerate(VOIGT_PAIRS):
            out[k, a] = e[k, i, j]
    return out


def tetragonal_voigt(c11, c33, c44, c66, c12, c13):
    cv = np.zeros((6, 6))
    cv[0, 0] = cv[1, 1] = c11
    cv[2, 2] = c33
    cv[3, 3] = cv[4, 4] = c44
    cv[5, 5] = c66
    cv[0, 1] = cv[1, 0] = c12
    cv[0, 2] = cv[2, 0] = cv[1, 2] = cv[2, 1] = c13
    return cv


def d_4mm(d15, d31, d33):
    """Strain coefficients (C/N), IEEE 176 convention, polar axis along z."""
    d = np.zeros((3, 6))
    d[0, 4] = d[1, 3] = d15
    d[2, 0] = d[2, 1] = d31
    d[2, 2] = d33
    return d


def ieee_e_and_eps_s(d, eps_t, cv_e):
    """e = d c^E and eps^S = eps^T - d e^T (IEEE Std 176, Voigt)."""
    e = d @ cv_e
    return e, eps_t - d @ e.T


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


def _shorted_tables(order):
    """1D integrals for potential functions q_n(s) = p_n(s) (1 - s^2), which
    vanish on s = +-1 (grounded electrodes), against p_m and against q_m."""
    nodes, weights = npleg.leggauss(order + 4)
    p = np.empty((order + 1, nodes.size))
    dp = np.empty_like(p)
    for n in range(order + 1):
        coef = np.zeros(n + 1)
        coef[n] = 1.0
        norm = np.sqrt((2 * n + 1) / 2.0)
        p[n] = norm * npleg.legval(nodes, coef)
        dp[n] = norm * npleg.legval(nodes, npleg.legder(coef))
    q = p * (1 - nodes ** 2)
    dq = dp * (1 - nodes ** 2) - 2 * nodes * p
    fun = {("p", 0): p, ("p", 1): dp, ("q", 0): q, ("q", 1): dq}
    tabs = {}
    for ka in ("p", "q"):
        for kb in ("p", "q"):
            for da in (0, 1):
                for db in (0, 1):
                    tabs[(ka, kb, da, db)] = (fun[(ka, da)] * weights) @ fun[(kb, db)].T
    return tabs


class Parallelepiped:
    """Rayleigh-Ritz problem for one sample shape and polynomial order N.

    shorted_z=True puts grounded electrodes on the two z faces: the potential
    basis is multiplied by (1 - (z/c)^2) so it vanishes there, and the
    constant potential stays in (it is no longer constant). floating_z=True
    puts floating electrodes there instead (see the module docstring).
    """

    def __init__(self, dims_m, order, shorted_z=False, floating_z=False):
        self.dims = np.asarray(dims_m, dtype=float)
        self.half = self.dims / 2.0
        self.order = order
        shorted_z = shorted_z or floating_z
        self.shorted_z = shorted_z
        self.floating_z = floating_z
        self.tables, self.edge = _legendre_tables(order)
        self.ztables = _shorted_tables(order) if shorted_z else None
        lmn = [(l, m, n) for l in range(order + 1)
               for m in range(order + 1 - l)
               for n in range(order + 1 - l - m)]
        self.lmn = np.array(lmn)
        par = self.lmn % 2
        keep_const = shorted_z

        self.blocks = []
        for gx in range(2):
            for gy in range(2):
                for gz in range(2):
                    comps, lams = [], []
                    for i in range(3):
                        pi = (self.lmn + np.eye(3, dtype=int)[i]) % 2
                        sel = np.where((pi[:, 0] == gx) & (pi[:, 1] == gy)
                                       & (pi[:, 2] == gz))[0]
                        comps.extend([i] * sel.size)
                        lams.extend(sel.tolist())
                    psi = np.where((par[:, 0] == gx) & (par[:, 1] == gy)
                                   & (par[:, 2] == 1 - gz)
                                   & ((self.lmn.sum(axis=1) > 0) | keep_const))[0]
                    self.blocks.append(((gx, gy, gz), np.array(comps),
                                        np.array(lams), psi))

    def _deriv(self, a, b, j, k, kinds=("p", "p")):
        """int d_j(f_a) d_k(f_b) dV for index arrays a, b. kinds says whether
        each side is a displacement ("p") or potential ("q") function; they
        differ only along z, and only when the z faces are shorted."""
        d = np.ones((a.size, b.size))
        for axis in range(3):
            da, db = int(axis == j), int(axis == k)
            if axis == 2 and self.shorted_z and "q" in kinds:
                t = self.ztables[(kinds[0], kinds[1], da, db)]
            else:
                t = self.tables[(da, db)]
            d *= t[np.ix_(self.lmn[a, axis], self.lmn[b, axis])]
        return d / (self.half[j] * self.half[k])

    def _stiffness(self, c4, comps, lams):
        gam = np.zeros((comps.size, comps.size))
        for j in range(3):
            for l in range(3):
                gam += c4[:, j, :, l][np.ix_(comps, comps)] * self._deriv(lams, lams, j, l)
        return gam

    def _coupling(self, e3, comps, lams, psi):
        """Omega[(i,a), m] = sum_jk e_kij int d_j(phi_a) d_k(psi_m) dV."""
        om = np.zeros((comps.size, psi.size))
        for j in range(3):
            for k in range(3):
                coef = e3[k][:, j][comps]          # e_kij for each row's i
                if np.any(coef):
                    om += coef[:, None] * self._deriv(lams, psi, j, k, ("p", "q"))
        return om

    def _dielectric(self, eps, psi):
        ph = np.zeros((psi.size, psi.size))
        for i in range(3):
            for j in range(3):
                if eps[i, j] != 0.0:
                    ph += eps[i, j] * self._deriv(psi, psi, i, j, ("q", "q"))
        return ph

    def face_mass(self, comps, lams, components=(0, 1, 2)):
        """int over both z faces of f_a f_b dA for the same displacement
        component, in the units of the other matrices (divided by abc)."""
        l, m, n = self.lmn[lams].T
        same = ((comps[:, None] == comps[None, :])
                & (l[:, None] == l[None, :]) & (m[:, None] == m[None, :])
                & np.isin(comps, components)[:, None])
        # p_n(+1) p_n'(+1) + p_n(-1) p_n'(-1)
        both = (1.0 + (-1.0) ** (n[:, None] + n[None, :]))
        return same * np.outer(self.edge[n], self.edge[n]) * both / self.half[2]

    def floating_vector(self, e3, comps, lams):
        """Coupling of each displacement function to the potential z/c:
        (1/c) sum_j e_3ij int d_j f_a dV, in the units of the other matrices."""
        l, m, n = self.lmn[lams].T
        mean = lambda k: np.where(k == 0, np.sqrt(2.0), 0.0)         # int p_k
        jump = lambda k: np.where(k % 2 == 1, 2.0 * self.edge[k], 0.0)  # p_k(1) - p_k(-1)
        grads = (jump(l) * mean(m) * mean(n) / self.half[0],
                 mean(l) * jump(m) * mean(n) / self.half[1],
                 mean(l) * mean(m) * jump(n) / self.half[2])
        om = np.zeros(lams.size)
        for j in range(3):
            om += e3[2][comps, j] * grads[j]
        return om / self.half[2]

    def solve(self, cv_pa, density, n_modes=None, vectors=False,
              e_voigt=None, eps_s=None, surface_mass=0.0,
              mass_components=(0, 1, 2)):
        """Sorted elastic modes (the six rigid-body modes are dropped).

        With e_voigt (3x6, C/m^2) and eps_s (3x3, F/m) the piezoelectric
        coupling is included; cv_pa is then the constant-field stiffness.
        surface_mass (kg/m^2) sits on each z face and moves with it in the
        displacement components listed in mass_components.
        """
        c4 = voigt_to_tensor(cv_pa)
        e3 = piezo_to_tensor(e_voigt) if e_voigt is not None else None
        freqs, groups, shapes = [], [], []
        for key, comps, lams, psi in self.blocks:
            gam = self._stiffness(c4, comps, lams)
            if e3 is not None and psi.size:
                om = self._coupling(e3, comps, lams, psi)
                if np.any(om):
                    ph = cho_factor(self._dielectric(eps_s, psi))
                    gam = gam + om @ cho_solve(ph, om.T)
                    gam = 0.5 * (gam + gam.T)
            if e3 is not None and self.floating_z:
                og = self.floating_vector(e3, comps, lams)
                if np.any(og):
                    gam = gam + np.outer(og, og) / (8.0 * eps_s[2, 2] / self.half[2] ** 2)
            mass = density * np.eye(lams.size)
            if surface_mass:
                mass = mass + surface_mass * self.face_mass(comps, lams, mass_components)
            if vectors:
                w2, vec = eigh(gam, mass)
            else:
                w2, vec = eigh(gam, mass, eigvals_only=True), None
            for idx, val in enumerate(w2):
                freqs.append(np.sqrt(max(val, 0.0)) / (2 * np.pi))
                groups.append(key)
                if vectors:
                    shapes.append((comps, lams, vec[:, idx]))
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
