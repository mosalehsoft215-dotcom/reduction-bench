"""Steady-state voltammetry at a disk electrode set in an insulating plane.

One-electron couple O + e- = R with Butler-Volmer kinetics, optionally
followed by R <-> Z in solution (forward kc, back kb, Z inert, D_Z = D_R).

Oblate spheroidal coordinates (xi >= 0, 0 <= eta <= 1),
    r = a sqrt((1 + xi^2)(1 - eta^2)),   z = a xi eta,
put the disk (radius a) on xi = 0, the insulating plane on eta = 0 and the
axis on eta = 1. With xi = tan(theta), steady diffusion with a first-order
decay lam becomes
    d/dxi[(1 + xi^2) dc/dxi] + d/deta[(1 - eta^2) dc/deta] = (lam a^2 / D)(xi^2 + eta^2) c,
and with xi = tan(theta) the first term is a plain second derivative in
theta. Pure diffusion is linear in theta, so a uniform theta grid carries the
far field exactly and the bulk condition sits at theta = pi/2 (infinity).
Finite volumes on (theta, eta) keep the limiting current 4 F D C a exact on
any grid.

The pair R, Z is carried as S = R + Z (diffusion only) and W = kc R - kb Z
(diffusion and decay at kc + kb), as in cvsim.py. Each species' response at
the disk nodes to a flux imposed there is precomputed once (a Green's
matrix), so each potential costs one small dense solve.

Units: cm, s, mol/cm^3, V, A. Currents follow the IUPAC convention.
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

FARADAY = 96485.33212
GAS = 8.314462618


def _operator(theta, eta, kappa2):
    """Sparse finite-volume operator for the free nodes (theta_0 .. theta_{N-1}
    by eta_0 .. eta_M); the node row theta_N = pi/2 is the bulk.
    Returns the matrix and the coupling of each free node to the bulk row."""
    n_t, n_e = theta.size - 1, eta.size
    t_lo = np.concatenate([[0.0], 0.5 * (theta[:-2] + theta[1:-1])])
    t_hi = 0.5 * (theta[:-1] + theta[1:])
    x_lo, x_hi = np.tan(t_lo), np.tan(t_hi)
    dxi = x_hi - x_lo
    e_lo = np.concatenate([[0.0], 0.5 * (eta[:-1] + eta[1:])])
    e_hi = np.concatenate([0.5 * (eta[:-1] + eta[1:]), [1.0]])
    deta = e_hi - e_lo

    idx = np.arange(n_t * n_e).reshape(n_t, n_e)
    rows, cols, vals = [], [], []
    diag = np.zeros((n_t, n_e))
    # theta faces, between i and i+1 (i + 1 = n_t is the bulk row)
    tcoef = deta[None, :] / np.diff(theta)[:, None]          # (n_t, n_e)
    diag += tcoef
    diag[1:] += tcoef[:-1]
    rows.append(idx[:-1].ravel()); cols.append(idx[1:].ravel()); vals.append(-tcoef[:-1].ravel())
    rows.append(idx[1:].ravel()); cols.append(idx[:-1].ravel()); vals.append(-tcoef[:-1].ravel())
    to_bulk = np.zeros((n_t, n_e))
    to_bulk[-1] = tcoef[-1]
    # eta faces, between j and j+1
    emid = 0.5 * (eta[:-1] + eta[1:])
    ecoef = (1.0 - emid ** 2)[None, :] * dxi[:, None] / np.diff(eta)[None, :]
    diag[:, :-1] += ecoef
    diag[:, 1:] += ecoef
    rows.append(idx[:, :-1].ravel()); cols.append(idx[:, 1:].ravel()); vals.append(-ecoef.ravel())
    rows.append(idx[:, 1:].ravel()); cols.append(idx[:, :-1].ravel()); vals.append(-ecoef.ravel())
    # first-order decay, integrated over each control volume
    if kappa2 > 0:
        diag += kappa2 * ((x_hi ** 3 - x_lo ** 3)[:, None] / 3.0 * deta[None, :]
                          + dxi[:, None] * ((e_hi ** 3 - e_lo ** 3) / 3.0)[None, :])
    rows.append(idx.ravel()); cols.append(idx.ravel()); vals.append(diag.ravel())
    mat = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                     shape=(n_t * n_e, n_t * n_e)).tocsc()
    return mat, to_bulk.ravel()


class Disk:
    """Green's matrices of one disk geometry for given D values and decay."""

    def __init__(self, radius, *, n_theta=160, n_eta=96, stretch=2.0):
        self.a = radius
        self.theta = np.linspace(0.0, 0.5 * np.pi, n_theta + 1)
        # nodes crowd toward the disk edge (eta = 0), where the current is
        # largest when kinetics are fast
        s = np.linspace(0.0, 1.0, n_eta + 1)
        self.eta = s ** stretch
        e_lo = np.concatenate([[0.0], 0.5 * (self.eta[:-1] + self.eta[1:])])
        e_hi = np.concatenate([0.5 * (self.eta[:-1] + self.eta[1:]), [1.0]])
        self.w = 0.5 * (e_hi ** 2 - e_lo ** 2)       # int eta d(eta) per disk node
        self._cache = {}

    def green(self, kappa2):
        """Disk-node concentrations per unit (scaled) production at each disk
        node, for decay kappa2 = lam a^2 / D."""
        key = round(kappa2, 14)
        if key not in self._cache:
            mat, _ = _operator(self.theta, self.eta, kappa2)
            lu = splu(mat)
            n_e = self.eta.size
            rhs = np.zeros((mat.shape[0], n_e))
            rhs[np.arange(n_e), np.arange(n_e)] = 1.0     # disk nodes are row theta_0
            self._cache[key] = lu.solve(rhs)[:n_e, :]
        return self._cache[key]

    def current(self, e, *, e0, k0, alpha, d_ox, d_red, c_ox, c_red,
                kc=0.0, kb=0.0, temp=298.15):
        """Steady-state current (A) at potentials e (no iR, no charging)."""
        if kc > 0 and c_red > 0:
            raise ValueError("a follow-up step needs c_red = 0 in the bulk")
        f = FARADAY / (GAS * temp)
        a = self.a
        wa_ox = (a / d_ox) * self.w
        wa_red = (a / d_red) * self.w
        g_diff = self.green(0.0)          # O and S: diffusion only
        if kc > 0:
            g_w = self.green((kc + kb) * a * a / d_red)
            h = (kb * g_diff + kc * g_w) / (kc + kb)
        else:
            h = g_diff
        m_ox = g_diff * wa_ox[None, :]
        m_red = h * wa_red[None, :]
        eye = np.eye(self.eta.size)
        out = np.empty(np.size(e))
        for n, ee in enumerate(np.atleast_1d(np.asarray(e, float))):
            eta_e = f * (ee - e0)
            e1 = k0 * np.exp(-alpha * eta_e)
            e2 = k0 * np.exp((1 - alpha) * eta_e)
            scale = 1.0 + e1 + e2
            lhs = (eye + e1 * m_ox + e2 * m_red) / scale
            rhs = np.full(self.eta.size, (e1 * c_ox - e2 * c_red) / scale)
            j = np.linalg.solve(lhs, rhs)
            out[n] = -FARADAY * 2.0 * np.pi * a * a * float(self.w @ j)
        return out


def _mass(theta, eta):
    """int (xi^2 + eta^2) dxi deta over each control volume (free nodes)."""
    t_lo = np.concatenate([[0.0], 0.5 * (theta[:-2] + theta[1:-1])])
    t_hi = 0.5 * (theta[:-1] + theta[1:])
    x_lo, x_hi = np.tan(t_lo), np.tan(t_hi)
    e_lo = np.concatenate([[0.0], 0.5 * (eta[:-1] + eta[1:])])
    e_hi = np.concatenate([0.5 * (eta[:-1] + eta[1:]), [1.0]])
    return ((x_hi ** 3 - x_lo ** 3)[:, None] / 3.0 * (e_hi - e_lo)[None, :]
            + (x_hi - x_lo)[:, None] * ((e_hi ** 3 - e_lo ** 3) / 3.0)[None, :]).ravel()


class _Transient:
    """One species on the disk grid, backward Euler: in scaled form
        (a^2/D)(1/dt + lam) M c' + K c' = (a^2/(D dt)) M c + bulk + injection,
    with c' the next time level and the injection at the disk nodes."""

    def __init__(self, theta, eta, a, d, lam, dt, bulk):
        scale = a * a / d
        mat, to_bulk = _operator(theta, eta, scale * (1.0 / dt + lam))
        self.lu = splu(mat)
        self.mt = scale / dt * _mass(theta, eta)
        self.bulk_rhs = to_bulk * bulk
        n, n_e = mat.shape[0], eta.size
        unit = np.zeros((n, n_e))
        unit[np.arange(n_e), np.arange(n_e)] = 1.0
        self.h = self.lu.solve(unit)            # field response to a unit injection at each disk node
        self.g = self.h[:n_e, :]                # its disk-node block
        self.c = np.full(n, bulk, float)

    def base(self):
        return self.lu.solve(self.mt * self.c + self.bulk_rhs)


def transient_cv(e_prog, dt, *, radius, e0, k0, alpha, d_ox, d_red, c_ox, c_red,
                 kc=0.0, kb=0.0, cdl=0.0, temp=298.15, n_theta=160, n_eta=96, stretch=2.0):
    """Current (A) at a disk of radius `radius` for the applied potentials
    e_prog[n] at t_n = n dt, starting from the uniform bulk at t = 0.

    Butler-Volmer at every disk node, R <-> Z in solution as S = R + Z and
    W' = (kc R - kb Z)/(kc + kb) (so with a follow-up step the bulk must hold
    no R), and double-layer charging cdl dE/dt (cdl in F for the whole disk).
    No iR drop: at a microdisk the currents are nA.
    """
    if kc > 0 and c_red > 0:
        raise ValueError("a follow-up step needs c_red = 0 in the bulk")
    f = FARADAY / (GAS * temp)
    a = radius
    theta = np.linspace(0.0, 0.5 * np.pi, n_theta + 1)
    eta = np.linspace(0.0, 1.0, n_eta + 1) ** stretch
    e_lo = np.concatenate([[0.0], 0.5 * (eta[:-1] + eta[1:])])
    e_hi = np.concatenate([0.5 * (eta[:-1] + eta[1:]), [1.0]])
    w = 0.5 * (e_hi ** 2 - e_lo ** 2)
    n_e = eta.size
    ox = _Transient(theta, eta, a, d_ox, 0.0, dt, c_ox)
    s_sp = _Transient(theta, eta, a, d_red, 0.0, dt, c_red)
    chem = kc > 0
    if chem:
        w_sp = _Transient(theta, eta, a, d_red, kc + kb, dt, 0.0)
        th_b, th_f = kb / (kc + kb), kc / (kc + kb)
        h_red = th_b * s_sp.g + th_f * w_sp.g
    else:
        h_red = s_sp.g
    wa_ox = (a / d_ox) * w
    wa_red = (a / d_red) * w
    m_ox = ox.g * wa_ox[None, :]
    m_red = h_red * wa_red[None, :]
    eye = np.eye(n_e)
    out = np.zeros(len(e_prog))
    for n in range(1, len(e_prog)):
        b_o = ox.base()
        b_s = s_sp.base()
        if chem:
            b_w = w_sp.base()
            r0 = th_b * b_s[:n_e] + b_w[:n_e]
        else:
            r0 = b_s[:n_e]
        eta_e = f * (e_prog[n] - e0)
        e1 = k0 * np.exp(-alpha * eta_e)
        e2 = k0 * np.exp((1 - alpha) * eta_e)
        sc = 1.0 + e1 + e2
        j = np.linalg.solve((eye + e1 * m_ox + e2 * m_red) / sc, (e1 * b_o[:n_e] - e2 * r0) / sc)
        ox.c = b_o - ox.h @ (wa_ox * j)
        s_sp.c = b_s + s_sp.h @ (wa_red * j)
        if chem:
            w_sp.c = b_w + th_f * (w_sp.h @ (wa_red * j))
        out[n] = -FARADAY * 2.0 * np.pi * a * a * float(w @ j) + cdl * (e_prog[n] - e_prog[n - 1]) / dt
    return out


def aoki_peak(p):
    """Aoki et al. (J. Electroanal. Chem. 171 (1984) 219): peak current of a
    reversible CV at a disk over 4 n F D C r, with p = r sqrt(n F v / (R T D))."""
    return 0.34 * np.exp(-0.66 * p) + 0.66 - 0.13 * np.exp(-11.0 / p) + 0.351 * p


def oldham_disk(e, *, e0, k0, alpha, d, c_ox, radius, temp=298.15):
    """Oldham's approximation (J. Electroanal. Chem. 313 (1991) 3) for a
    quasi-reversible steady-state wave at a disk, equal D, reduction of O."""
    f = FARADAY / (GAS * temp)
    eta_e = f * (np.asarray(e, float) - e0)
    theta = 1.0 + np.exp(eta_e)
    kap = np.pi * k0 * radius / (4.0 * d) * np.exp(-alpha * eta_e)
    ilim = -4.0 * FARADAY * d * c_ox * radius
    return ilim / (theta + (np.pi / kap) * (2 * kap * theta + 3 * np.pi) / (4 * kap * theta + 3 * np.pi ** 2))
