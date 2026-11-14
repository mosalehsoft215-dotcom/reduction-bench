"""Rotation-torsion Hamiltonian of a planar molecule with two C3v methyl tops
(principal axis method, PAM) and its rotational line strengths.

H = A Pa^2 + B Pb^2 + C Pc^2 + H_dist(Watson A, I^r)
    + sum_ij F_ij (p_i - rho_i . P)(p_j - rho_j . P)
    + V31/2 (1 - cos 3 a1) + V32/2 (1 - cos 3 a2) + V12 cos 3 a1 cos 3 a2

The kinetic part follows from the inertia of the whole molecule:
r_ij = d_ij Ia_i - Ia_i Ia_j sum_g lam_ig lam_jg / I_g,  F = (h / 8 pi^2) r^-1,
rho_ig = lam_ig Ia_i / I_g, with lam_i the unit vector of top i in the
principal frame. F_12 is the kinetic coupling of the two tops.

Molecule-fixed axes a=z, b=x, c=y; both top axes lie in the ab plane, so all
matrices are real. Torsional basis exp(i m1 a1) exp(i m2 a2), m_i = 3 k_i + s_i,
|k_i| <= kmax. Symmetry species (s1, s2) = (0,0), (0,1), (1,0), (1,1), (1,2).
Energies in MHz.

Two routes to the same eigenvalues:
  full_matrix - the whole |J K> x |m1> x |m2> product basis (validation only)
  Solver      - the torsional problem is diagonalised at every K first; the
                rotation-torsion matrix is then set up in the |J K> x |t_K>
                basis of the lowest nkeep torsional states of each K
"""
import math
from functools import lru_cache

import numpy as np

SPECIES = ((0, 0), (0, 1), (1, 0), (1, 1), (1, 2))


@lru_cache(maxsize=None)
def rot_ops(J):
    """Px, Py^2, Pz for one J (molecule-fixed; P+ lowers K)."""
    n = 2 * J + 1
    K = np.arange(-J, J + 1)
    JJ = J * (J + 1)
    Pp = np.zeros((n, n))
    Pm = np.zeros((n, n))
    for i, k in enumerate(K):
        if k - 1 >= -J:
            Pp[i - 1, i] = math.sqrt(JJ - k * (k - 1))
        if k + 1 <= J:
            Pm[i + 1, i] = math.sqrt(JJ - k * (k + 1))
    Px = 0.5 * (Pp + Pm)
    Py2 = -(Pp - Pm) @ (Pp - Pm) / 4.0
    Pz = np.diag(K.astype(float))
    return Px, Py2, Pz


def rigid_block(J, p):
    """Rigid rotor + quartic Watson A reduction (I^r)."""
    Px, Py2, Pz = rot_ops(J)
    Px2 = Px @ Px
    Pz2 = Pz @ Pz
    JJ = J * (J + 1.0)
    n = 2 * J + 1
    H = p["A"] * Pz2 + p["B"] * Px2 + p["C"] * Py2
    Pm2 = Px2 - Py2
    H = H - p["DJ"] * JJ ** 2 * np.eye(n) - p["DJK"] * JJ * Pz2 - p["DK"] * Pz2 @ Pz2 \
        - 2.0 * p["dJ"] * JJ * Pm2 - p["dK"] * (Pz2 @ Pm2 + Pm2 @ Pz2)
    return H


def tor_grid(s1, s2, kmax):
    m1 = 3.0 * np.arange(-kmax, kmax + 1) + s1
    m2 = 3.0 * np.arange(-kmax, kmax + 1) + s2
    return np.repeat(m1, m2.size), np.tile(m2, m1.size), m1.size


def potential(p, n):
    I = np.eye(n)
    c = 0.5 * (np.eye(n, k=1) + np.eye(n, k=-1))          # cos 3a
    return 0.5 * p["V31"] * np.kron(I - c, I) + 0.5 * p["V32"] * np.kron(I, I - c) \
        + p["V12"] * np.kron(c, c)


def full_matrix(J, s1, s2, p, kmax):
    """Whole product basis; index = iK * nt + (i1 * n + i2)."""
    Px, Py2, Pz = rot_ops(J)
    nr = 2 * J + 1
    M1, M2, n = tor_grid(s1, s2, kmax)
    nt = M1.size
    F, rho = p["F"], p["rho"]
    It, Ir = np.eye(nt), np.eye(nr)
    D = []
    for i, M in enumerate((M1, M2)):
        Q = rho[i, 0] * Pz + rho[i, 1] * Px
        D.append(np.kron(Ir, np.diag(M)) - np.kron(Q, It))
    H = np.kron(rigid_block(J, p), It) + np.kron(Ir, potential(p, n))
    for i in range(2):
        for j in range(2):
            H += F[i, j] * (D[i] @ D[j])
    return H


class Solver:
    """Two-step diagonalisation for one symmetry species (s1, s2).

    Stage 1: for every K, h_K = sum_ij F_ij (m_i - rho_ia K)(m_j - rho_ja K) + V
    is diagonalised; the lowest nkeep states are kept, and their order is made
    consistent from K to K by maximum overlap, so that column t is the same
    torsional state at every K.
    Stage 2: H(J) in the |J K> x |t_K> basis, with the b-axis terms
    -sum_i R_i (D_i Pb + Pb D_i) + (rho_b^T F rho_b) Pb^2 (R_i = sum_j F_ij rho_jb,
    D_i = p_i - rho_ia Pa) coupling K to K +- 1, K +- 2.
    """

    def __init__(self, p, s1, s2, Jmax, kmax=8, nkeep=50):
        self.p, self.Jmax, self.nkeep = p, Jmax, nkeep
        M1, M2, n = tor_grid(s1, s2, kmax)
        F, rho = p["F"], p["rho"]
        V = potential(p, n)
        self.U, self.e = {}, {}

        def solve(K):
            d1 = M1 - rho[0, 0] * K
            d2 = M2 - rho[1, 0] * K
            h = V + np.diag(F[0, 0] * d1 * d1 + F[1, 1] * d2 * d2 + 2.0 * F[0, 1] * d1 * d2)
            w, U = np.linalg.eigh(h)
            return w[:nkeep], U[:, :nkeep]

        self.e[0], self.U[0] = solve(0)
        for step in (1, -1):
            for K in range(step, step * (Jmax + 1), step):
                w, U = solve(K)
                prev = self.U[K - step]
                ov = np.abs(prev.T @ U)
                order = np.full(nkeep, -1)
                free = set(range(nkeep))
                for t in np.argsort(-ov.max(axis=1)):
                    c = max(free, key=lambda j: ov[t, j])
                    order[t] = c
                    free.discard(c)
                U = U[:, order]
                # a consistent sign for every state
                sg = np.sign(np.sum(prev * U, axis=0))
                sg[sg == 0] = 1.0
                self.U[K], self.e[K] = U * sg, w[order]
        self.Rb = F @ rho[:, 1]
        self.Xb2 = float(rho[:, 1] @ F @ rho[:, 1])
        self._S, self._D = {}, {}
        for K in range(-Jmax, Jmax + 1):
            for dK in (1, 2):
                Kp = K + dK
                if Kp > Jmax:
                    continue
                self._S[(Kp, K)] = self.U[Kp].T @ self.U[K]
                if dK == 1:
                    op = self.Rb[0] * (2.0 * M1 - rho[0, 0] * (K + Kp)) + \
                        self.Rb[1] * (2.0 * M2 - rho[1, 0] * (K + Kp))
                    self._D[(Kp, K)] = self.U[Kp].T @ (op[:, None] * self.U[K])

    def overlap(self, Kp, K):
        if Kp == K:
            return np.eye(self.nkeep)
        if Kp > K:
            return self._S[(Kp, K)]
        return self._S[(K, Kp)].T

    def matrix(self, J):
        Px, Py2, Pz = rot_ops(J)
        n = 2 * J + 1
        nk = self.nkeep
        O = rigid_block(J, self.p) + self.Xb2 * (Px @ Px)
        H = np.zeros((n * nk, n * nk))
        for a in range(n):
            K = a - J
            sl = slice(a * nk, (a + 1) * nk)
            H[sl, sl] = np.diag(self.e[K] + O[a, a])
            for dK in (1, 2):
                b = a + dK
                if b >= n:
                    continue
                slp = slice(b * nk, (b + 1) * nk)
                blk = O[b, a] * self._S[(K + dK, K)]
                if dK == 1:
                    blk = blk - Px[b, a] * self._D[(K + 1, K)]
                H[slp, sl] = blk
                H[sl, slp] = blk.T
        return H

    def levels(self, J, nlev):
        """Lowest nlev levels: energies, eigenvectors, dominant torsional state."""
        w, U = np.linalg.eigh(self.matrix(J))
        w, U = w[:nlev], U[:, :nlev]
        n = 2 * J + 1
        wt = (U.reshape(n, self.nkeep, nlev) ** 2).sum(axis=0)
        return w, U, np.argmax(wt, axis=0), wt.max(axis=0)


def rank_labels(J):
    """(Ka, Kc) of the energy ranks 0..2J of one torsional state."""
    out = []
    for r in range(2 * J + 1):
        tau = -J + r
        Ka = (tau + J + 1) // 2
        out.append((Ka, Ka - tau))
    return out


@lru_cache(maxsize=None)
def wigner3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0 or j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    f = math.factorial
    tri = f(j1 + j2 - j3) * f(j1 - j2 + j3) * f(-j1 + j2 + j3) / f(j1 + j2 + j3 + 1)
    pref = math.sqrt(tri * f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2)
                     * f(j3 + m3) * f(j3 - m3))
    s = 0.0
    for k in range(max(0, j2 - j3 - m1, j1 - j3 + m2), min(j1 + j2 - j3, j1 - m1, j2 + m2) + 1):
        s += (-1) ** k / (f(k) * f(j3 - j2 + k + m1) * f(j3 - j1 + k - m2)
                          * f(j1 + j2 - j3 - k) * f(j1 - k - m1) * f(j2 - k + m2))
    return (-1) ** (j1 - j2 - m3) * pref * s


@lru_cache(maxsize=None)
def dircos(Jp, J):
    """C_q[K', K] = sqrt((2J+1)(2J'+1)) (-1)^K' (J' 1 J; -K' q K), q = -1, 0, 1."""
    out = {}
    for q in (-1, 0, 1):
        M = np.zeros((2 * Jp + 1, 2 * J + 1))
        for i, Kp in enumerate(range(-Jp, Jp + 1)):
            for j, K in enumerate(range(-J, J + 1)):
                if Kp == K + q:
                    M[i, j] = math.sqrt((2 * J + 1) * (2 * Jp + 1)) * (-1) ** Kp \
                        * wigner3j(Jp, 1, J, -Kp, q, K)
        out[q] = M
    return out


def dipole_rot(Jp, J, mu_a, mu_b):
    C = dircos(Jp, J)
    s2 = math.sqrt(2.0)
    return mu_a * C[0] - (mu_b / s2) * C[1] + (mu_b / s2) * C[-1]


def strengths(solver, Jl, Ul, Ju, Uu, mu_a, mu_b):
    """mu^2 S (Debye^2) between the eigenvector sets of J = Jl and J = Ju."""
    O = dipole_rot(Ju, Jl, mu_a, mu_b)
    nk = solver.nkeep
    nl, nu = 2 * Jl + 1, 2 * Ju + 1
    T = np.zeros((nu * nk, nl * nk))
    for a in range(nl):
        K = a - Jl
        for dK in (-1, 0, 1):
            b = K + dK + Ju
            if 0 <= b < nu and O[b, a] != 0.0:
                T[b * nk:(b + 1) * nk, a * nk:(a + 1) * nk] = O[b, a] * solver.overlap(K + dK, K)
    return (Uu.T @ T @ Ul) ** 2
