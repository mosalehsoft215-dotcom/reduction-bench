"""Checks of the forward model (authoring/provenance/pam.py) and of the
reference Hamiltonian (solution/mwfit.py). Writes forward_model_checks.txt.

1. angular momentum algebra
2. the two-step solver against the whole product basis
3. convergence in the torsional basis (kmax) and in the kept torsional states (nkeep)
4. the two-dimensional torsional problem at J = 0 against an independent
   finite-difference solver on a grid in (a1, a2), with F12 and V12
5. time reversal: species (0,1)/(0,2), (1,1)/(2,2) and (1,2)/(2,1) are degenerate
6. exchanging the two tops maps (1,0) onto (0,1)
7. line strengths of the two-step route against the whole product basis
8. the reference Hamiltonian (explicit matrix elements, its own code) against
   the generator, both conformers, all five species

Run from the task root: python authoring/evidence/forward_model_checks.py
"""
import importlib.util
import math
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.sparse import diags, identity, kron as skron
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "authoring", "provenance"))
import generate_spectra as gen  # noqa: E402
import molecule  # noqa: E402
import pam  # noqa: E402

CM = gen.CM
OUT = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    OUT.append(line)


def params(conf):
    return gen.hamiltonian_params(conf, molecule.inertial(conf))


def ground_freqs(p, Jmax, kmax, nkeep):
    """Ground-manifold transition frequencies J -> J, J+1 of all five species."""
    out = []
    for s in pam.SPECIES:
        S = pam.Solver(p, s[0], s[1], Jmax, kmax, nkeep)
        E = {J: S.levels(J, 2 * J + 1)[0] for J in range(Jmax + 1)}
        for J in range(Jmax):
            for Jp in (J, J + 1):
                d = E[Jp][None, :] - E[J][:, None]
                out.append(d.ravel())
    return np.concatenate(out)


t0 = time.time()
say("Forward-model checks, two-top PAM (%s)" % time.strftime("%Y-%m-%d"))
say("")

# 1 ------------------------------------------------------------------
worst = 0.0
for J in range(0, 14):
    Px, Py2, Pz = pam.rot_ops(J)
    worst = max(worst, np.abs(Px @ Px + Py2 + Pz @ Pz - J * (J + 1) * np.eye(2 * J + 1)).max())
say("1. Px^2 + Py^2 + Pz^2 = J(J+1), J <= 13: max deviation %.1e" % worst)

# 2 ------------------------------------------------------------------
say("2. two-step solver (kmax 6, nkeep 50) against the whole product basis (kmax 6),")
say("   lowest two torsional manifolds, J = 0..5:")
for conf in ("cis", "trans"):
    p = params(conf)
    worst = 0.0
    for s in pam.SPECIES:
        S = pam.Solver(p, s[0], s[1], 5, 6, 50)
        for J in range(6):
            n = 2 * (2 * J + 1)
            w2 = S.levels(J, n)[0]
            w = np.linalg.eigvalsh(pam.full_matrix(J, s[0], s[1], p, 6))[:n]
            worst = max(worst, np.abs(w - w2).max())
    say("   %-5s all five species: max |difference| %.1e MHz" % (conf, worst))

# 3 ------------------------------------------------------------------
say("3. convergence of the ground-state transition frequencies (J <= 10, five species):")
for conf in ("cis", "trans"):
    p = params(conf)
    ref = ground_freqs(p, 10, 10, 60)
    row = []
    for k in (5, 6, 7, 8):
        row.append("kmax %d: %.1e" % (k, np.abs(ground_freqs(p, 10, k, 60) - ref).max()))
    say("   %-5s (nkeep 60) against kmax 10: %s MHz" % (conf, ", ".join(row)))
    ref = ground_freqs(p, 10, 8, 70)
    row = []
    for nk in (20, 24, 30, 40, 50):
        row.append("nkeep %d: %.1e" % (nk, np.abs(ground_freqs(p, 10, 8, nk) - ref).max()))
    say("   %-5s (kmax 8) against nkeep 70: %s MHz" % (conf, ", ".join(row)))
say("   the generator uses kmax %d, nkeep %d; the reference fits use kmax 6, nkeep 24 and a"
    % (gen.KMAX, gen.NKEEP))
say("   final pass with kmax 7, nkeep 30")


# 4 ------------------------------------------------------------------
def fd_levels(p, N, nev):
    """J = 0: -sum_ij F_ij d_i d_j + V on an N x N periodic grid, 8th-order
    central differences, all nine (s1, s2) combinations at once."""
    h = 2.0 * np.pi / N
    c1 = np.array([1 / 280, -4 / 105, 1 / 5, -4 / 5, 0, 4 / 5, -1 / 5, 4 / 105, -1 / 280]) / h
    c2 = np.array([-1 / 560, 8 / 315, -1 / 5, 8 / 5, -205 / 72, 8 / 5, -1 / 5, 8 / 315, -1 / 560]) / h ** 2

    def circ(c):
        off = np.arange(-4, 5)
        r, cidx, v = [], [], []
        for i in range(N):
            for o, val in zip(off, c):
                if val != 0.0:
                    r.append(i)
                    cidx.append((i + o) % N)
                    v.append(val)
        from scipy.sparse import csr_matrix
        return csr_matrix((v, (r, cidx)), shape=(N, N))
    D1, D2 = circ(c1), circ(c2)
    I = identity(N, format="csr")
    F = p["F"]
    a = np.arange(N) * h
    A1, A2 = np.meshgrid(a, a, indexing="ij")
    V = 0.5 * p["V31"] * (1 - np.cos(3 * A1)) + 0.5 * p["V32"] * (1 - np.cos(3 * A2))         + p["V12"] * np.cos(3 * A1) * np.cos(3 * A2)
    # p_i = -i d_i  ->  p_i p_j = -d_i d_j
    H = (-F[0, 0] * skron(D2, I) - F[1, 1] * skron(I, D2) - 2.0 * F[0, 1] * skron(D1, D1)
         + diags(V.ravel())).tocsc()
    w = eigsh(H, k=nev, sigma=V.min(), which="LM", return_eigenvectors=False)
    return np.sort(w)


say("4. J = 0 torsional levels against a finite-difference solver (8th order, periodic grid):")
for conf in ("cis", "trans"):
    p = params(conf)
    spec = []
    for s1 in range(3):
        for s2 in range(3):
            S = pam.Solver(p, s1, s2, 0, 10, 12)
            spec.extend(S.e[0][:6])
    spec = np.sort(spec)[:24]
    for N in (96, 132):
        fd = fd_levels(p, N, 24)
        d = np.abs(fd - spec) / CM
        say("   %-5s N = %3d: lowest 24 levels (all nine species, 0 to %.1f cm-1) agree to %.1e cm-1 max"
            % (conf, N, (spec[-1] - spec[0]) / CM, d.max()))

# 5 ------------------------------------------------------------------
p = params("cis")
worst = 0.0
for a, b in (((0, 1), (0, 2)), ((1, 1), (2, 2)), ((1, 2), (2, 1)), ((1, 0), (2, 0))):
    Sa = pam.Solver(p, a[0], a[1], 8, 7, 30)
    Sb = pam.Solver(p, b[0], b[1], 8, 7, 30)
    for J in range(9):
        worst = max(worst, np.abs(Sa.levels(J, 2 * J + 1)[0] - Sb.levels(J, 2 * J + 1)[0]).max())
say("5. time reversal, cis, J <= 8: (0,1)=(0,2), (1,1)=(2,2), (1,2)=(2,1), (1,0)=(2,0) to %.1e MHz" % worst)

# 6 ------------------------------------------------------------------
q = dict(p)
q["rho"] = p["rho"][::-1].copy()
q["F"] = p["F"][::-1, ::-1].copy()
q["V31"], q["V32"] = p["V32"], p["V31"]
worst = 0.0
for a, b in (((1, 0), (0, 1)), ((1, 1), (1, 1)), ((1, 2), (2, 1))):
    Sa = pam.Solver(p, a[0], a[1], 8, 7, 30)
    Sb = pam.Solver(q, b[0], b[1], 8, 7, 30)
    for J in range(9):
        worst = max(worst, np.abs(Sa.levels(J, 2 * J + 1)[0] - Sb.levels(J, 2 * J + 1)[0]).max())
say("6. tops exchanged (rho, F, V3 swapped), cis, J <= 8: species map onto each other to %.1e MHz" % worst)

# 7 ------------------------------------------------------------------
worst, big = 0.0, 0.0
for conf in ("cis", "trans"):
    p = params(conf)
    mu = molecule.dipole(conf)
    for s in ((0, 0), (1, 1), (1, 2)):
        S = pam.Solver(p, s[0], s[1], 4, 6, 50)
        for J in range(3):
            for Jp in (J, J + 1):
                wl, Ul, _, _ = S.levels(J, 2 * J + 1)
                wu, Uu, _, _ = S.levels(Jp, 2 * Jp + 1)
                S2 = pam.strengths(S, J, Ul, Jp, Uu, mu[0], mu[1])
                # whole product basis: eigenvectors and the dipole operator directly
                Hl = pam.full_matrix(J, s[0], s[1], p, 6)
                Hu = pam.full_matrix(Jp, s[0], s[1], p, 6)
                el, vl = np.linalg.eigh(Hl)
                eu, vu = np.linalg.eigh(Hu)
                nt = (13) ** 2
                O = np.kron(pam.dipole_rot(Jp, J, mu[0], mu[1]), np.eye(nt))
                Sf = (vu[:, :2 * Jp + 1].T @ O @ vl[:, :2 * J + 1]) ** 2
                worst = max(worst, np.abs(S2 - Sf).max())
                big = max(big, Sf.max())
say("7. line strengths, two-step against whole product basis (J <= 3, species 00, 11, 12):")
say("   max |difference| %.1e D^2 (largest strength %.2f D^2)" % (worst, big))

# 8 ------------------------------------------------------------------
spec = importlib.util.spec_from_file_location("mwfit", os.path.join(ROOT, "solution", "mwfit.py"))
argv = sys.argv
sys.argv = [argv[0], "/nonexistent", "/nonexistent/out.json"]
mw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mw)
sys.argv = argv
say("8. the reference Hamiltonian (solution/mwfit.py) against the generator, J <= 10,")
say("   same A, B, C, distortion, V3, V12 and axis directions (kmax 8, nkeep 50):")
for conf in ("cis", "trans"):
    d = molecule.inertial(conf)
    th = [math.degrees(math.atan2(l[1], l[0])) for l in d["lam"]]
    dist = gen.DIST_KHZ[conf]
    theta = np.array([d["A"], d["B"], d["C"]] + [dist[k] * 1e-3 for k in ("DJ", "DJK", "DK", "dJ", "dK")]
                     + [gen.V31_CM[conf], gen.V32_CM[conf], gen.V12_CM[conf]] + th)
    F2, rho2 = mw.pam2_constants(theta)
    p = params(conf)
    say("   %-5s F and rho from the reference (I_alpha 3.146) against the generator's"
        " geometry (I_alpha %.6f): %.2e MHz, %.1e" % (conf, d["Ialpha"][0], np.abs(F2 - p["F"]).max(),
                                                   np.abs(rho2 - p["rho"]).max()))
    p = dict(p, F=F2, rho=rho2)
    worst = 0.0
    for s in pam.SPECIES:
        S = pam.Solver(p, s[0], s[1], 10, 8, 50)
        Er = mw.pam2_species(theta, s, 10, 8, 50)
        for J in range(11):
            worst = max(worst, np.abs(S.levels(J, 2 * J + 1)[0] - Er[J]).max())
    say("   %-5s all five species: max |E_reference - E_generator| %.1e MHz" % (conf, worst))

say("")
say("time %.0f s" % (time.time() - t0))
with open(os.path.join(HERE, "forward_model_checks.txt"), "w", newline="\n") as fh:
    fh.write("\n".join(OUT) + "\n")
