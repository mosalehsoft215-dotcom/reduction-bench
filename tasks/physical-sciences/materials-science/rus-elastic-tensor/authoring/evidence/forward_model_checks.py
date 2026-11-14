"""Sanity checks on the forward model used by the data generator.

Run from the bundle root:  python authoring/evidence/forward_model_checks.py
Results from the last run are kept in forward_model_checks.txt next to this file.
"""

import os
import sys

import numpy as np
from scipy.linalg import eigh

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "provenance"))
import rusmodel as rm  # noqa: E402


def monomial_spectrum(dims_m, cv_pa, rho, order, n_modes):
    """Plain Visscher monomial basis, no symmetry blocking. Written separately
    from rusmodel on purpose, as an independent cross-check."""
    a, b, c = np.asarray(dims_m) / 2.0
    lmn = [(l, m, n) for l in range(order + 1) for m in range(order + 1 - l)
           for n in range(order + 1 - l - m)]
    nb = len(lmn)

    def f(p):
        return 0.0 if (p < 0 or p % 2) else 2.0 / (p + 1)

    half = (a, b, c)
    vol = a * b * c
    mass = np.zeros((nb, nb))
    dd = np.zeros((3, 3, nb, nb))
    for s, (l1, m1, n1) in enumerate(lmn):
        for t, (l2, m2, n2) in enumerate(lmn):
            p1 = (l1, m1, n1)
            p2 = (l2, m2, n2)
            mass[s, t] = vol * f(l1 + l2) * f(m1 + m2) * f(n1 + n2)
            for j in range(3):
                for k in range(3):
                    if p1[j] == 0 or p2[k] == 0:
                        continue
                    e = [p1[x] + p2[x] for x in range(3)]
                    e[j] -= 1
                    e[k] -= 1
                    dd[j, k, s, t] = (p1[j] * p2[k] * vol / (half[j] * half[k])
                                      * f(e[0]) * f(e[1]) * f(e[2]))
    c4 = rm.voigt_to_tensor(cv_pa)
    gam = np.zeros((3 * nb, 3 * nb))
    emat = np.zeros((3 * nb, 3 * nb))
    for i in range(3):
        emat[i * nb:(i + 1) * nb, i * nb:(i + 1) * nb] = rho * mass
        for k in range(3):
            blk = np.zeros((nb, nb))
            for j in range(3):
                for l in range(3):
                    blk += c4[i, j, k, l] * dd[j, l]
            gam[i * nb:(i + 1) * nb, k * nb:(k + 1) * nb] = blk
    w2 = eigh(gam, emat, eigvals_only=True)
    fr = np.sqrt(np.clip(w2, 0, None)) / (2 * np.pi)
    return np.sort(fr)[6:6 + n_modes]


def main():
    out = []
    say = out.append

    # --- a tetragonal crystal on a generic box ---------------------------
    dims = np.array([1.53e-3, 1.21e-3, 0.87e-3])
    rho = 7900.0
    cv = rm.tetragonal_voigt(171e9, 140e9, 52e9, 45e9, 64e9, 83e9)

    say("1) Legendre-block model vs independent monomial model, N = 8")
    ref = rm.Parallelepiped(dims, 8).solve(cv, rho, n_modes=40)["freq"]
    mono = monomial_spectrum(dims, cv, rho, 8, 40)
    rel = np.abs(ref / mono - 1)
    say(f"   max |rel diff| over 40 modes = {rel.max():.2e}")

    say("2) convergence of the lowest 60 modes against N = 28")
    top = rm.Parallelepiped(dims, 28).solve(cv, rho, n_modes=60)["freq"]
    for order in (8, 10, 12, 14, 16, 18, 20, 22, 24):
        fr = rm.Parallelepiped(dims, order).solve(cv, rho, n_modes=60)["freq"]
        r = fr / top - 1
        say(f"   N={order:2d}: max rel error modes 1-30 = {np.abs(r[:30]).max():.2e}, "
            f"modes 31-60 = {np.abs(r[30:]).max():.2e}")

    say("3) slender isotropic bar (10 x 1 x 1 mm) against beam theory")
    e_mod, nu, rho_s = 200e9, 0.3, 7800.0
    mu = e_mod / (2 * (1 + nu))
    lam = e_mod * nu / ((1 + nu) * (1 - 2 * nu))
    iso = rm.tetragonal_voigt(lam + 2 * mu, lam + 2 * mu, mu, mu, lam, lam)
    length, side = 10e-3, 1e-3
    bar = rm.Parallelepiped([side, side, length], 18).solve(iso, rho_s, vectors=True)
    fr = bar["freq"]
    grp = bar["group"]
    c_bar = np.sqrt(e_mod / rho_s)
    f_long = c_bar / (2 * length)
    f_tors = np.sqrt(mu / rho_s) / (2 * length) * np.sqrt(0.140577 * 6)
    f_flex = (4.730041 ** 2 / (2 * np.pi)) * np.sqrt(e_mod * side ** 2 / 12
                                                     / rho_s) / length ** 2
    # longitudinal fundamental: u_z odd in z, even in x and y -> group (0,0,0)
    long_modes = [f for f, g in zip(fr, grp) if g == (0, 0, 0)]
    # torsion about z: u_x ~ -y z-odd?  fundamental torsion is twist odd in z
    # -> u_x ~ y*z: (l,m,n)=(0,1,1) for i=0 -> group (1,1,1)
    tors_modes = [f for f, g in zip(fr, grp) if g == (1, 1, 1)]
    flex_modes = sorted(fr)[:2]
    say(f"   longitudinal: RR {long_modes[0]:.1f} Hz vs c/2L {f_long:.1f} Hz "
        f"(ratio {long_modes[0] / f_long:.5f}; Love correction ~ -0.07%)")
    say(f"   torsion:      RR {tors_modes[0]:.1f} Hz vs Saint-Venant {f_tors:.1f} Hz "
        f"(ratio {tors_modes[0] / f_tors:.5f})")
    say(f"   1st flexure:  RR {flex_modes[0]:.1f} Hz vs Euler-Bernoulli {f_flex:.1f} Hz "
        f"(ratio {flex_modes[0] / f_flex:.5f}; Timoshenko shear and rotary "
        "inertia take off about 3.5% at L/a = 10)")

    say("4) 90 degree rotation about z with Lx and Ly swapped must not change the spectrum")
    ortho = np.zeros((6, 6))
    ortho[0, 0], ortho[1, 1], ortho[2, 2] = 150e9, 190e9, 130e9
    ortho[3, 3], ortho[4, 4], ortho[5, 5] = 40e9, 55e9, 70e9
    ortho[0, 1] = ortho[1, 0] = 60e9
    ortho[0, 2] = ortho[2, 0] = 75e9
    ortho[1, 2] = ortho[2, 1] = 50e9
    rz = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    f1 = rm.Parallelepiped(dims, 12).solve(ortho, rho, n_modes=40)["freq"]
    f2 = rm.Parallelepiped(dims[[1, 0, 2]], 12).solve(
        rm.rotate_voigt(ortho, rz), rho, n_modes=40)["freq"]
    say(f"   max rel diff = {np.abs(f1 / f2 - 1).max():.2e}")
    rx = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])
    f3 = rm.Parallelepiped(dims[[0, 2, 1]], 12).solve(
        rm.rotate_voigt(ortho, rx), rho, n_modes=40)["freq"]
    say(f"   same test about x (Ly, Lz swapped): max rel diff = "
        f"{np.abs(f1 / f3 - 1).max():.2e}")

    say("5) 45 degree rotation of a 4/mmm tensor about c")
    s = 1 / np.sqrt(2)
    r45 = np.array([[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]])
    rot = rm.rotate_voigt(cv, r45) / 1e9
    c11, c12, c66 = 171.0, 64.0, 45.0
    say(f"   C11' {rot[0, 0]:.4f} (formula {(c11 + c12) / 2 + c66:.4f}), "
        f"C12' {rot[0, 1]:.4f} (formula {(c11 + c12) / 2 - c66:.4f}), "
        f"C66' {rot[5, 5]:.4f} (formula {(c11 - c12) / 2:.4f})")
    say(f"   largest off-pattern element after rotation: "
        f"{np.abs(rot[[0, 1, 2, 3, 4], 5]).max():.2e} GPa")

    text = "\n".join(out)
    print(text)
    with open(os.path.join(os.path.dirname(__file__), "forward_model_checks.txt"),
              "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
