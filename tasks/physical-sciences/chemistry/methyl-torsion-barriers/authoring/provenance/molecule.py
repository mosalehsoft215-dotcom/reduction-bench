"""Geometry of the two planar conformers of 3,4-dimethylthiophene-2-carbaldehyde
and the inertial quantities of the two-top rotation-torsion Hamiltonian.

Atoms: S1 C2 C3 C4 C5 (ring), aldehyde C6 O7 H8 on C2, methyl C9 on C3
(H12 H13 H14), methyl C10 on C4 (H15 H16 H17), ring hydrogen H11 on C5.
Top 1 = the C3 methyl (next to the aldehyde), top 2 = the C4 methyl.
"""
import numpy as np

_H = 6.62607015e-34
_U = 1.66053906660e-27
BCONV = _H / (8.0 * np.pi ** 2 * _U * 1e-20) / 1e6  # 505379.0 MHz u A^2

MASS = {"H": 1.00782503207, "C": 12.0, "13C": 13.00335483778,
        "O": 15.99491461956, "S": 31.97207100, "34S": 33.96786690}


def _place(p_ref, p_ang, r, theta_deg, side):
    u = (p_ang - p_ref)
    u = u / np.linalg.norm(u)
    th = np.radians(theta_deg)
    c, s = np.cos(th), side * np.sin(th)
    v = np.array([c * u[0] - s * u[1], s * u[0] + c * u[1], 0.0])
    return p_ref + r * v


GEOM = {
    "cis": dict(
        r_S1C2=1.7295, r_C2C3=1.3835, r_C3C4=1.4330, r_C4C5=1.3625,
        a_S1C2C3=111.25, a_C2C3C4=111.65, a_C3C4C5=112.60,
        r_C2C6=1.4535, a_S1C2C6=120.35,
        r_C6O7=1.2165, a_C2C6O7=124.20, r_C6H8=1.1115, a_C2C6H8=114.70,
        r_C3C9=1.5010, a_C2C3C9=124.10,
        r_C4C10=1.5005, a_C3C4C10=125.30,
        r_C5H11=1.0775, a_S1C5H11=120.10,
        r_CH3=1.090560, a_HCH=108.20,
        dihedral_O=0.0,
    ),
    "trans": dict(
        r_S1C2=1.7275, r_C2C3=1.3855, r_C3C4=1.4315, r_C4C5=1.3630,
        a_S1C2C3=110.85, a_C2C3C4=112.15, a_C3C4C5=112.40,
        r_C2C6=1.4595, a_S1C2C6=117.55,
        r_C6O7=1.2145, a_C2C6O7=125.00, r_C6H8=1.1105, a_C2C6H8=114.20,
        r_C3C9=1.5040, a_C2C3C9=125.60,
        r_C4C10=1.5010, a_C3C4C10=124.80,
        r_C5H11=1.0775, a_S1C5H11=120.10,
        r_CH3=1.090560, a_HCH=108.20,
        dihedral_O=180.0,
    ),
}

METHYL_C = (8, 9)
RING_C = (2, 3)
METHYL_H = ((11, 12, 13), (14, 15, 16))
ATOM_NAMES = ["S1", "C2", "C3", "C4", "C5", "C6", "O7", "H8", "C9", "C10",
              "H11", "H12", "H13", "H14", "H15", "H16", "H17"]


def _methyl(Cr, Cm, r_ch, a_hch):
    axis = (Cm - Cr) / np.linalg.norm(Cm - Cr)
    cos_hch = np.cos(np.radians(a_hch))
    cth = np.sqrt((2.0 * cos_hch + 1.0) / 3.0)
    sth = np.sqrt(1.0 - cth ** 2)
    perp1 = np.array([-axis[1], axis[0], 0.0])
    perp2 = np.array([0.0, 0.0, 1.0])
    out = []
    for phi in (0.0, 120.0, 240.0):
        ph = np.radians(phi)
        d = cth * axis + sth * (np.cos(ph) * perp1 + np.sin(ph) * perp2)
        out.append(Cm + r_ch * d)
    return out


def build(conformer):
    g = GEOM[conformer]
    S1 = np.zeros(3)
    C2 = np.array([g["r_S1C2"], 0.0, 0.0])
    C3 = _place(C2, S1, g["r_C2C3"], g["a_S1C2C3"], +1)
    C4 = _place(C3, C2, g["r_C3C4"], g["a_C2C3C4"], +1)
    C5 = _place(C4, C3, g["r_C4C5"], g["a_C3C4C5"], +1)
    centre = (S1 + C2 + C3 + C4 + C5) / 5.0

    def outside(p_ref, p_ang, r, theta):
        a = _place(p_ref, p_ang, r, theta, +1)
        b = _place(p_ref, p_ang, r, theta, -1)
        return a if np.linalg.norm(a - centre) > np.linalg.norm(b - centre) else b

    C6 = outside(C2, S1, g["r_C2C6"], g["a_S1C2C6"])
    o_a = _place(C6, C2, g["r_C6O7"], g["a_C2C6O7"], +1)
    o_b = _place(C6, C2, g["r_C6O7"], g["a_C2C6O7"], -1)
    s_side = +1 if np.linalg.norm(o_a - S1) < np.linalg.norm(o_b - S1) else -1
    o_side = s_side if g["dihedral_O"] == 0.0 else -s_side
    O7 = _place(C6, C2, g["r_C6O7"], g["a_C2C6O7"], o_side)
    H8 = _place(C6, C2, g["r_C6H8"], g["a_C2C6H8"], -o_side)
    C9 = outside(C3, C2, g["r_C3C9"], g["a_C2C3C9"])
    C10 = outside(C4, C3, g["r_C4C10"], g["a_C3C4C10"])
    H11 = outside(C5, S1, g["r_C5H11"], g["a_S1C5H11"])
    h1 = _methyl(C3, C9, g["r_CH3"], g["a_HCH"])
    h2 = _methyl(C4, C10, g["r_CH3"], g["a_HCH"])
    xyz = np.array([S1, C2, C3, C4, C5, C6, O7, H8, C9, C10, H11] + h1 + h2)
    elem = ["S", "C", "C", "C", "C", "C", "O", "H", "C", "C", "H"] + ["H"] * 6
    return xyz, elem


def inertial(conformer, subst=None):
    """Principal-axis quantities; lam[i] is the unit vector of top i
    (ring C -> methyl C) in the principal frame (a, b, c)."""
    xyz, elem = build(conformer)
    masses = np.array([MASS[e] for e in elem])
    if subst:
        for i, lab in subst.items():
            masses[i] = MASS[lab]
    com = (masses[:, None] * xyz).sum(0) / masses.sum()
    x = xyz - com
    I = np.zeros((3, 3))
    for m, r in zip(masses, x):
        I += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    w, V = np.linalg.eigh(I)
    pa = x @ V
    if pa[0, 0] < 0:
        V[:, 0] *= -1
    pa = x @ V
    if pa[METHYL_C[0], 1] < 0:
        V[:, 1] *= -1
    V[:, 2] = np.cross(V[:, 0], V[:, 1])
    pa = x @ V
    A, B, C = BCONV / w
    lam, Ia = [], []
    for t in range(2):
        axis = xyz[METHYL_C[t]] - xyz[RING_C[t]]
        axis = axis / np.linalg.norm(axis)
        lam.append(V.T @ axis)
        Ial = 0.0
        for i in METHYL_H[t]:
            d = xyz[i] - xyz[METHYL_C[t]]
            dp = d - np.dot(d, axis) * axis
            Ial += masses[i] * np.dot(dp, dp)
        Ia.append(Ial)
    lam = np.array(lam)
    Ia = np.array(Ia)
    rho = lam * Ia[:, None] / w[None, :]
    r = np.diag(Ia) - np.einsum("ig,jg,g->ij", lam, lam, 1.0 / w) * np.outer(Ia, Ia)
    Fm = BCONV * np.linalg.inv(r)
    ang = np.degrees(np.arccos(np.abs(lam[:, 0])))
    return dict(A=A, B=B, C=C, I=w, coords=pa, lam=lam, Ialpha=Ia, rho=rho, r=r,
                F=Fm, angle=ang, masses=masses, Pcc=0.5 * (w[0] + w[1] - w[2]))


def dipole(conformer, subst=None):
    xyz, elem = build(conformer)
    info = inertial(conformer, subst)
    masses = info["masses"]
    com = (masses[:, None] * xyz).sum(0) / masses.sum()
    x = xyz - com
    V = np.linalg.lstsq(x, info["coords"], rcond=None)[0]

    def unit(a, b):
        d = x[b] - x[a]
        return d / np.linalg.norm(d)
    ring_centre = x[[0, 1, 2, 3, 4]].mean(0)
    u_ring = ring_centre - x[0]
    u_ring /= np.linalg.norm(u_ring)
    mu = (2.45 * unit(5, 6) + 0.35 * unit(1, 5) + 0.52 * u_ring
          + 0.30 * unit(2, 8) + 0.30 * unit(3, 9))
    return mu @ V


if __name__ == "__main__":
    for conf in ("cis", "trans"):
        d = inertial(conf)
        print(conf, "A B C = %.3f %.3f %.3f" % (d["A"], d["B"], d["C"]),
              "kappa = %.4f" % ((2 * d["B"] - d["A"] - d["C"]) / (d["A"] - d["C"])),
              "Pcc = %.4f" % d["Pcc"])
        print("  lam =", np.round(d["lam"], 5).tolist(), "angles", np.round(d["angle"], 3))
        print("  Ialpha", np.round(d["Ialpha"], 5), "rho", np.round(d["rho"], 6).tolist())
        print("  F (MHz)", np.round(d["F"], 2).tolist(), " F12/F11 = %.5f" % (d["F"][0, 1] / d["F"][0, 0]))
        print("  mu", np.round(dipole(conf), 3))
        xyz, _ = build(conf)
        print("  S1-C5 %.4f  C9-C10 %.4f  H(Me1)-H(Me2) min %.3f" % (
            np.linalg.norm(xyz[4] - xyz[0]), np.linalg.norm(xyz[8] - xyz[9]),
            min(np.linalg.norm(xyz[i] - xyz[j]) for i in (11, 12, 13) for j in (14, 15, 16))))
