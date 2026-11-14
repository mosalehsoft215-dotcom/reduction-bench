"""RBS forward model (4He on a layered target), shared by the generator and the reference.

Slab method. Each layer is cut into thin slabs. For every slab the incoming energy is tracked with
the tabulated stopping (Bragg's rule). For every isotope of every element in the slab (natural
abundances, each with its own kinematic factor), particles backscattered at the slab's front and
back faces (K*E) are carried out through the slabs above; the exit map E_in -> E_out and its Bohr
straggling variance are built outward slab by slab on an energy grid. The slab yield
Q*Omega*sigma(E)*N_i*dx/cos(alpha) is spread uniformly between the two exit energies and
broadened by a Gaussian combining detector resolution and straggling.
"Sn119" is tin from a 119Sn-enriched target (its own isotope list, the stopping of Sn). Hydrogen cannot
backscatter 4He (M2 < M1): it enters the stopping (Bragg's rule) and the straggling only.
"""
import numpy as np
from scipy.special import erf

M1, Z1 = 4.002602, 2
ELEM = {"O": (15.999, 8), "Al": (26.982, 13), "Si": (28.086, 14), "Ag": (107.868, 47), "Sn": (118.71, 50),
        "Sn119": (118.896, 50), "H": (1.008, 1)}
BASE = {"Sn119": "Sn"}            # isotopically enriched tin has the stopping of tin
# natural isotopic composition (mass u, abundance)
ISO = {"O": [(15.99491, 0.99757), (16.99913, 0.00038), (17.99916, 0.00205)],
       "Al": [(26.98154, 1.0)],
       "Si": [(27.97693, 0.92223), (28.97649, 0.04685), (29.97377, 0.03092)],
       "Ag": [(106.90509, 0.51839), (108.90476, 0.48161)],
       "Sn119": [(115.90174, 0.004), (116.90295, 0.005), (117.90161, 0.026), (118.90331, 0.926), (119.90220, 0.039)],
       "Sn": [(111.90482, 0.0097), (113.90278, 0.0066), (114.90334, 0.0034), (115.90174, 0.1454), (116.90295, 0.0768),
              (117.90161, 0.2422), (118.90331, 0.0859), (119.90220, 0.3258), (121.90344, 0.0463), (123.90527, 0.0579)]}
ISOTOPES = True                    # False: one kinematic factor per element at the mean atomic mass
E2 = 1.439964e-10                 # e^2/(4 pi eps0) in keV*cm
BOHR = 0.26 * Z1 ** 2 * 1e-3      # keV^2 per (Z2 * 1e15 at/cm2)


def kin(M2, theta):
    t = np.radians(theta)
    return ((np.sqrt(M2 ** 2 - (M1 * np.sin(t)) ** 2) + M1 * np.cos(t)) / (M1 + M2)) ** 2


def rutherford(E, M2, Z2, theta):
    """lab-frame Rutherford cross section (cm^2/sr) with Andersen screening; E in keV"""
    t = np.radians(theta)
    root = np.sqrt(M2 ** 2 - (M1 * np.sin(t)) ** 2)
    sig = (Z1 * Z2 * E2 / (4 * E)) ** 2 * 4 * (root + M2 * np.cos(t)) ** 2 / (M2 * np.sin(t) ** 4 * root)
    ecm = E * M2 / (M1 + M2)
    v1 = 0.04873 * Z1 * Z2 * np.sqrt(Z1 ** (2 / 3) + Z2 ** (2 / 3))
    tcm = t + np.arcsin(M1 / M2 * np.sin(t))
    F = (1 + 0.5 * v1 / ecm) ** 2 / (1 + v1 / ecm + (v1 / (2 * ecm * np.sin(tcm / 2))) ** 2) ** 2
    return sig * F


def isotopes(el):
    if ISOTOPES:
        m, a = zip(*ISO[el])
        return np.array(m), np.array(a)
    return np.array([ELEM[el][0]]), np.array([1.0])


def load_stopping(path):
    """file: first line header 'E_keV O Al Si Ag Sn'; values in eV per 1e15 atoms/cm2"""
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    head = lines[0].split()
    arr = np.loadtxt(lines[1:])
    return {el: (arr[:, 0], arr[:, i]) for i, el in enumerate(head) if i > 0}


def s_layer(table, comp, E):
    """Bragg's rule, keV per 1e15 atoms/cm2 (linear interpolation of the table)"""
    return sum(x * np.interp(E, *table[BASE.get(el, el)]) for el, x in comp.items()) * 1e-3


def _ierf(u):
    return u * erf(u) + np.exp(-u * u) / np.sqrt(np.pi)


def spectrum(layers, geom, table, channels, cal, fwhm_det, qomega, slab=2.0, sub_slab=10.0, e_floor=450.0):
    """layers: [(comp {el: atomic fraction}, areal density 1e15 at/cm2), ...] from the top; the last
    entry is the substrate (its density is ignored). geom: (E0 keV, theta, alpha, beta) degrees.
    cal: (offset keV, gain keV/channel). Returns expected counts per channel."""
    E0, theta, alpha, beta = geom
    ca, cb = np.cos(np.radians(alpha)), np.cos(np.radians(beta))
    off, gain = cal
    Ech = off + gain * np.asarray(channels, float)
    counts = np.zeros(Ech.size)
    grid = np.linspace(50.0, E0 * 1.001, 700)
    out_map = grid.copy()                 # exit energy for a particle leaving the current depth at energy grid
    out_var = np.zeros_like(grid)         # straggling variance accumulated on the way out
    iso = {}
    slabs = []
    for li, (comp, nt) in enumerate(layers):
        for el in comp:
            if el == "H":                 # no backscattering from hydrogen (M2 < M1): stopping and straggling only
                continue
            if el not in iso:
                m, a = isotopes(el)
                iso[el] = (kin(m, theta), m, a)
        if li == len(layers) - 1:
            slabs += [(comp, sub_slab)] * 100000
        else:
            n = max(1, int(np.ceil(nt / slab)))
            slabs += [(comp, nt / n)] * n
    E, var_in = E0, 0.0
    for comp, dx in slabs:
        zbar = sum(x * ELEM[el][1] for el, x in comp.items())
        dE = s_layer(table, comp, E) * dx / ca
        # exit map through this slab (built outward): leaving from the back face of this slab
        loss = s_layer(table, comp, grid) * dx / cb
        new_map = np.interp(grid - loss, grid, out_map, left=0.0)
        new_var = np.interp(grid - loss, grid, out_var, left=0.0) + BOHR * zbar * dx / cb
        Emid = E - 0.5 * dE
        v_in = var_in + BOHR * zbar * 0.5 * dx / ca
        done = True
        for el, x in comp.items():
            if el == "H":
                continue
            K, M2, ab = iso[el]
            Z2 = ELEM[el][1]
            ef = np.interp(K * E, grid, out_map)
            eb = np.interp(K * (E - dE), grid, new_map)
            vf = np.interp(K * E, grid, out_var)
            vb = np.interp(K * (E - dE), grid, new_var)
            if ef.max() > e_floor:
                done = False
            ok = eb >= 60
            if not ok.any():
                continue
            K, M2, ab, ef, eb, vf, vb = K[ok], M2[ok], ab[ok], ef[ok], eb[ok], vf[ok], vb[ok]
            y = qomega * rutherford(Emid, M2, Z2, theta) * x * ab * dx * 1e15 / ca
            s = np.sqrt((fwhm_det / 2.3548) ** 2 + K * K * v_in + 0.5 * (vf + vb))
            w = np.maximum(ef - eb, 1e-6)
            r2 = (s * np.sqrt(2.0))[:, None]
            a, b = (Ech + 0.5 * gain - eb[:, None]) / r2, (Ech - 0.5 * gain - eb[:, None]) / r2
            c, d = (Ech + 0.5 * gain - ef[:, None]) / r2, (Ech - 0.5 * gain - ef[:, None]) / r2
            counts += ((y / (2 * w))[:, None] * r2 * (_ierf(a) - _ierf(b) - _ierf(c) + _ierf(d))).sum(0)
        out_map, out_var = new_map, new_var
        var_in += BOHR * zbar * dx / ca
        E -= dE
        if done or E < 200:
            break
    return counts
