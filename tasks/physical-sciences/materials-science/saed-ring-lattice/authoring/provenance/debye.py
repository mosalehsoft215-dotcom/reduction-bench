"""Kinematic electron scattering of one whole crystallite, orientation averaged (Debye scattering equation).

A spherical crystallite of diameter D is cut from the ideal lattice about a cation site.  Its orientation-averaged
intensity at scattering vector s = 1/d (1/A) is
    I(s) = sum_i f_i(s)^2 + 2 sum_{i<j} f_i(s) f_j(s) sin(2 pi s r_ij) / (2 pi s r_ij),
evaluated with pair-distance histograms per species pair.  f(s) = sum_k a_k exp(-b_k s^2) from
scattering_factors.json; mixed cation sites use the composition-weighted factor.
"""
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist, pdist

SITES = {
    "rocksalt": ([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
                 [[.5, 0, 0], [0, .5, 0], [0, 0, .5], [.5, .5, .5]]),
    "fluorite": ([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
                 [[x, y, z] for x in (.25, .75) for y in (.25, .75) for z in (.25, .75)]),
}


def factor(table, composition, s):
    """Composition-weighted electron scattering factor (A) at s (1/A)."""
    s2 = np.asarray(s) ** 2
    out = 0.0
    for element, share in composition.items():
        a, b = table[element]["a"], table[element]["b"]
        out = out + share * sum(ai * np.exp(-bi * s2) for ai, bi in zip(a, b))
    return out


def cluster(kind, a, diameter_A):
    cation, anion = (np.array(x, dtype=float) for x in SITES[kind])
    n = int(diameter_A / (2 * a)) + 2
    cells = np.array(list(itertools.product(range(-n, n + 1), repeat=3)), dtype=float)
    cat = ((cells[:, None, :] + cation[None]).reshape(-1, 3)) * a
    an = ((cells[:, None, :] + anion[None]).reshape(-1, 3)) * a
    r = diameter_A / 2
    return cat[np.linalg.norm(cat, axis=1) <= r + 1e-9], an[np.linalg.norm(an, axis=1) <= r + 1e-9]


def histograms(cat, an, width=0.002):
    out = []
    for kind, d in (("cc", pdist(cat)), ("aa", pdist(an)), ("ca", cdist(cat, an).ravel())):
        h, e = np.histogram(d, bins=np.arange(0.0, d.max() + 2 * width, width))
        keep = h > 0
        out.append((kind, 0.5 * (e[1:] + e[:-1])[keep], h[keep].astype(float)))
    return out


def intensity(s, cat, an, f_cat, f_an, hist=None):
    """Orientation-averaged intensity of the crystallite (f_cat, f_an evaluated at s)."""
    hist = hist or histograms(cat, an)
    total = len(cat) * f_cat ** 2 + len(an) * f_an ** 2
    weights = {"cc": 2 * f_cat * f_cat, "aa": 2 * f_an * f_an, "ca": 2 * f_cat * f_an}
    for kind, r, h in hist:
        x = 2 * np.pi * np.outer(s, r)
        total = total + weights[kind] * (np.sinc(x / np.pi) @ h)
    return total


def load_table(path):
    return json.loads(Path(path).read_text())
