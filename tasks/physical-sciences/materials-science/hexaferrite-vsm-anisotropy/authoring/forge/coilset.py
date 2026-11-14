"""Response of the magnetometer's pick-up coils to an extended, uniformly magnetised body.

Lab axes: X along the applied field (the coil axes), Z vertical (the vibration direction), Y the third. Each pick-up
coil acts as one circular turn of its mean radius, axis along X, centred where instrument.json puts it, with a
winding sense of +1 or -1. By reciprocity the signal of a moment component m_X at r is m_X * G(r), where
G = dB_X/dZ is the field gradient the coil set would make per unit current; a uniformly magnetised body therefore
gives m_X times the volume average of G. The day's moment factor is fixed with the reference sphere at the
coil-set centre, so a body's reported moment carries the extra factor <G>_body / G(0).
"""
import math

import numpy as np
from scipy.special import ellipe, ellipk

MU0 = 4e-7 * math.pi


def loop_bz(a, rho, z):
    """Axial field of a unit-current circular turn (radius a, centre at the origin, axis z) at (rho, z)."""
    rho = np.maximum(np.asarray(rho, float), 1e-12)
    z = np.asarray(z, float)
    s = (a + rho) ** 2 + z ** 2
    k2 = 4.0 * a * rho / s
    near = (a - rho) ** 2 + z ** 2
    return MU0 / (2.0 * math.pi * np.sqrt(s)) * (ellipk(k2) + (a * a - rho * rho - z * z) / near * ellipe(k2))


def coil_bx(points_mm, coils):
    """B_X (tesla per ampere) of the coil set at lab points (n, 3) in mm (lengths in mm, scale irrelevant)."""
    p = np.atleast_2d(points_mm)
    out = np.zeros(len(p))
    for c in coils:
        q = p - np.asarray(c["centre_mm"], float)
        out += c["sense"] * loop_bz(c["radius_mm"], np.hypot(q[:, 1], q[:, 2]), q[:, 0])
    return out


def gradient(points_mm, coils, step=1e-3):
    p = np.atleast_2d(points_mm).astype(float)
    up, dn = p.copy(), p.copy()
    up[:, 2] += step
    dn[:, 2] -= step
    return (coil_bx(up, coils) - coil_bx(dn, coils)) / (2.0 * step)


def _gauss(n, half):
    x, w = np.polynomial.legendre.leggauss(n)
    return x * half, w * half


def body_points(shape, lab, substrate_mm=None, n=10):
    """Quadrature points (lab mm) and weights of the body. lab maps lab axis 'X','Y','Z' -> sample axis."""
    to_lab = {v: k for k, v in lab.items()}          # sample axis -> lab axis
    idx = {"X": 0, "Y": 1, "Z": 2}
    kind = shape["kind"]
    if kind == "sphere":
        return np.zeros((1, 3)), np.ones(1)          # mean value property: a ball averages G to G(centre)
    if kind in ("box", "film"):
        e = dict(shape["edges_mm"])
        offset = {}
        if kind == "film":
            e[shape["normal"]] = 0.0
            if substrate_mm:
                offset[shape["normal"]] = 0.5 * substrate_mm
        grids = []
        for ax in "xyz":
            if e[ax] > 0:
                x, w = _gauss(n, e[ax] / 2.0)
            else:
                x, w = np.zeros(1), np.ones(1)
            grids.append((ax, x + offset.get(ax, 0.0), w))
        P = np.zeros((n ** 3, 3))
        X = np.meshgrid(*[g[1] for g in grids], indexing="ij")
        W = np.einsum("i,j,k->ijk", *[g[2] for g in grids]).ravel()
        P = np.zeros((W.size, 3))
        for (ax, _, _), coords in zip(grids, X):
            P[:, idx[to_lab[ax]]] = coords.ravel()
        return P, W
    if kind == "disc":
        R, t, axis = shape["diameter_mm"] / 2.0, shape["thickness_mm"], shape["axis"]
        r, wr = _gauss(n, R / 2.0)
        r, wr = r + R / 2.0, wr * (r + R / 2.0)
        phi = np.linspace(0.0, 2 * math.pi, 4 * n, endpoint=False)
        z, wz = _gauss(n, t / 2.0)
        others = [ax for ax in "xyz" if ax != axis]
        Rr, Pp, Zz = np.meshgrid(r, phi, z, indexing="ij")
        W = np.einsum("i,j,k->ijk", wr, np.full(phi.size, 2 * math.pi / phi.size), wz).ravel()
        P = np.zeros((W.size, 3))
        P[:, idx[to_lab[axis]]] = Zz.ravel()
        P[:, idx[to_lab[others[0]]]] = (Rr * np.cos(Pp)).ravel()
        P[:, idx[to_lab[others[1]]]] = (Rr * np.sin(Pp)).ravel()
        return P, W
    raise ValueError(kind)


def size_factor(shape, lab, coils, substrate_mm=None, n=10):
    """<G>_body / G(0) for the body mounted as `lab` (lab axis -> sample axis)."""
    P, W = body_points(shape, lab, substrate_mm, n)
    g0 = gradient(np.zeros((1, 3)), coils)[0]
    return float(np.sum(W * gradient(P, coils)) / np.sum(W) / g0)


def mallinson(radius, dx, dz):
    """Four coils on the pole caps: centres (+-dx, 0, +-dz); upper and lower pairs in series opposition."""
    return [{"centre_mm": [sx * dx, 0.0, sz * dz], "radius_mm": radius, "sense": sz}
            for sx in (1, -1) for sz in (1, -1)]
