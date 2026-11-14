"""Exact potential of point currents on a thin sheet that fills a round wafer with a primary flat.

The film covers {|z| < a, Im z > -h}: a disk of radius a cut by the flat y = -h, h = sqrt(a^2 - (L/2)^2), on an
insulating substrate, so no current crosses the boundary (Neumann condition).  The Moebius map w = (z - A)/(z - B),
with A and B the two ends of the flat, sends the domain to a wedge with its apex at w = 0 (the flat becomes the
negative real axis, the arc another ray); zeta = (w e^{-i phi0})^(pi/alpha) opens the wedge of angle alpha to the
upper half-plane, where a source at zeta0 has the Neumann Green's function -(1/2 pi)(ln|zeta - zeta0| +
ln|zeta - conj zeta0|).  Point sources are conformally invariant up to constants that cancel between two needles.
"""
import cmath
import math

import numpy as np


class FlatWafer:
    def __init__(self, radius, flat_length):
        self.a = radius
        c = flat_length / 2
        self.h = math.sqrt(radius ** 2 - c ** 2)
        self.A, self.B = complex(-c, -self.h), complex(c, -self.h)
        # the arc's direction in the w plane: any arc point maps onto one ray
        w_arc = self._mobius(complex(0.0, radius))
        self.phi_arc = cmath.phase(w_arc)
        # the domain lies between the flat (arg w = pi) and the arc ray; the centre fixes which side
        w0 = self._mobius(0j)
        lo, hi = sorted((self.phi_arc % (2 * math.pi), math.pi))
        inside = lo < cmath.phase(w0) % (2 * math.pi) < hi
        if inside:
            self.start, self.alpha = lo, hi - lo
        else:
            self.start, self.alpha = hi, 2 * math.pi - (hi - lo)

    def _mobius(self, z):
        return (z - self.A) / (z - self.B)

    def zeta(self, z):
        w = self._mobius(np.asarray(z, dtype=complex))
        arg = (np.angle(w) - self.start) % (2 * np.pi)
        return np.abs(w) ** (np.pi / self.alpha) * np.exp(1j * arg * np.pi / self.alpha)

    def pair_potential(self, point, z_in, z_out):
        """phi(point) for unit current entering at z_in and leaving at z_out, unit sheet resistance."""
        p, s_in, s_out = self.zeta(point), self.zeta(z_in), self.zeta(z_out)
        term = lambda s: np.log(np.abs(p - s)) + np.log(np.abs(p - np.conj(s)))
        return (term(s_out) - term(s_in)) / (2 * np.pi)
