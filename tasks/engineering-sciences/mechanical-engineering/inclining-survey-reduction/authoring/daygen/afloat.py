"""Floating equilibrium of a ship with its weights, tank liquids and hanging loads."""
import numpy as np
from scipy.optimize import fsolve

from hull import Box, liquid_group


class Ship:
    def __init__(self, hull, tanks, appendages):
        self.hull = hull
        self.boxes = {t["name"]: Box(t) for t in tanks}
        self.tanks = {t["name"]: t for t in tanks}
        self.app = appendages               # [(volume m3, (x, y, z))], always immersed

    def groups(self, vols):
        """tanks joined through an open cross-connection share one free surface"""
        seen, out = set(), []
        for name in vols:
            if name in seen:
                continue
            grp = [name] + [m for m in self.tanks[name].get("joined", [])]
            seen.update(grp)
            out.append((grp, sum(vols[g] for g in grp)))
        return out

    def weight_and_centre(self, plane, items, vols, hanging):
        a, bx, cy = plane
        n = np.array([-bx, -cy, 1.0]) / np.sqrt(1 + bx * bx + cy * cy)
        m = sum(i[0] for i in items)
        mom = sum(np.array(i[1:4]) * i[0] for i in items) if items else np.zeros(3)
        for grp, vol in self.groups(vols):
            rho = self.tanks[grp[0]]["density"]
            _, v, c = liquid_group([self.boxes[g] for g in grp], vol, bx, cy)
            m += rho * v
            mom = mom + rho * v * c
        for w, s, length in hanging:      # hangs along the true vertical below its point of suspension
            p = np.array(s) - length * n
            m += w
            mom = mom + w * p
        return m, mom / m

    def buoyancy(self, plane, rho):
        v, c = self.hull.immersed(*plane)
        mom = v * c
        for va, ca in self.app:
            v += va
            mom = mom + va * np.array(ca)
        return rho * v, mom / v

    def residual(self, plane, rho, items, vols, hanging, m_heel=0.0):
        mw, g = self.weight_and_centre(plane, items, vols, hanging)
        mb, b = self.buoyancy(plane, rho)
        _, bx, cy = plane
        # an outside heeling moment (wind on the hull, t.m, to starboard) acts like G moved by m_heel / mass
        return np.array([(mb - mw) / mw * 100.0,
                         (g[0] + bx * g[2]) - (b[0] + bx * b[2]),
                         (g[1] + cy * g[2]) + m_heel / mw - (b[1] + cy * b[2])])

    def float_(self, rho, items, vols, hanging=(), guess=None, m_heel=0.0):
        if guess is None:
            guess = np.array([self.hull.D * 0.4, 0.0, 0.0])
        sol, info, ier, msg = fsolve(lambda q: self.residual(q, rho, items, vols, hanging, m_heel), guess,
                                     full_output=True, xtol=1e-12, epsfcn=1e-10)
        r = self.residual(sol, rho, items, vols, hanging, m_heel)
        if np.max(np.abs(r)) > 1e-6:
            raise RuntimeError(f"no equilibrium: {r} {msg}")
        return sol
