"""Hull forms of the synthetic ships and their exact hydrostatics.

Frame: x forward from the aft perpendicular, y to starboard, z up from the baseline (underside of the
keel). A waterplane is z = a + bx*x + cy*y in this frame; the true vertical is along (-bx, -cy, 1).
"""
import numpy as np


class Hull:
    def __init__(self, p, nx=1600, ny=900):
        self.p = p
        self.L, self.B, self.D = p["lbp"], p["beam"], p["depth"]
        self.R = p["bilge_radius"]
        self.x0 = p["stern_tip"]            # aftmost point of the hull (deck level), x < 0
        self.x1 = p["stem_tip"]             # foremost point (deck level), x > L
        self.xs = np.linspace(self.x0, self.x1, nx + 1)
        self.xs = 0.5 * (self.xs[1:] + self.xs[:-1])
        self.dx = (self.x1 - self.x0) / nx
        ys = np.linspace(-self.B / 2, self.B / 2, ny + 1)
        self.ys = 0.5 * (ys[1:] + ys[:-1])
        self.dy = self.B / ny
        X, Y = np.meshgrid(self.xs, self.ys, indexing="ij")
        self.X, self.Y = X, Y
        self.zb = self._bottom(X, np.abs(Y))          # lowest hull point of each column (inf outside)
        self.inside = np.isfinite(self.zb)
        self.open = []                                 # spaces open to the sea
        for s in p.get("open_spaces", []):
            self.open.append(s)
        self.mask = self.inside.copy()
        for s in self.open:
            if s["kind"] == "moonpool":
                m = (np.abs(X - s["x"]) <= s["length"] / 2) & (np.abs(Y - s["y"]) <= s["breadth"] / 2)
                self.mask &= ~m
        self.zb = np.where(self.mask, self.zb, self.D)

    # ---- form -------------------------------------------------------------------------------
    def mid_half_breadth(self, z):
        B2, R = self.B / 2, self.R
        zc = np.clip(z, 0.0, R)
        return np.where(z >= R, B2, B2 - R + np.sqrt(np.maximum(R * R - (R - zc) ** 2, 0.0)))

    def ends(self, z):
        """aft and forward ends of the waterline at height z"""
        p, D = self.p, self.D
        t = np.clip(z / D, 0.0, 1.0)
        fwd = p["stem_foot"] + (self.x1 - p["stem_foot"]) * t ** p["stem_shape"]
        aft = p["stern_foot"] + (self.x0 - p["stern_foot"]) * np.minimum(1.0, z / p["stern_rise"]) ** p["stern_shape"]
        return aft, fwd

    def half_breadth(self, x, z):
        p, D = self.p, self.D
        x, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(z, float))
        bm = self.mid_half_breadth(z)
        aft, fwd = self.ends(z)
        t = np.clip(z / D, 0.0, 1.0)
        xa, xf = p["body_aft"], p["body_fwd"]
        uf = np.clip((x - xf) / np.maximum(fwd - xf, 1e-9), 0.0, None)
        ua = np.clip((xa - x) / np.maximum(xa - aft, 1e-9), 0.0, None)
        pf = p["entrance_exp"][0] + (p["entrance_exp"][1] - p["entrance_exp"][0]) * t
        pa = p["run_exp"][0] + (p["run_exp"][1] - p["run_exp"][0]) * t
        f = np.where(uf >= 1.0, 0.0, 1.0 - np.minimum(uf, 1.0) ** pf)
        a = np.where(ua >= 1.0, 0.0, 1.0 - np.minimum(ua, 1.0) ** pa)
        hb = bm * np.where(x > xf, f, 1.0) * np.where(x < xa, a, 1.0)
        return np.where((x > fwd) | (x < aft), 0.0, hb)

    def _bottom(self, X, Yabs):
        """lowest z at which the hull reaches |y| (half-breadth grows with z), by bisection"""
        lo = np.zeros_like(X)
        hi = np.full_like(X, self.D)
        top = self.half_breadth(X, hi)
        out = top < Yabs
        for _ in range(48):
            mid = 0.5 * (lo + hi)
            ok = self.half_breadth(X, mid) >= Yabs
            hi = np.where(ok, mid, hi)
            lo = np.where(ok, lo, mid)
        zb = np.where(self.half_breadth(X, np.zeros_like(X)) >= Yabs, 0.0, hi)
        return np.where(out, np.inf, zb)

    # ---- hydrostatics -----------------------------------------------------------------------
    def immersed(self, a, bx, cy):
        """volume (m3) and centroid of the hull below the waterplane, open spaces excluded"""
        zw = a + bx * self.X + cy * self.Y
        top = np.minimum(zw, self.D)
        h = np.where(self.mask, np.clip(top - self.zb, 0.0, None), 0.0)
        zlo, zhi = self.zb, self.zb + h
        hz = h * 0.5 * (zlo + zhi)
        for s in self.open:
            if s["kind"] == "tunnel":            # transverse cylinder through the hull
                d = np.sqrt(np.maximum(s["radius"] ** 2 - (self.X - s["x"]) ** 2, 0.0))
                g0 = np.maximum(s["z"] - d, zlo)
                g1 = np.minimum(s["z"] + d, zhi)
                g = np.where(self.mask & (d > 0), np.clip(g1 - g0, 0.0, None), 0.0)
                h = h - g
                hz = hz - g * 0.5 * (g0 + g1)
        dA = self.dx * self.dy
        v = h.sum() * dA
        if v <= 0:
            return 0.0, np.zeros(3)
        c = np.array([(h * self.X).sum(), (h * self.Y).sum(), hz.sum()]) * dA / v
        return v, c

    def offsets(self, stations, waterlines):
        return self.half_breadth(np.asarray(stations)[:, None], np.asarray(waterlines)[None, :])


class Box:
    """a rectangular tank; liquid surface parallel to the waterplane"""

    def __init__(self, t, n=160):
        self.t = t
        xs = np.linspace(t["x_aft"], t["x_fwd"], n + 1)
        ys = np.linspace(t["y_port"], t["y_stbd"], n + 1)
        self.X, self.Y = np.meshgrid(0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1]), indexing="ij")
        self.dA = (t["x_fwd"] - t["x_aft"]) * (t["y_stbd"] - t["y_port"]) / n / n
        self.H = t["z_top"] - t["z_bottom"]

    def fill(self, a, bx, cy):
        h = np.clip(a + bx * self.X + cy * self.Y - self.t["z_bottom"], 0.0, self.H)
        return h

    def volume(self, a, bx, cy):
        return self.fill(a, bx, cy).sum() * self.dA

    def level_for(self, vol, bx, cy):
        lo = self.t["z_bottom"] - abs(bx) * 300 - abs(cy) * 60
        hi = self.t["z_top"] + abs(bx) * 300 + abs(cy) * 60
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if self.volume(mid, bx, cy) < vol:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def liquid_group(boxes, vol, bx, cy):
    """liquid of volume vol shared by tanks connected below its surface: one common free surface"""
    lo = min(b.t["z_bottom"] for b in boxes) - 30
    hi = max(b.t["z_top"] for b in boxes) + 30
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if sum(b.volume(mid, bx, cy) for b in boxes) < vol:
            lo = mid
        else:
            hi = mid
    a = 0.5 * (lo + hi)
    m = np.zeros(3)
    v = 0.0
    for b in boxes:
        h = b.fill(a, bx, cy)
        v += h.sum() * b.dA
        m += np.array([(h * b.X).sum(), (h * b.Y).sum(), (h * (b.t["z_bottom"] + h / 2)).sum()]) * b.dA
    return a, v, m / v
