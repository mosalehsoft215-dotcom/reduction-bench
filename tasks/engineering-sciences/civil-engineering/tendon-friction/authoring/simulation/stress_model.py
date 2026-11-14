"""Tendon forces and strand extensions for stressing campaigns on continuous post-tensioned box girders.

Every duct is built as a curve in space from the drawings:
  - the girder axis in plan is a chain of straights and circular arcs (left or right) by station, the station being
    the horizontal distance along the axis from anchorage A;
  - the soffit at the axis follows the vertical alignment: level, or a crest or sag curve (a parabola in station of
    the given radius between two stations, zero grade before it, constant grade after it);
  - a web duct lies in the mid-plane of its web; the webs lean outward by web_tan horizontal per unit height and their
    mid-planes meet the soffit at web_y0 on either side of the axis, so a duct at height h above the soffit sits
    web_y0 + h web_tan from the axis, on the side of its web; a bottom-slab duct keeps its lateral offset;
  - heights above the soffit are vertical; between two control points the height is one parabola in station.
The strands do not run along the duct centre: under tension they take the shortest path the duct leaves them, a
taut line inside the tube of radius e (the eccentricity of the strand bundle in its duct) around the duct centre,
centred at both anchor heads. Where the duct bends they bear on the inside of the bend; where the bend reverses, or
a straight runs into a bend, they cross the duct on a straight chord and turn less than the duct does.
The friction angle theta(s) is the turn of the tangent of that strand path from anchorage A, accumulated over its arc
length s, plus the local kink of a duct damaged during concreting (an angle at one station that no drawing shows).
Friction: P(s) = P0 exp(-(mu theta + K s)). Stressing from one end: the strand moves where the jacking curve exceeds
the force already in it. Lock-off: the wedges draw in; near the anchor the force then follows the mirrored friction
curve. The ram travel also contains the elongation, without friction, of the strand inside the jack between the
anchor head and the gripping wedges.
"""
import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import splu

DS = 0.02                       # station step (m)
STRAIGHT = 6.0                  # straight length behind each anchor (m)
PATH_EVERY = 5                  # the strand path is solved on every 5th station (0.1 m)


def slope_knots(spans, slopes):
    """knots (station, slope) of the profile of a draped tendon over any number of spans: straight for STRAIGHT
    metres behind each anchor at slope sa, slope si at the inflection points of the end spans, sm at the inflection
    points of the inner spans, zero slope at the low points and over the piers"""
    sa, si, sm = slopes
    ends = np.concatenate([[0.0], np.cumsum(spans)])
    xs, ss = [0.0, STRAIGHT, 0.4 * spans[0], 0.9 * spans[0], ends[1]], [-sa, -sa, 0.0, si, 0.0]
    for i in range(1, len(spans) - 1):
        a, length = ends[i], spans[i]
        xs += [a + 0.1 * length, a + 0.5 * length, a + 0.9 * length, ends[i + 1]]
        ss += [-sm, 0.0, sm, 0.0]
    a, length = ends[-2], spans[-1]
    xs += [a + 0.1 * length, ends[-1] - 0.4 * length, ends[-1] - STRAIGHT, ends[-1]]
    ss += [-si, 0.0, sa, sa]
    return np.array(xs), np.array(ss)


def knot_names(n_spans):
    names = ["anchor A", "end of straight", "low point", "inflection", "pier 1"]
    for i in range(1, n_spans - 1):
        names += ["inflection", "low point", "inflection", f"pier {i + 1}"]
    return names + ["inflection", "low point", "end of straight", "anchor B"]


def cumtrapz(f, x):
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


def plan_curvature(plan, st):
    """signed plan curvature on the station grid from a chain of segments (length_m, radius_m or None, +1 left/-1 right)"""
    k, start = np.zeros_like(st), 0.0
    for length, radius, turn in plan:
        if radius:
            k[(st >= start) & (st < start + length)] = turn / radius
        start += length
    return k


def soffit_level(vertical, st):
    """level of the soffit at the axis (m) along the station grid"""
    if not vertical:
        return np.zeros_like(st)
    a, b, r = vertical["from"], vertical["to"], vertical["radius"]
    sign = -1.0 if vertical["kind"] == "crest" else 1.0
    u = np.clip(st, a, b) - a
    return sign * (u ** 2 / (2.0 * r) + np.where(st > b, (st - b) * (b - a) / r, 0.0))


def heights(spans, slopes, h0):
    st = np.arange(0.0, float(sum(spans)) + DS / 2, DS)
    slope = np.interp(st, *slope_knots(spans, slopes))
    return st, h0 + cumtrapz(slope, st)


def curve(br, td):
    """station grid, 3-D points (n, 3) of the duct centre"""
    st, h = heights(br["spans"], td["slopes"], td["h0"])
    psi = cumtrapz(plan_curvature(br["plan"], st), st)
    px, py = cumtrapz(np.cos(psi), st), cumtrapz(np.sin(psi), st)
    if td["family"].startswith("web"):
        y = td["side"] * (br["section"]["web_y0"] + h * br["section"]["web_tan"])
    else:
        y = np.full_like(st, td["y"])
    z = soffit_level(br["vertical"], st) + h
    return st, np.stack([px - y * np.sin(psi), py + y * np.cos(psi), z], 1)


def strand_path(c, ecc, rho=0.05, tol=1e-10, iters=20000):
    """shortest path through the tube of radius ecc around the polyline c (n, 3), ends at the centre: each node moves
    in the plane normal to the duct; the sum of squared chord lengths is minimised by ADMM with the disc constraint"""
    n = len(c)
    t = np.gradient(c, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    n1 = np.cross(np.array([0.0, 0.0, 1.0]), t)
    n1 /= np.linalg.norm(n1, axis=1)[:, None]
    n2 = np.cross(t, n1)
    m = n - 2
    i = np.arange(n - 1)
    rows, cols, vals = [], [], []
    for k in range(3):
        r = 3 * i + k
        nxt = (i + 1 >= 1) & (i + 1 <= n - 2)
        prv = (i >= 1) & (i <= n - 2)
        for a, use, j, sign in ((n1, nxt, i, 1.0), (n2, nxt, i, 1.0), (n1, prv, i - 1, -1.0), (n2, prv, i - 1, -1.0)):
            node = j + 1 if sign > 0 else j + 1
            col = 2 * j + (0 if a is n1 else 1)
            rows.append(r[use]); cols.append(col[use]); vals.append(sign * a[node[use], k])
    b = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * (n - 1), 2 * m))
    g = b.T @ np.diff(c, axis=0).ravel()
    lu = splu(((b.T @ b) + rho * sparse.identity(2 * m)).tocsc())
    z, y = np.zeros(2 * m), np.zeros(2 * m)
    for it in range(iters):
        u = lu.solve(-g + rho * (z - y))
        w = (u + y).reshape(-1, 2)
        r = np.linalg.norm(w, axis=1)
        zn = (w * np.minimum(1.0, ecc / np.maximum(r, 1e-15))[:, None]).ravel()
        y += u - zn
        done = it > 50 and np.abs(zn - z).max() < tol and np.abs(u - zn).max() < tol
        z = zn
        if done:
            break
    off = np.zeros((n, 2))
    off[1:-1] = z.reshape(-1, 2)
    return c + off[:, :1] * n1 + off[:, 1:] * n2


def turned(pts):
    """arc length and turned angle along a polyline (the turn between two chords is placed at their common node)"""
    d = np.diff(pts, axis=0)
    seg = np.linalg.norm(d, axis=1)
    t = d / seg[:, None]
    turn = np.arctan2(np.linalg.norm(np.cross(t[1:], t[:-1]), axis=1), np.sum(t[1:] * t[:-1], 1))
    return np.concatenate([[0.0], np.cumsum(seg)]), np.concatenate([[0.0, 0.0], np.cumsum(turn)])


def geometry(br, td, kink=None):
    """station grid, arc length s from A and turned angle theta from A along the strands: on the duct centre where
    the bundle has no eccentricity, otherwise along the strand path solved on every PATH_EVERY-th station and
    carried to the station grid by linear interpolation"""
    st, pts = curve(br, td)
    ecc = br.get("ducts", {}).get("ecc", 0.0)
    if ecc > 0.0:
        sub = np.arange(0, len(st), PATH_EVERY)
        if sub[-1] != len(st) - 1:
            sub = np.append(sub, len(st) - 1)
        s_sub, th_sub = turned(strand_path(pts[sub], ecc))
        s, theta = np.interp(st, st[sub], s_sub), np.interp(st, st[sub], th_sub)
    else:
        s, theta = turned(pts)
    if kink:                                  # (station, angle): a local deviation of a damaged duct
        theta = theta + kink[1] * (st >= kink[0])
    return st, s, theta


def jack_from(force, s, theta, mu, k, end):
    if end == "A":
        return force * np.exp(-(mu * theta + k * s))
    return force * np.exp(-(mu * (theta[-1] - theta) + k * (s[-1] - s)))


def integrate(f, s):
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(s)))


def stress(existing, force, s, theta, mu, k, end):
    """force profile after jacking at `end` to `force`, and the strand drawn out of the duct (N m, divide by EA)"""
    new = np.maximum(existing, jack_from(force, s, theta, mu, k, end))
    return new, integrate(new - existing, s)


def lock_off(force_profile, draw_in, ea, s, end):
    """force profile after the wedges draw in by `draw_in` (m) at `end`: near the anchor the force follows the
    mirrored friction curve p_i^2 / p up to the node i where the strand pulled back equals the draw-in"""
    p = force_profile if end == "A" else force_profile[::-1]
    u = s if end == "A" else (s[-1] - s)[::-1]

    slack = cumtrapz(p, u) - p ** 2 * cumtrapz(1.0 / p, u)       # strand pulled back if the curve turned at i
    hit = np.nonzero(slack[1:] >= draw_in * ea)[0]
    i = hit[0] + 1 if hit.size else p.size - 1
    out = p.copy()
    out[: i + 1] = p[i] ** 2 / p[: i + 1]
    return out if end == "A" else out[::-1]
