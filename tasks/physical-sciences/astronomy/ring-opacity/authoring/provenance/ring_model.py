"""Forward model of the ring occultation, shared by the generator and the reference.

Sky plane: xi (east) and eta (north), km at the distance of the body, measured from the body's
ephemeris position. A station sees the star move on a straight line, (xi0 + vxi t, eta0 + veta t).
Ring: its mid-line is a circle of radius a about a centre (xic, etac) in a plane seen at opening angle B;
the apparent major axis has position angle P (north through east). Ring-plane polar coordinates of a sky
point: x along the major axis, y along the minor axis, rho = hypot(x, y / sin B), theta = atan2(y / sin B, x).
Radial width W(theta) = w0 [1 - q cos(theta - theta_p)], edges at a -+ W / 2.

Direct beam (star_flux). Where the line of sight crosses the ring a fraction T(theta) = exp(-tau_los(theta))
of the starlight is left in the direct beam. The band is gray with straight edges; d = (rho - r_edge) sin B /
sqrt(sin^2 B cos^2 theta + sin^2 theta) is the distance of the star from an edge at right angles to it in
the sky, F = sqrt(lambda D / 2), amplitude 1 - (1 - sqrt T) [Q(d_in / F) - Q(d_out / F)] with
Q(u) = [(1/2 + C(u)) + i (1/2 + S(u))] / (1 + i). The flux is |amplitude|^2 averaged over the passband
and over the exposure, which starts at the time stamp.

Diffracted light (halo). A particle of radius r, opaque and far larger than the wavelength, stops the light
that falls on it and diffracts as much again into the pattern of a disc of its size. Along a line of
sight that meets x = tau / sin B particle cross-sections on average, being stopped and being diffracted
are two independent Poisson processes of mean x each: exp(-2x) of the light is left in the direct beam and
exp(-2x) x^k / k! has been diffracted k times without being stopped. One diffraction turns the ray by an
angle drawn from the pattern of the disc, [2 J1(z) / z]^2 with z = 2 pi r alpha / lambda, the radius
being drawn in proportion to the cross-section the particles of that size offer; k diffractions convolve
the pattern k times with itself. At the distance D an angle alpha is a displacement alpha D, so a station
receives from every piece of the ring the light that it has turned by the right angle:
    halo(p) = sum_k  integral over the ring  exp(-2x) x^k / k!  P_k(|p - p'|)  dA(p')
with P_k the k-fold pattern per unit area of the sky plane. The Fourier transform of the single pattern
is the overlap of two discs, 2/pi [acos(v) - v sqrt(1 - v^2)] with v = s lambda D / (2 r) at spatial
frequency s; the transform of the k-fold pattern is its k-th power.
"""
import numpy as np
from scipy.special import factorial, fresnel, j0

AU_KM = 149597870.7
GL_X, GL_W = np.polynomial.legendre.leggauss(128)                # exposure average of the reference, see README
GL_X, GL_W = (GL_X + 1) / 2, GL_W / 2                      # nodes and weights on [0, 1]
ORDERS = 12                                                # diffractions kept; x^13 / 13! is below 1e-6 for x < 2
REACH_KM = 700.0                                           # the patterns are followed this far from a ring element


def q_edge(u):
    s, c = fresnel(u)
    return ((0.5 + c) + 1j * (0.5 + s)) / (1 + 1j)


def ring_plane(xi, eta, centre, sin_b, pa):
    dx, dy = xi - centre[0], eta - centre[1]
    x = dx * np.sin(pa) + dy * np.cos(pa)
    y = (-dx * np.cos(pa) + dy * np.sin(pa)) / sin_b
    return np.hypot(x, y), np.arctan2(y, x)


def width(theta, w0, q, theta_p):
    return w0 * (1 - q * np.cos(theta - theta_p))


def star_flux(t, station, ring, lam_km, lam_w, dist_km, nodes=None, rows=None):
    """direct beam of a point source (unit flux outside the ring) at the time stamps t of one station.
    station: dict(xi0, eta0, vxi, veta, exposure); ring: dict(centre, a, sin_b, pa, w0, q, theta_p, tau_los)
    where tau_los(theta) gives the line-of-sight optical depth of the direct beam at ring longitude theta.
    nodes: (positions in the exposure from 0 to 1, weights); default 128-point Gauss-Legendre.
    rows: None, or a boolean array over t: the samples that get the diffraction integral; the others get
    geometric optics (1 outside the ring, T inside). The generator never uses it."""
    gx, gw = (GL_X, GL_W) if nodes is None else nodes
    te = t[:, None] + station["exposure"] * gx[None, :]
    xi = station["xi0"] + station["vxi"] * te
    eta = station["eta0"] + station["veta"] * te
    sb = ring["sin_b"]
    rho, th = ring_plane(xi, eta, ring["centre"], sb, ring["pa"])
    w = width(th, ring["w0"], ring["q"], ring["theta_p"])
    g = sb / np.sqrt(sb ** 2 * np.cos(th) ** 2 + np.sin(th) ** 2)
    d_in, d_out = (rho - ring["a"] + w / 2) * g, (rho - ring["a"] - w / 2) * g
    tau = ring["tau_los"](th)
    if rows is None:
        rows = np.ones(t.size, bool)
    out = np.where((d_in > 0) & (d_out < 0), np.exp(-tau), 1.0) @ gw
    d_in, d_out, block = d_in[rows], d_out[rows], 1 - np.exp(-tau[rows] / 2)
    flux = np.zeros_like(d_in)
    for lam, wt in zip(lam_km, lam_w):
        f = np.sqrt(lam * dist_km / 2)
        flux += wt * np.abs(1 - block * (q_edge(d_in / f) - q_edge(d_out / f))) ** 2
    out[rows] = flux @ gw
    return out


def near_edges(t, station, ring, reach_km):
    """the samples whose exposure comes within reach_km of an edge of the ring (distance in the sky)"""
    te = t[:, None] + station["exposure"] * np.array([0.0, 0.5, 1.0])[None, :]
    sb = ring["sin_b"]
    rho, th = ring_plane(station["xi0"] + station["vxi"] * te, station["eta0"] + station["veta"] * te, ring["centre"], sb, ring["pa"])
    w = width(th, ring["w0"], ring["q"], ring["theta_p"])
    g = sb / np.sqrt(sb ** 2 * np.cos(th) ** 2 + np.sin(th) ** 2)
    d_in, d_out = (rho - ring["a"] + w / 2) * g, (rho - ring["a"] - w / 2) * g
    return (np.minimum(np.abs(d_in), np.abs(d_out)).min(1) < reach_km) | (np.ptp(np.sign(d_in), axis=1) > 0) | (np.ptp(np.sign(d_out), axis=1) > 0)


def disc_overlap(v):
    v = np.clip(v, 0.0, 1.0)
    return 2 / np.pi * (np.arccos(v) - v * np.sqrt(1 - v * v))


def patterns(r_min_km, r_max_km, slope, lam_km, lam_w, dist_km, orders=ORDERS, n_rho=1400, n_s=1600, n_r=48):
    """P_k(rho) for k = 1..orders, per km^2 of the sky plane, averaged over the passband.
    The particles have radii from r_min to r_max with dn/dr ~ r^-slope; a diffraction happens at a
    particle of radius r in proportion to r^2 dn/dr. Returns the rho grid (km) and an array (orders, n_rho)."""
    rho = REACH_KM * (np.arange(n_rho) / (n_rho - 1.0)) ** 2
    x, w = np.polynomial.legendre.leggauss(n_r)
    lr = 0.5 * (np.log(r_max_km) + np.log(r_min_km)) + 0.5 * (np.log(r_max_km) - np.log(r_min_km)) * x
    r = np.exp(lr)
    wr = w * r ** (3.0 - slope)                                 # r^2 dn/dr dr = r^(3 - slope) d ln r
    wr /= wr.sum()
    x, w = np.polynomial.legendre.leggauss(n_s)
    sig = r_max_km * (x + 1)                                    # sigma = s lambda D, from 0 to 2 r_max (km)
    ws = r_max_km * w
    m = (wr[:, None] * disc_overlap(sig[None, :] / (2 * r[:, None]))).sum(0)
    powers = np.array([m ** k for k in range(1, orders + 1)]) * (2 * np.pi * ws * sig)[None, :]
    out = np.zeros((orders, n_rho))
    for lam, wt in zip(lam_km, lam_w):
        ld = lam * dist_km
        out += wt / ld ** 2 * (powers @ j0(2 * np.pi * np.outer(sig, rho / ld)))
    return rho, out


def ring_elements(ring, n_theta=1440, n_r=6):
    """sky positions (km), sky areas (km^2) and line-of-sight particle counts x of the pieces of the ring"""
    sb, pa = ring["sin_b"], ring["pa"]
    th = (np.arange(n_theta) + 0.5) * 2 * np.pi / n_theta
    w = width(th, ring["w0"], ring["q"], ring["theta_p"])
    gx, gw = np.polynomial.legendre.leggauss(n_r)
    rr = ring["a"] + w[:, None] * gx[None, :] / 2
    area = sb * rr * w[:, None] * gw[None, :] / 2 * (2 * np.pi / n_theta)
    x, y = rr * np.cos(th)[:, None], rr * np.sin(th)[:, None] * sb
    xi = ring["centre"][0] + x * np.sin(pa) - y * np.cos(pa)
    eta = ring["centre"][1] + x * np.cos(pa) + y * np.sin(pa)
    count = np.repeat((ring["tau_los"](th) / 2)[:, None], n_r, axis=1)
    return xi.ravel(), eta.ravel(), area.ravel(), count.ravel()


def halo(xi, eta, ring, pat, n_theta=1440, n_r=6, orders=None, share=None, chunk=400):
    """light diffracted by the ring particles that reaches the sky-plane points (xi, eta), in units of the
    unocculted flux of the star. pat: the output of patterns(). orders: number of diffractions kept (default:
    all in pat). share: replaces the Poisson shares exp(-2x) x^k / k!; a function of the particle count x of
    each ring element that returns an array (orders, elements). Only calibrate.py uses the last two."""
    rho, pk = pat
    n_k = pk.shape[0] if orders is None else orders
    ex, ey, area, count = ring_elements(ring, n_theta, n_r)
    near = np.hypot(ex - np.mean(xi), ey - np.mean(eta)) < REACH_KM + np.hypot(xi - np.mean(xi), eta - np.mean(eta)).max()
    ex, ey, area, count = ex[near], ey[near], area[near], count[near]
    if share is None:
        k = np.arange(1, n_k + 1)
        src = np.exp(-2 * count)[None, :] * count[None, :] ** k[:, None] / factorial(k)[:, None] * area[None, :]
    else:
        src = share(count) * area[None, :]
        n_k = src.shape[0]
    xi, eta = np.atleast_1d(xi).astype(float), np.atleast_1d(eta).astype(float)
    out = np.zeros(xi.size)
    for i in range(0, xi.size, chunk):
        d = np.hypot(xi[i:i + chunk, None] - ex[None, :], eta[i:i + chunk, None] - ey[None, :])
        for j in range(n_k):
            out[i:i + chunk] += np.interp(d, rho, pk[j], right=0.0) @ src[j]
    return out
