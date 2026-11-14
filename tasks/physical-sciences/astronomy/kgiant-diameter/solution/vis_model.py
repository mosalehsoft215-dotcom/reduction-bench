"""What the interferometer measures on a resolved star with a faint companion far from it (generator and reference).

Sky: uniform disc of diameter theta (mas); point companion offset by dra (east) and ddec (north), in mas, with
flux ratio f (companion / giant, the same at all wavelengths). Baseline (u, v) in metres, wavenumber s = 1 / lambda.

Instrument, as described in observing_notes.txt:
  * each telescope feeds a single-mode fibre. The mode is a Gaussian matched to the unobstructed aperture for the
    highest coupling (mode radius w with (R / w)^2 = q, exp(q) = 1 + 2 q, coupling 81.45 %). A star at angle a from
    the axis of the fibre enters with amplitude A(a) / A(0),
        A(a) = integral over the aperture of exp(-r^2 / w^2) J0(2 pi a r s) r dr,
    and with power eta = [A(a) / A(0)]^2 (coupling());
  * no atmospheric dispersion corrector: the image of a star at wavenumber s is raised toward the zenith by the
    refraction R(s, z) of the air at the site. The tip-tilt loop and the fibre injection hold the K-band (2.2 um)
    image of the observed star on the axis of every fibre. At 1 m the K-band image of the giant and its companion is
    one blend, and the centroiding tip-tilt sensor holds its centroid, which lies ck (dra, ddec) from the giant with
    ck = fK / (1 + fK) (fK the K-band flux ratio). So at wavenumber s the giant sits -ck (dra, ddec) + d (sin q, cos q)
    from the axis, d = R(s, z) - R(2.2 um, z) toward the zenith (position angle on the sky = parallactic angle q),
    and the companion at (1 - ck) (dra, ddec) + d (sin q, cos q) (Atmosphere). The fringe phase is not changed: a
    plane-stratified atmosphere keeps the component of the wave vector along the ground, so (u, v) and the
    companion's offset in the fringe term are those of the sky above the air;
  * coherent flux of a baseline in one frame, in units of the flux of the giant on the axis of a fibre:
        C = eta_g 2 J1(x) / x + f eta_c exp(-2 pi i s (u dra + v ddec)),   x = pi theta B s,   B = hypot(u, v);
    flux of a telescope: eta_g + f eta_c. The same for all six telescopes (same refraction and same guiding);
  * a spectral channel is a top-hat in wavenumber: the detector adds C and the telescope flux over the channel;
  * an observation lasts int_time seconds while the Earth turns: the reduction averages |C|^2 and the triple
    products C12 C23 C31 over the frames. Squared visibility = mean |C|^2 / (mean flux of a telescope)^2,
    closure phase = argument of the mean triple product.
The file gives (u, v) at the middle of an observation; a row also needs how they move (track()) and, for the
refraction, the hour angle and declination of date of the star at that moment (keys "ha", "dec", radians).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0, j1, jv

MAS = np.pi / 180 / 3600 / 1000          # radians per mas
SIDEREAL = 2 * np.pi / 86164.0905        # rotation of the Earth, rad / s
APERTURE = 1.0                           # telescope diameter, m
MODE_Q = 1.2564312086261697              # (aperture radius / mode radius)^2 at the highest coupling
GUIDE_UM = 2.2                           # wavelength of the image held on the fibre axis (K band)


def nodes(n):
    """Gauss-Legendre nodes and weights on [-1/2, 1/2]; one node at 0 for n = 1"""
    if n == 1:
        return np.zeros(1), np.ones(1)
    x, w = leggauss(n)
    return x / 2, w / 2


def disc(theta, b, s):
    x = np.pi * theta * MAS * b * s
    return np.where(x > 1e-9, 2 * j1(x) / np.where(x > 1e-9, x, 1.0), 1.0)


def ld_disc(theta, b, s, u_ld):
    """visibility of a linearly limb-darkened disc I(mu) = 1 - u (1 - mu), limb-darkened diameter theta (mas)"""
    x = np.maximum(np.pi * theta * MAS * b * s, 1e-9)
    v = ((1 - u_ld) * j1(x) / x + u_ld * np.sqrt(np.pi / 2) * jv(1.5, x) / x ** 1.5) / ((1 - u_ld) / 2 + u_ld / 3)
    return np.where(np.pi * theta * MAS * b * s > 1e-9, v, 1.0)


def precession_matrix(mjd):
    """rotation from the mean equator and equinox of J2000 to those of the date (IAU 1976 angles)"""
    t = (mjd - 51544.5) / 36525.0
    a = np.radians(np.array([2306.2181 * t + 0.30188 * t ** 2 + 0.017998 * t ** 3,
                             2306.2181 * t + 1.09468 * t ** 2 + 0.018203 * t ** 3,
                             2004.3109 * t - 0.42665 * t ** 2 - 0.041833 * t ** 3]) / 3600.0)
    zeta, z, th = a
    rz = lambda g: np.array([[np.cos(g), -np.sin(g), 0], [np.sin(g), np.cos(g), 0], [0, 0, 1]])
    ry = lambda g: np.array([[np.cos(g), 0, -np.sin(g)], [0, 1, 0], [np.sin(g), 0, np.cos(g)]])
    return rz(z) @ ry(th) @ rz(zeta)


def _unit(ra, dec):
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def radec_of_date(ra2000, dec2000, mjd):
    """right ascension and declination (radians) on the mean equator and equinox of the date"""
    q = precession_matrix(mjd) @ _unit(ra2000, dec2000)
    return np.arctan2(q[1], q[0]) % (2 * np.pi), np.arcsin(q[2])


def _offset(p0, p1):
    """east and north offset (radians, tangent plane) of direction p1 from direction p0"""
    ra, dec = np.arctan2(p0[1], p0[0]), np.arcsin(p0[2])
    e = np.array([-np.sin(ra), np.cos(ra), 0.0])
    n = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])
    return p1 @ e / (p1 @ p0), p1 @ n / (p1 @ p0)


def offsets_to_date(east, north, ra2000, dec2000, mjd):
    """companion offset (mas) given on the J2000 frame -> the same offset on the frame of the date"""
    p0 = _unit(ra2000, dec2000)
    e = np.array([-np.sin(ra2000), np.cos(ra2000), 0.0])
    n = np.array([-np.sin(dec2000) * np.cos(ra2000), -np.sin(dec2000) * np.sin(ra2000), np.cos(dec2000)])
    p1 = p0 + east * MAS * e + north * MAS * n
    P = precession_matrix(mjd)
    de, dn = _offset(P @ p0, P @ (p1 / np.linalg.norm(p1)))
    return de / MAS, dn / MAS


def offsets_to_j2000(east, north, ra2000, dec2000, mjd):
    """companion offset (mas) measured on the frame of the date -> the J2000 frame"""
    P = precession_matrix(mjd)
    q0 = P @ _unit(ra2000, dec2000)
    ra, dec = np.arctan2(q0[1], q0[0]), np.arcsin(q0[2])
    e = np.array([-np.sin(ra), np.cos(ra), 0.0])
    n = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])
    q1 = q0 + east * MAS * e + north * MAS * n
    de, dn = _offset(P.T @ q0, P.T @ (q1 / np.linalg.norm(q1)))
    return de / MAS, dn / MAS


def gmst(mjd):
    """Greenwich mean sidereal time (radians) at UTC mjd (UT1 = UTC; IAU 1982 expression)"""
    d = np.asarray(mjd, float) - 51544.5
    t = d / 36525.0
    deg = 280.46061837 + 360.98564736629 * d + 0.000387933 * t ** 2 - t ** 3 / 38710000.0
    return np.radians(deg % 360.0)


def coupling(rho, s, n=48):
    """light of a star rho mas off the axis that enters the fibre, relative to a star on the axis"""
    x, w = leggauss(n)
    x, w = (x + 1) / 2, w / 2
    g = w * x * np.exp(-MODE_Q * x ** 2)
    z = np.pi * APERTURE * np.asarray(rho, float) * MAS * np.asarray(s, float)
    return ((g * j0(z[..., None] * x)).sum(-1) / g.sum()) ** 2


# ---- refraction ------------------------------------------------------------------------------------------------
def saturation_pressure(t_c):
    """saturation vapour pressure over water, Pa (IAPWS form used by Ciddor 1996 and NIST)"""
    T = t_c + 273.15
    return np.exp(1.2378847e-5 * T ** 2 - 1.9121316e-2 * T + 33.93711047 - 6.3431645e3 / T)


def ciddor(lam_um, t_c, p_pa, rh, xc=450.0):
    """refractivity n - 1 of moist air (Ciddor 1996, Appl. Opt. 35, 1566) at vacuum wavelength lam_um"""
    s2 = 1.0 / np.asarray(lam_um, float) ** 2
    T = t_c + 273.15
    r_as = 1e-8 * (5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2))
    r_vs = 1.022e-8 * (295.235 + 2.6422 * s2 - 0.032380 * s2 ** 2 + 0.004028 * s2 ** 3)
    ma = 0.0289635 + 1.2011e-8 * (xc - 400.0)
    r_axs = r_as * (1 + 5.34e-7 * (xc - 450.0))
    R = 8.314472
    rho_axs = 101325.0 * ma / (0.9995922115 * R * 288.15)
    rho_vs = 0.00985938
    f = 1.00062 + 3.14e-8 * p_pa + 5.6e-7 * t_c ** 2
    xv = rh * f * saturation_pressure(t_c) / p_pa
    zm = (1 - p_pa / T * (1.58123e-6 - 2.9331e-8 * t_c + 1.1043e-10 * t_c ** 2 + (5.707e-6 - 2.051e-8 * t_c) * xv
                          + (1.9898e-4 - 2.376e-6 * t_c) * xv ** 2) + (p_pa / T) ** 2 * (1.83e-11 - 0.765e-8 * xv ** 2))
    rho_a = p_pa * ma * (1 - xv) / (zm * R * T)
    rho_v = xv * p_pa * 0.018015 / (zm * R * T)
    return rho_a / rho_axs * r_axs + rho_v / rho_vs * r_vs


def zenith_parallactic(ha, dec, lat):
    """zenith distance and parallactic angle (position angle of the direction to the zenith), radians"""
    cz = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(ha)
    q = np.arctan2(np.sin(ha), np.tan(lat) * np.cos(dec) - np.sin(dec) * np.cos(ha))
    return np.arccos(np.clip(cz, -1, 1)), q


class Atmosphere:
    """refraction at the site (Stone 1996, PASP 108, 1051, with true zenith distance) and the displacement of the
    image at wavenumber s from the image at the guiding wavelength, toward the zenith"""

    def __init__(self, lat, height_m, t_c, p_hpa, rh, guide_um=GUIDE_UM, index=ciddor, sign=1.0, mid=False, giant=True,
                 ck=0.0):
        """index(lambda_um, t_c, p_pa, rh) -> n - 1. Options for alternative analyses: sign -1 moves the images away
        from the zenith, mid evaluates the refraction at the middle of each observation only, giant False leaves the
        giant on the axis at every wavenumber (only the companion moves); ck is the offset of the guided K-band
        centroid from the giant in units of the companion's offset."""
        self.lat, self.h, self.t, self.p, self.rh = lat, height_m, t_c, 100.0 * p_hpa, rh
        self.guide, self.index, self.sign, self.mid, self.giant, self.ck = guide_um, index, sign, mid, giant, ck
        self.kappa = 1 + 0.005302 * np.sin(lat) ** 2 - 0.00000583 * np.sin(2 * lat) ** 2 - 0.000000315 * height_m
        self.beta = 0.001254 * (t_c + 273.15) / 273.15

    def refraction(self, s, z):
        g = self.index(1e6 / np.asarray(s, float), self.t, self.p, self.rh)
        tz = np.tan(z)
        return self.kappa * g * ((1 - self.beta) * tz - (self.beta - g / 2) * tz ** 3)

    def offset(self, ha, dec, s):
        """east, north displacement (mas) of the image at wavenumber s; ha and dec broadcast against s"""
        z, q = zenith_parallactic(ha, dec, self.lat)
        d = self.sign * (self.refraction(s, z) - self.refraction(1e6 / self.guide, z)) / MAS
        return d * np.sin(q), d * np.cos(q)


def track(tr, tau):
    """(u, v) of each row tau seconds from the middle of its observation: arrays (rows, len(tau)).
    u = a sin(Wt) + b cos(Wt) and v = vc + d sin(Wt) + e cos(Wt); up and vp are du/dt / W and dv/dt / W."""
    c, s = np.cos(SIDEREAL * tau), np.sin(SIDEREAL * tau)
    u = tr["u"][:, None] * c + tr["up"][:, None] * s
    v = tr["vc"][:, None] + (tr["v"] - tr["vc"])[:, None] * c + tr["vp"][:, None] * s
    return u, v


def couplings(tr, s, tau, dra, ddec, lobe=coupling, atm=None):
    """coupled power of the giant and of the companion, arrays (rows, instants, wavenumbers). Without atm the giant
    is on the axis at every wavenumber. Rows of the same observation and channel are computed once."""
    rows = s.shape[0]
    if atm is None:
        g = lobe(np.hypot(dra, ddec), s)
        return np.ones((rows, 1, s.shape[1])), g[:, None, :]
    key = np.c_[tr["ha"], tr["dec"], s[:, 0]]
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    tau = np.asarray(tau, float)
    tau = tau[first] if tau.ndim == 2 else tau[None, :]
    ha = tr["ha"][first][:, None, None] + (0.0 if atm.mid else SIDEREAL * tau[:, :, None])
    su = s[first][:, None, :]
    de, dn = atm.offset(ha, tr["dec"][first][:, None, None], su)
    de, dn = np.broadcast_to(de, np.broadcast_shapes(de.shape, su.shape)), np.broadcast_to(dn, np.broadcast_shapes(dn.shape, su.shape))
    sb = np.broadcast_to(su, de.shape)
    ck = atm.ck
    eg = lobe(np.hypot(de - ck * dra, dn - ck * ddec), sb) if atm.giant else np.ones(de.shape)
    ec = lobe(np.hypot((1 - ck) * dra + de, (1 - ck) * ddec + dn), sb)
    return eg[inv], ec[inv]


def coherent(tr, s_lo, s_hi, int_time, p, kt=16, ks=12, lobe=coupling, star_vis=None, atm=None):
    """channel sums of the coherent flux at kt instants of each observation (rows, kt), channel sums of the flux of
    a telescope (rows, kt), weights of the instants. kt = 1: baseline at rest; ks = 1: one wavenumber per channel."""
    theta, dra, ddec, f = p
    tx, tw = nodes(kt)
    sx, sw = nodes(ks)
    u, v = track(tr, tx * int_time)
    s = ((s_lo + s_hi) / 2)[:, None] + sx * (s_hi - s_lo)[:, None]
    eg, ec = couplings(tr, s, tx * int_time, dra, ddec, lobe, atm)
    phase = 2 * np.pi * (u * dra + v * ddec) * MAS
    sv = star_vis or (lambda b, s_: disc(theta, b, s_))
    star = (sw * eg * sv(np.hypot(u, v)[:, :, None], s[:, None, :])).sum(-1)
    comp = (sw * f * ec * np.exp(-1j * phase[:, :, None] * s[:, None, :])).sum(-1)
    flux = np.broadcast_to((sw * (eg + f * ec)).sum(-1), star.shape)
    return star + comp, flux, tw


def vis2(d, p, **kw):
    c, flux, tw = coherent(d, d["s_lo"], d["s_hi"], d["T"], p, **kw)
    return (tw * np.abs(c) ** 2).sum(-1) / (tw * flux).sum(-1) ** 2


KEYS = ("u", "v", "up", "vp", "vc")


def closure(d, p, **kw):
    """closure phase in degrees of the triangle with baselines 1, 2 and 3 = -(1 + 2)"""
    one, two = {k: d[k + "1"] for k in KEYS}, {k: d[k + "2"] for k in KEYS}
    three = {k: -one[k] - two[k] for k in KEYS}
    for leg in (one, two, three):
        leg.update({k: d[k] for k in ("ha", "dec") if k in d})
    a, _, tw = coherent(one, d["s_lo"], d["s_hi"], d["T"], p, **kw)
    b = coherent(two, d["s_lo"], d["s_hi"], d["T"], p, **kw)[0]
    c = coherent(three, d["s_lo"], d["s_hi"], d["T"], p, **kw)[0]
    return np.degrees(np.angle((tw * a * b * c).sum(-1)))
