"""A second reduction written apart from the reference (it imports nothing from solution/ or sim/).

    python3 other_analysis.py <night folder> <result json>

Differences from solution/reduce_night.py on purpose:
  * coupled amplitude from scipy.integrate.quad of the truncated Gaussian pupil field (not a fixed quadrature);
  * coupling and its image-motion averages taken at the middle of each channel (2 x 2 moments per telescope and
    channel); channel sums of the fringe terms in closed form (top-hat in wavenumber -> sinc), disc visibility on
    16 nodes;
  * image motion averaged in polar coordinates (Gauss-Laguerre in r^2, 12 angles);
  * zenith direction from the parallactic angle; refraction from the notes' law;
  * calibrator test on second differences in time of each calibrator's raw closure phases; transfer factor and
    closure offset as quadratics in hours from the middle of the night;
  * companion search on closure phases and squared visibilities together (linear templates, 1.5 mas grid, then
    0.15 mas), fit with trust-region least squares;
  * true sky: own nutation (four largest terms), Earth velocity from the change of the Sun's geocentric position,
    first-order aberration, apparent sidereal time; offsets fitted directly as barycentric J2000 offsets;
  * standard-air table and scans turned to vacuum with the notes' refractivity law.
"""
import csv
import datetime
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.laguerre import laggauss
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.special import j0, j1

MAS = math.pi / 648e6
PHI, LAMBDA_E = math.radians(34.2244), math.radians(-118.0570)
POS = {"S1": (0.0, 0.0, 0.0), "S2": (-5.748, -33.586, 0.637), "E1": (125.333, 305.928, -5.919),
       "E2": (70.396, 269.718, -2.856), "W1": (-175.073, 216.327, -10.786), "W2": (-69.086, 199.342, 0.471)}
TEL = list(POS)
STOP_R = 200.0
LK = 2.20

# ------------------------------------------------------------------- coupling: A(z) = int_0^1 exp(-q r^2) J0(z r) r dr
_Q = 1.2564312086261697
_Z = np.linspace(0, 14, 2801)
_A = np.array([quad(lambda r: math.exp(-_Q * r * r) * j0(z * r) * r, 0, 1, limit=200)[0] for z in _Z])
_A /= _A[0]


def amp(d_mas, s):
    return np.interp(math.pi * 1.0 * np.asarray(d_mas) * MAS * np.asarray(s), _Z, _A)


# ------------------------------------------------------------------- polar quadrature for a 2-D Gaussian
_LX, _LW = laggauss(8)
_ANG = np.arange(12) * 2 * math.pi / 12


def jitter_nodes(sigma):
    r = sigma * np.sqrt(2 * _LX)
    dx = (r[:, None] * np.cos(_ANG)[None, :]).ravel()
    dy = (r[:, None] * np.sin(_ANG)[None, :]).ravel()
    w = (_LW[:, None] * np.ones(12)[None, :] / 12).ravel()
    return dx, dy, w / w.sum()


# ------------------------------------------------------------------- astronomy
def gmst_rad(mjd):
    return math.radians((280.46061837 + 360.98564736629 * (mjd - 51544.5)) % 360)


def prec(mjd):
    T = (mjd - 51544.5) / 36525
    z1 = math.radians((2306.2181 * T + 0.30188 * T * T + 0.017998 * T ** 3) / 3600)
    z2 = math.radians((2306.2181 * T + 1.09468 * T * T + 0.018203 * T ** 3) / 3600)
    th = math.radians((2004.3109 * T - 0.42665 * T * T - 0.041833 * T ** 3) / 3600)
    cz1, sz1, cz2, sz2, ct, st = math.cos(z1), math.sin(z1), math.cos(z2), math.sin(z2), math.cos(th), math.sin(th)
    return np.array([[cz1 * ct * cz2 - sz1 * sz2, -sz1 * ct * cz2 - cz1 * sz2, -st * cz2],
                     [cz1 * ct * sz2 + sz1 * cz2, -sz1 * ct * sz2 + cz1 * cz2, -st * sz2],
                     [cz1 * st, -sz1 * st, ct]])


def vec(ra, de):
    return np.array([math.cos(de) * math.cos(ra), math.cos(de) * math.sin(ra), math.sin(de)])


def radec(v):
    return math.atan2(v[1], v[0]) % (2 * math.pi), math.asin(v[2] / np.linalg.norm(v))


def tangent(ra, de, v):
    """east, north (mas) of direction v in the tangent plane at (ra, de), gnomonic"""
    p = vec(ra, de)
    e = np.array([-math.sin(ra), math.cos(ra), 0])
    n = np.cross(p, e)
    return (v @ e) / (v @ p) / MAS, (v @ n) / (v @ p) / MAS


def from_tangent(ra, de, x, y):
    p = vec(ra, de)
    e = np.array([-math.sin(ra), math.cos(ra), 0])
    n = np.cross(p, e)
    v = p + x * MAS * e + y * MAS * n
    return v / np.linalg.norm(v)


def eq_baseline(t1, t2):
    de, dn, du = (POS[t2][i] - POS[t1][i] for i in range(3))
    return np.array([-math.sin(PHI) * dn + math.cos(PHI) * du, de, math.cos(PHI) * dn + math.sin(PHI) * du])


def uv_at(B, H, de):
    return (np.sin(H) * B[0] + np.cos(H) * B[1],
            -math.sin(de) * np.cos(H) * B[0] + math.sin(de) * np.sin(H) * B[1] + math.cos(de) * B[2])


def lift_dir(H, de):
    """zenith distance and the (east, north) unit vector towards the zenith, from the parallactic angle"""
    q = math.atan2(math.sin(H), math.tan(PHI) * math.cos(de) - math.sin(de) * math.cos(H))
    z = math.acos(math.sin(PHI) * math.sin(de) + math.cos(PHI) * math.cos(de) * math.cos(H))
    return z, np.array([math.sin(q), math.cos(q)])


def nm1(lam, p, t):
    return 1e-6 * p / t * (77.53 + 0.4395 / lam ** 2 + 0.0059 / lam ** 4)


def vac(lam_air):
    """standard-air wavelength (um) -> vacuum wavelength"""
    lam = lam_air
    for _ in range(4):
        lam = lam_air * (1 + nm1(lam, 1013.25, 288.15))
    return lam


def nut_small(mjd):
    """nutation from its four largest terms (IAU 1980), rad; mean obliquity"""
    T = (mjd + 0.0008 - 51544.5) / 36525
    om = math.radians(125.04452 - 1934.136261 * T)
    ls = math.radians(280.4665 + 36000.7698 * T)
    lm = math.radians(218.3165 + 481267.8813 * T)
    dpsi = (-17.20 * math.sin(om) - 1.32 * math.sin(2 * ls) - 0.23 * math.sin(2 * lm) + 0.21 * math.sin(2 * om)) / 3600
    deps = (9.20 * math.cos(om) + 0.57 * math.cos(2 * ls) + 0.10 * math.cos(2 * lm) - 0.09 * math.cos(2 * om)) / 3600
    eps0 = 23.4392911 - 0.0130042 * T
    return math.radians(dpsi), math.radians(deps), math.radians(eps0)


def sun_geo(mjd):
    """geocentric Sun, ecliptic of date (au), low precision"""
    T = (mjd - 51544.5) / 36525
    M = math.radians(357.52911 + 35999.05029 * T)
    lon = math.radians(280.46646 + 36000.76983 * T + (1.914602 - 0.004817 * T) * math.sin(M) + 0.019993 * math.sin(2 * M))
    e = 0.016708634 - 0.000042037 * T
    nu = M + math.radians((1.914602 - 0.004817 * T) * math.sin(M) + 0.019993 * math.sin(2 * M))
    r = 1.000001018 * (1 - e * e) / (1 + e * math.cos(nu))
    return np.array([r * math.cos(lon), r * math.sin(lon), 0.0])


def earth_beta(mjd):
    """Earth's velocity / c on the mean equator of date, from the change of the Sun's geocentric position"""
    h = 0.01
    v = -(sun_geo(mjd + h) - sun_geo(mjd - h)) / (2 * h)            # au / day
    _, _, eps0 = nut_small(mjd)
    v = np.array([v[0], v[1] * math.cos(eps0), v[1] * math.sin(eps0)])
    return v * 149597870.7e3 / 86400 / 299792458.0


def true_frame(mjd):
    """rotation J2000 -> true equator and equinox of date, local apparent sidereal time, observer velocity / c"""
    dpsi, deps, eps0 = nut_small(mjd)
    eps = eps0 + deps
    c, s_ = math.cos, math.sin

    def R1(a):
        return np.array([[1, 0, 0], [0, c(a), s_(a)], [0, -s_(a), c(a)]])

    def R3(a):
        return np.array([[c(a), s_(a), 0], [-s_(a), c(a), 0], [0, 0, 1]])
    Nm = R1(-eps) @ R3(-dpsi) @ R1(eps0)
    last = gmst_rad(mjd) + dpsi * c(eps) + LAMBDA_E
    beta = Nm @ earth_beta(mjd) + 465.1 * c(PHI) / 299792458.0 * np.array([-s_(last), c(last), 0])
    return Nm @ prec(mjd + 0.0008), last, beta


# ------------------------------------------------------------------- reading
def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def sexa(txt, hours):
    neg = txt.strip().startswith("-")
    a, b, c = (float(x) for x in re.findall(r"[0-9.]+", txt)[:3])
    v = (a + b / 60 + c / 3600) * (15 if hours else 1)
    return math.radians(-v if neg else v)


def target_file(path):
    info = {}
    for line in open(path):
        if ":" not in line or line.startswith("#"):
            continue
        k, v = (x.strip() for x in line.split(":", 1))
        if k == "target_ra_j2000":
            info["ra"] = sexa(v, True)
        elif k == "target_dec_j2000":
            info["dec"] = sexa(v, False)
        elif k == "companion_k_flux_ratio_percent":
            info["fk"] = float(v) / 100
        elif k.startswith("calibrator_"):
            d = {}
            for part in v.split(","):
                key, val = part.strip().split(" ", 1)
                d[key] = val
            info[k[-1]] = (float(d["limb_darkened_diameter_mas"]), float(d["linear_limb_darkening_h"]))
    return info


def log_refocus(path):
    out = []
    for line in open(path):
        if line.rstrip().endswith("spectrograph camera refocused"):
            stamp = datetime.datetime.strptime(line[:19], "%Y-%m-%d %H:%M:%S")
            out.append((stamp - datetime.datetime(1858, 11, 17)).total_seconds() / 86400)
    return sorted(out)


# ------------------------------------------------------------------- the reduction
class Night:
    def __init__(self, folder):
        f = Path(folder)
        self.info = target_file(f / "target.txt")
        self.V = rows(f / "vis2.csv")
        self.C = rows(f / "closure.csv")
        self.ref = log_refocus(f / "night_log.txt")
        self.W = np.array([[float(r["mjd"]), float(r["pressure_hpa"]), float(r["temperature_k"])] for r in rows(f / "weather.csv")])
        self.TT = {(r["star"], r["obs"], r["tel"]): float(r["tt_rms_mas"]) for r in rows(f / "ao_telemetry.csv")}
        scans = {}
        for r in rows(f / "wavecal.csv"):
            scans.setdefault(r["scan"], []).append((float(r["mjd"]), float(r["lambda_measured_um"]) / float(r["lambda_table_um"])))
        self.kscan = {}
        for sc in scans.values():
            seg = self.seg(sc[0][0])
            ratio = np.array([x[1] for x in sc])
            self.kscan.setdefault(seg, []).append(float(np.sum(ratio) / len(ratio)))
        self.nseg = len(self.ref) + 1
        self.free = [g for g in range(self.nseg) if g not in self.kscan]
        self.tmid = 0.5 * (min(float(r["mjd"]) for r in self.V) + max(float(r["mjd"]) for r in self.V))

    def seg(self, mjd):
        return sum(1 for t in self.ref if t <= mjd)

    def kfac(self, mjd, free):
        g = self.seg(mjd)
        if g in self.kscan:
            return float(np.mean(self.kscan[g]))
        return free[self.free.index(g)]

    def wx(self, mjd):
        return np.interp(mjd, self.W[:, 0], self.W[:, 1]), np.interp(mjd, self.W[:, 0], self.W[:, 2])


def ld_vis2(theta, u, b, s_lo, s_hi):
    s = s_lo + (np.arange(16) + 0.5) / 16 * (s_hi - s_lo)
    x = math.pi * theta * MAS * b * s
    v = ((1 - u) * j1(x) / x + u * (np.sin(x) / x ** 3 - np.cos(x) / x ** 2)) / ((1 - u) / 2 + u / 3)
    return float(np.mean(v)) ** 2


def calibrate(N, free):
    hrs = lambda m: (m - N.tmid) * 24
    # which calibrators are single: second differences in time of raw closure phases
    bad = {}
    for star in ("A", "B"):
        ser = {}
        for r in N.C:
            if r["star"] == star:
                ser.setdefault((r["tel1"], r["tel2"], r["tel3"], r["lambda_lo_um"]), []).append(
                    (float(r["mjd"]), float(r["closure_deg"]), float(r["closure_err_deg"])))
        chi, n = 0.0, 0
        for v in ser.values():
            v.sort()
            y = np.degrees(np.unwrap(np.radians([x[1] for x in v])))
            e = np.array([x[2] for x in v])
            d2 = y[2:] - 2 * y[1:-1] + y[:-2]
            chi += float(np.sum(d2 ** 2 / (e[2:] ** 2 + 4 * e[1:-1] ** 2 + e[:-2] ** 2)))
            n += len(d2)
        bad[star] = chi / n
    use = [s for s in ("A", "B") if bad[s] < 3.0]
    tv, tc = {}, {}
    for r in N.V:
        if r["star"] in use:
            k = N.kfac(float(r["mjd"]), free)
            th, u = N.info[r["star"]]
            b = math.hypot(float(r["u_m"]), float(r["v_m"]))
            mod = ld_vis2(th, u, b, 1e6 / vac(float(r["lambda_hi_um"]) * k), 1e6 / vac(float(r["lambda_lo_um"]) * k))
            tv.setdefault((r["tel1"], r["tel2"], r["lambda_lo_um"]), []).append(
                (hrs(float(r["mjd"])), float(r["vis2"]) / mod, float(r["vis2_err"]) / mod))
    for r in N.C:
        if r["star"] in use:
            tc.setdefault((r["tel1"], r["tel2"], r["tel3"], r["lambda_lo_um"]), []).append(
                (hrs(float(r["mjd"])), float(r["closure_deg"]), float(r["closure_err_deg"])))

    def qfit(pts, t, unwrap=False):
        pts = sorted(pts)
        x = np.array([p[0] for p in pts])
        y = np.array([p[1] for p in pts])
        if unwrap:
            y = np.degrees(np.unwrap(np.radians(y)))
        w = 1 / np.array([p[2] for p in pts])
        A = np.vstack([np.ones_like(x), x, x * x]).T
        c, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
        cov = np.linalg.inv((A * w[:, None]).T @ (A * w[:, None]))
        a = np.array([1, t, t * t])
        return float(a @ c), float(math.sqrt(a @ cov @ a))

    sv, sc = [], []
    for r in N.V:
        if r["star"] == "SCI":
            T, eT = qfit(tv[(r["tel1"], r["tel2"], r["lambda_lo_um"])], hrs(float(r["mjd"])))
            y = float(r["vis2"]) / T
            sv.append((r, y, math.hypot(float(r["vis2_err"]) / T, y * eT / T)))
    for r in N.C:
        if r["star"] == "SCI":
            O, eO = qfit(tc[(r["tel1"], r["tel2"], r["tel3"], r["lambda_lo_um"])], hrs(float(r["mjd"])), True)
            sc.append((r, (float(r["closure_deg"]) - O + 180) % 360 - 180, math.hypot(float(r["closure_err_deg"]), eO)))
    return sv, sc, bad, use


class Forward:
    def __init__(self, N, sv, sc, nt=10):
        self.N = N
        ra, de = radec(prec(N.tmid) @ vec(N.info["ra"], N.info["dec"]))
        self.obs = {}
        tn, tw = np.polynomial.legendre.leggauss(nt)
        for r, _, _ in sv:
            self.obs.setdefault(r["obs"], {"mjd": float(r["mjd"]), "T": float(r["int_time_s"])})
        for o, d in self.obs.items():
            m = d["mjd"]
            R, last, beta = true_frame(m)
            pg = R @ vec(N.info["ra"], N.info["dec"]) + beta
            pg /= np.linalg.norm(pg)
            rao, deo = radec(pg)
            d["R"], d["beta"], d["pg"] = R, beta, pg
            lasts = last + tn / 2 * d["T"] * 7.292115855e-5
            d["tw"] = tw / 2
            d["B"] = {}
            for i, a in enumerate(TEL):
                for b in TEL[i + 1:]:
                    X, Y, Z = eq_baseline(a, b)
                    d["B"][(a, b)] = np.c_[X * np.cos(lasts) - Y * np.sin(lasts), X * np.sin(lasts) + Y * np.cos(lasts),
                                           np.full_like(lasts, Z)]
            d["z"], d["up"] = lift_dir(last - rao, deo)
            d["wx"] = N.wx(m)
            d["sig"] = {t: N.TT[("SCI", o, t)] for t in TEL}
        self.sv, self.sc = sv, sc

    def mom(self, d, e, n, s, tel):
        p, t = d["wx"]
        lift = (nm1(1e6 / s, p, t) - nm1(LK, p, t)) * math.tan(d["z"]) / MAS
        w = self.N.info["fk"] / (1 + self.N.info["fk"]) if math.hypot(e, n) <= STOP_R else 0.0
        g = lift * d["up"] - w * np.array([e, n])
        c = g + np.array([e, n])
        dx, dy, wq = jitter_nodes(d["sig"][tel])
        ag = amp(np.hypot(g[0] + dx, g[1] + dy), s)
        ac = amp(np.hypot(c[0] + dx, c[1] + dy), s)
        return np.array([[wq @ (ag * ag), wq @ (ag * ac)], [wq @ (ag * ac), wq @ (ac * ac)]])

    def terms(self, d, a, b, theta, pc, f, s_lo, s_hi):
        """channel sums of giant and companion terms at the time nodes, baseline a-b (true sky)"""
        B = d["B"][(a, b)] if (a, b) in d["B"] else -d["B"][(b, a)]
        bb = np.sqrt(np.maximum((B * B).sum(1) - (B @ d["pg"]) ** 2, 0))
        s = s_lo + (np.arange(16) + 0.5) / 16 * (s_hi - s_lo)
        x = math.pi * theta * (1 - d["beta"] @ d["pg"]) * MAS * bb[:, None] * s[None, :]
        giant = np.mean(2 * j1(x) / x, 1)
        ph = B @ (pc - d["pg"])
        comp = f * np.exp(-2j * math.pi * 0.5 * (s_lo + s_hi) * ph) * np.sinc((s_hi - s_lo) * ph)
        return np.vstack([giant + 0j, comp])                        # (2, T)

    def model(self, theta, e2000, n2000, f, free):
        N = self.N
        cache = {}
        out_v, out_c = [], []
        c0 = from_tangent(N.info["ra"], N.info["dec"], e2000, n2000)
        for d in self.obs.values():
            pc = d["R"] @ c0 + d["beta"]
            d["pc"] = pc / np.linalg.norm(pc)
            ra_, de_ = radec(d["pg"])
            d["off"] = tangent(ra_, de_, d["pc"])
        for r, _, _ in self.sv:
            d = self.obs[r["obs"]]
            k = N.kfac(d["mjd"], free)
            s_lo, s_hi = 1e6 / vac(float(r["lambda_hi_um"]) * k), 1e6 / vac(float(r["lambda_lo_um"]) * k)
            key = (r["obs"], r["lambda_lo_um"])
            if key not in cache:
                cache[key] = {t: self.mom(d, *d["off"], 0.5 * (s_lo + s_hi), t) for t in TEL}
            M = cache[key]
            a, b = r["tel1"], r["tel2"]
            X = self.terms(d, a, b, theta, d["pc"], f, s_lo, s_hi)
            P = M[a] * M[b]
            num = np.einsum("kl,kt,lt->t", P, X, np.conj(X)).real @ d["tw"]
            fa = M[a][0, 0] + f * M[a][1, 1]
            fb = M[b][0, 0] + f * M[b][1, 1]
            out_v.append(num / (fa * fb))
        for r, _, _ in self.sc:
            d = self.obs[r["obs"]]
            k = N.kfac(d["mjd"], free)
            s_lo, s_hi = 1e6 / vac(float(r["lambda_hi_um"]) * k), 1e6 / vac(float(r["lambda_lo_um"]) * k)
            key = (r["obs"], r["lambda_lo_um"])
            if key not in cache:
                cache[key] = {t: self.mom(d, *d["off"], 0.5 * (s_lo + s_hi), t) for t in TEL}
            M = cache[key]
            a, b, c = r["tel1"], r["tel2"], r["tel3"]
            X12 = self.terms(d, a, b, theta, d["pc"], f, s_lo, s_hi)
            X23 = self.terms(d, b, c, theta, d["pc"], f, s_lo, s_hi)
            X31 = self.terms(d, c, a, theta, d["pc"], f, s_lo, s_hi)
            tp = np.einsum("ik,ij,jk,it,jt,kt->t", M[a], M[b], M[c], X12, X23, X31) @ d["tw"]
            out_c.append(math.degrees(np.angle(tp)))
        return np.array(out_v), np.array(out_c)


def search(F, theta):
    """chi-square gain of a companion term (linear in its strength) over a grid, closure phases + vis2"""
    sv, sc = F.sv, F.sc
    N = F.N
    U = np.array([float(r["u_m"]) for r, _, _ in sv])
    Vv = np.array([float(r["v_m"]) for r, _, _ in sv])
    k = np.array([N.kfac(float(r["mjd"]), [N.kscan[min(N.kscan)][0]] * len(N.free)) for r, _, _ in sv])
    s = 2e6 / ((np.array([float(r["lambda_lo_um"]) for r, _, _ in sv]) + np.array([float(r["lambda_hi_um"]) for r, _, _ in sv])) * k)
    y = np.array([x[1] for x in sv])
    e = np.array([x[2] for x in sv])
    b = np.hypot(U, Vv)
    xx = math.pi * theta * MAS * b * s
    D = 2 * j1(xx) / xx
    scale = np.sum(y * D * D / e ** 2) / np.sum(D ** 4 / e ** 2)
    res = y / scale - D * D
    short = b < 125
    U1 = np.array([float(r["u1_m"]) for r, _, _ in sc])
    V1 = np.array([float(r["v1_m"]) for r, _, _ in sc])
    U2 = np.array([float(r["u2_m"]) for r, _, _ in sc])
    V2 = np.array([float(r["v2_m"]) for r, _, _ in sc])
    kc = np.array([N.kfac(float(r["mjd"]), [N.kscan[min(N.kscan)][0]] * len(N.free)) for r, _, _ in sc])
    scl = 2e6 / ((np.array([float(r["lambda_lo_um"]) for r, _, _ in sc]) + np.array([float(r["lambda_hi_um"]) for r, _, _ in sc])) * kc)
    yc = np.radians(np.array([x[1] for x in sc]))
    ec = np.radians(np.array([x[2] for x in sc]))
    legs = [(U1, V1), (U2, V2), (-U1 - U2, -V1 - V2)]
    Dl = []
    for uu, vv in legs:
        xq = math.pi * theta * MAS * np.hypot(uu, vv) * scl
        Dl.append(2 * j1(xq) / xq)
    sel = (np.abs(Dl[0]) > 0.15) & (np.abs(Dl[1]) > 0.15) & (np.abs(Dl[2]) > 0.15) & \
          (np.maximum.reduce([np.hypot(*l) for l in legs]) < 200)

    def gain(pts, use_all):
        g = np.zeros(len(pts))
        for i0 in range(0, len(pts), 1500):
            p = pts[i0:i0 + 1500]
            m = slice(None) if use_all else short
            ph = (np.outer(p[:, 0], U[m]) + np.outer(p[:, 1], Vv[m])) * MAS * s[m]
            t = 2 * D[m] * np.cos(2 * math.pi * ph)
            num = (t * res[m] / e[m] ** 2).sum(1)
            den = (t * t / e[m] ** 2).sum(1)
            if not use_all:                      # coarse grid: short-baseline squared visibilities only
                g[i0:i0 + 1500] = np.where(num > 0, num * num / den, 0)
                continue
            tc = 0
            for (uu, vv), dl in zip(legs, Dl):
                phc = (np.outer(p[:, 0], uu[sel]) + np.outer(p[:, 1], vv[sel])) * MAS * scl[sel]
                tc = tc - np.sin(2 * math.pi * phc) / dl[sel]
            num += (tc * yc[sel] / ec[sel] ** 2).sum(1)
            den += (tc * tc / ec[sel] ** 2).sum(1)
            g[i0:i0 + 1500] = np.where(num > 0, num * num / den, 0)
        return g

    ax = np.arange(-330, 330.1, 1.0)
    gx, gy = np.meshgrid(ax, ax)
    rr = np.hypot(gx, gy)
    pts = np.c_[gx[(rr > 8) & (rr < 330)], gy[(rr > 8) & (rr < 330)]]
    g = gain(pts, False)
    best = []
    for i in np.argsort(g)[::-1]:
        if all(np.hypot(*(pts[i] - q)) > 4 and np.hypot(*(pts[i] + q)) > 4 for q in best):
            best.append(pts[i])
        if len(best) == 8:
            break
    loc = np.arange(-2.25, 2.26, 0.15)
    lx, ly = np.meshgrid(loc, loc)
    lp = np.c_[lx.ravel(), ly.ravel()]
    out = []
    for q0 in best:
        for q in (q0, -q0):                      # the squared visibilities cannot tell a place from its mirror
            gg = gain(q + lp, True)
            j = int(np.argmax(gg))
            out.append((float(gg[j]), q + lp[j]))
    out.sort(key=lambda t: -t[0])
    return out[:3]


def reduce(folder):
    N = Night(folder)
    free = [N.kscan[max(g for g in N.kscan if g < fg)][0] if any(g < fg for g in N.kscan) else N.kscan[min(N.kscan)][0]
            for fg in N.free]
    p = None
    for _ in range(2 if N.free else 1):
        sv, sc, bad, use = calibrate(N, free)
        F = Forward(N, sv, sc)
        y2 = np.array([x[1] for x in sv]); e2 = np.array([x[2] for x in sv])
        y3 = np.array([x[1] for x in sc]); e3 = np.array([x[2] for x in sc])

        def res(q):
            m2, m3 = F.model(q[0], q[1], q[2], q[3], list(q[4:]))
            return np.r_[(y2 - m2) / e2, ((y3 - m3 + 180) % 360 - 180) / e3]

        if p is None:
            theta = 1.0
            b = np.array([math.hypot(float(r["u_m"]), float(r["v_m"])) for r, _, _ in sv])
            s = np.array([2e6 / (float(r["lambda_lo_um"]) + float(r["lambda_hi_um"])) for r, _, _ in sv])
            theta = least_squares(lambda t: (y2 - t[1] * (2 * j1(math.pi * t[0] * MAS * b * s) / (math.pi * t[0] * MAS * b * s)) ** 2) / e2,
                                  [1.0, 1.0]).x[0]
            cands = search(F, theta)
            # stage 1: only squared visibilities on baselines under 150 m (both mirror images, three strengths)
            short = [x for x in sv if math.hypot(float(x[0]["u_m"]), float(x[0]["v_m"])) < 150]
            Fs = Forward(N, short, [])
            ys = np.array([x[1] for x in short]); es = np.array([x[2] for x in short])
            trials = []
            m0 = float(np.mean([d["mjd"] for d in F.obs.values()]))
            P0 = prec(m0)
            ra0, de0 = radec(P0 @ vec(N.info["ra"], N.info["dec"]))
            for g, q in cands:
                # the search uses the files' (u, v): mean equator of the date -> J2000
                q = np.array(tangent(N.info["ra"], N.info["dec"], P0.T @ from_tangent(ra0, de0, q[0], q[1])))
                for sgn in (1,):
                    best1 = None
                    for fs in (0.01, 0.03, 0.08):
                        r1 = least_squares(lambda z: (ys - Fs.model(z[0], z[1], z[2], z[3], free)[0]) / es,
                                           [theta, sgn * q[0], sgn * q[1], fs], method="trf",
                                           x_scale=[1e-3, 1e-2, 1e-2, 1e-4], max_nfev=30)
                        if best1 is None or r1.cost < best1.cost:
                            best1 = r1
                    x0 = np.r_[best1.x, free]
                    trials.append((float(np.sum(res(x0) ** 2)), x0))
            trials.sort(key=lambda t: t[0])
            fits = [least_squares(res, x0, method="trf", x_scale=np.r_[1e-3, 1e-2, 1e-2, 1e-4, [1e-4] * len(free)])
                    for _, x0 in trials[:2]]
            p = min(fits, key=lambda r: r.cost).x
        else:
            p = least_squares(res, p, method="trf", x_scale=np.r_[1e-3, 1e-2, 1e-2, 1e-4, [1e-4] * len(free)]).x
        free = list(p[4:])
        chi = float(np.sum(res(p) ** 2))
    x, y = p[1], p[2]                       # fitted directly as the barycentric offset on J2000
    return {"diameter_mas": float(p[0]), "flux_ratio_percent": float(100 * p[3]), "separation_mas": float(math.hypot(x, y)),
            "position_angle_deg": float(math.degrees(math.atan2(x, y)) % 360),
            "fit": {"chi_square": chi, "measurements": len(y2) + len(y3), "calibrators_used": use,
                    "free_wavelength_factors": free, "closure_second_difference_chi2": bad}}


if __name__ == "__main__":
    out = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("diameter_mas", "flux_ratio_percent", "separation_mas", "position_angle_deg")}))
