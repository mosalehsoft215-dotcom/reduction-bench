#!/usr/bin/env python3
"""A second reduction of a DLTS diode, written apart from solution/trap_reduce.py.

    python second_reduction.py <diode folder> <answer.json>

Differences in method from the reference:
  - the junction capacitance comes from the complex impedance: Z = 1/(G + jwC) of the meter reading, and the
    contact resistance R = Re Z is taken off before the junction admittance is formed again;
  - the donor profile is read point by point from the C-V curve, N(W) = 2 / (q eps A^2 d(1/C^2)/dV) on a
    smoothed derivative, and the profile below the first point is held at that point's value;
  - rates come from curve_fit of C(t) directly through the depletion relation solved for W(t) (no 1/C^2 form);
  - the kind of centre is chosen from the bias series by the rise of the single-exponential rate with bias;
  - the Arrhenius line is unweighted; N_T is found by root-finding the bias relation for the t = 0 width.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.signal import savgol_filter

q, kB, eps = 1.602176634e-19, 8.617333262e-5, 11.7 * 8.8541878128e-14
GL_X, GL_W = np.polynomial.legendre.leggauss(64)


def Nc(T):
    return 2.86e19 * (T / 300.0) ** 1.5


def Vth(T):
    return 2.30e7 * (T / 300.0) ** 0.5


def load(path):
    raw = np.genfromtxt(path, delimiter=",", names=True)
    return {k: np.atleast_1d(raw[k]) for k in raw.dtype.names}


def split(tab):
    cut = np.flatnonzero(np.diff(tab["t_s"]) < 0) + 1
    edges = np.concatenate([[0], cut, [len(tab["t_s"])]])
    return [{k: v[a:b] for k, v in tab.items()} for a, b in zip(edges[:-1], edges[1:])]


class Diode:
    def __init__(self, folder):
        f = Path(folder)
        self.sheet = json.loads((f / "diode.json").read_text())
        self.A = self.sheet["area_cm2"]
        self.w = 2 * math.pi * self.sheet["meter"]["frequency_hz"]
        self.vr = self.sheet["scan"]["reverse_bias_v"]
        self.vp = self.sheet["scan"]["pulse_v"]
        self.cv = load(f / "cv_300k.csv")
        self.scan = split(load(f / "scan.csv"))
        self.bias = split(load(f / "bias.csv"))
        self.cap = split(load(f / "capture.csv"))
        self.chk = load(f / "checks.csv")

    def cj(self, c_pf, g_us):
        y = np.asarray(g_us) * 1e-6 + 1j * self.w * np.asarray(c_pf) * 1e-12
        z = 1.0 / y
        yj = 1.0 / (1j * z.imag)
        return yj.imag / self.w

    def ok(self, rec):
        cm = rec["C_pF"] * 1e-12
        return bool(np.all(cm > 0) and np.max(rec["G_uS"] * 1e-6 / (self.w * cm)) < 1.5)


class Profile:
    def __init__(self, W, N, phi, level=0.45, NT=0.0, donor=False):
        self.W, self.N, self.phi, self.level, self.NT, self.donor = W, N, phi, level, NT, donor

    def nd(self, x):
        return np.interp(x, self.W, self.N)

    def nfree(self, x=3e-4):
        return self.nd(x) - (0 if self.donor else self.NT)

    def ef(self, T, x=3e-4):
        return kB * T * math.log(Nc(T) / self.nfree(x))

    def lam(self, T, x):
        d = self.level - self.ef(T, x)
        return math.sqrt(2 * eps * d / (q * self.nd(x))) if d > 0 else 0.0

    def mom(self, a, b, extra=0.0):
        """integral_a^b x (N_d(x) + extra) dx by Gauss-Legendre."""
        if b <= a:
            return 0.0
        x = 0.5 * (b - a) * GL_X + 0.5 * (a + b)
        return 0.5 * (b - a) * float(np.sum(GL_W * x * (self.nd(x) + extra)))

    def pot_steady(self, T, W):
        x1 = max(W - self.lam(T, W), 0.0)
        e = self.NT if self.donor else 0.0
        fl = 0.0 if self.donor else -self.NT
        return q / eps * (self.mom(0, x1, e) + self.mom(x1, W, fl))

    def W_at(self, T, V):
        rhs = self.phi - self.ef(T) - V - kB * T
        return brentq(lambda W: self.pot_steady(T, W) - rhs, 1e-7, 6e-3, xtol=1e-12)


def profile_from_cv(d, level=0.45, NT=0.0, donor=False):
    v = d.cv["bias_v"]
    C = d.cj(d.cv["C_pF"], d.cv["G_uS"])
    W = eps * d.A / C
    inv = 1.0 / C ** 2
    slope = savgol_filter(inv, 9, 2, deriv=1, delta=v[1] - v[0])
    N = -2.0 / (q * eps * d.A ** 2 * slope)
    order = np.argsort(W)
    W, N, v = W[order], N[order], v[order]
    # the lambda layer makes the C-V profile read N_d - (lambda / W) N_T for acceptor-like traps
    if NT > 0 and not donor:
        P0 = Profile(W, N, 0.8, level, NT, donor)
        N = N + np.array([P0.lam(300.0, w) / w for w in W]) * NT
    # barrier from the zero-bias point
    P = Profile(W, N, 0.8, level, NT, donor)
    phi = P.pot_steady(300.0, W[np.argmin(abs(v))]) + v[np.argmin(abs(v))] + kB * 300 + P.ef(300.0)
    return Profile(W, N, phi, level, NT, donor)


def window(P, T, Winf, Vp):
    Wp = P.W_at(T, Vp)
    return max(Wp - P.lam(T, Wp), 0.0), Winf - P.lam(T, Winf)


def pf_mult(P, T, x, x2, Winf):
    if not P.donor:
        return np.ones_like(x)
    e = P.NT
    F = np.array([q / eps * (quad_n(P, xi, x2, e) + quad_n(P, x2, Winf, 0.0)) for xi in x])
    return np.exp(np.sqrt(q * F / (math.pi * eps)) / (kB * T))


def quad_n(P, a, b, extra):
    if b <= a:
        return 0.0
    x = 0.5 * (b - a) * GL_X + 0.5 * (a + b)
    return 0.5 * (b - a) * float(np.sum(GL_W * (P.nd(x) + extra)))


def transient_model(P, T, Vr, Vp, d):
    """Returns f(t, e0, Winf_um, NT_eff) giving C(t) for the window filled at t = 0."""
    def model(t, le0, Winf_um, s):
        Winf = Winf_um * 1e-4
        x1, x2 = window(P, T, Winf, Vp)
        xg = 0.5 * (x2 - x1) * GL_X + 0.5 * (x1 + x2)
        m = pf_mult(P, T, xg, x2, Winf)
        occ = np.exp(-np.outer(t, 10 ** le0 * m))                # (nt, ng)
        S = 0.5 * (x2 - x1) * (occ * (GL_W * xg)[None, :]).sum(axis=1)   # integral x f dx
        S0 = (x2 * x2 - x1 * x1) / 2.0
        # W(t)^2 - Winf^2 = 2 s S(t) / N_b with s the trap density scale, solved with the local profile
        Nb = P.nd(Winf) - (0 if P.donor else s * 1e12)
        W = np.sqrt(Winf ** 2 + 2.0 * s * 1e12 * S / Nb)
        return eps * d.A / W * 1e12
    return model


def fit_rec(d, P, rec, Vr):
    keep = rec["t_s"] >= 6e-5
    t = rec["t_s"][keep]
    C = d.cj(rec["C_pF"][keep], rec["G_uS"][keep]) * 1e12
    T = float(rec["T_K"].mean())
    m = transient_model(P, T, Vr, d.vp, d)
    Winf0 = eps * d.A / (C[-1] * 1e-12) * 1e4
    best = None
    smax = 0.9 * float(P.nd(Winf0 * 1e-4)) / 1e12
    for le in np.arange(-0.5, 5.0, 0.5):
        try:
            p, cov = curve_fit(m, t, C, p0=[le, Winf0, min(10.0, 0.5 * smax)], maxfev=400,
                               bounds=([-3.0, 0.5 * Winf0, 0.0], [7.0, 2.0 * Winf0, smax]))
        except Exception:
            continue
        r = float(np.sum((m(t, *p) - C) ** 2))
        if best is None or r < best[1]:
            best = (p, r, cov)
    if best is None:
        return None
    p, r, cov = best
    amp = abs(m(np.array([0.0]), *p)[0] - m(np.array([1e3]), *p)[0])
    return {"T": T, "e0": 10 ** p[0], "Winf": p[1] * 1e-4, "s": p[2], "amp": amp}


def sample_t_map(d, P):
    """Sample temperature from the absolute forward voltage at 10 uA through the thermionic law. Ideality (from the
    two-current difference) and barrier are fixed near room temperature; the ohmic drop uses the meter's R."""
    Ts = d.chk["T_K"]
    R = []
    for rec in d.scan:
        z = 1.0 / (rec["G_uS"][-30:] * 1e-6 + 1j * d.w * rec["C_pF"][-30:] * 1e-12)
        R.append(float(np.mean(z.real)))
    R = np.interp(Ts, [float(r["T_K"].mean()) for r in d.scan], R)
    warm = (Ts > 285) & (Ts < 305)
    n = float(np.mean((d.chk["VF_100uA_v"] - d.chk["VF_10uA_v"] - 9e-5 * R)[warm] / (kB * Ts[warm] * math.log(10))))
    # barrier from the same check near room temperature, where the sample cannot sit away from the stage
    vj_w = (d.chk["VF_10uA_v"] - 1e-5 * R)[warm]
    phi = float(np.mean(vj_w / n - kB * Ts[warm] * np.log(1e-5 / (d.A * 112.0 * Ts[warm] ** 2))))
    out = []
    for T0, v, r in zip(Ts, d.chk["VF_10uA_v"], R):
        vj = v - 1e-5 * r
        f = lambda T: n * kB * T * math.log(1e-5 / (d.A * 112.0 * T * T * math.exp(-phi / (kB * T))) + 1.0) - vj
        try:
            out.append(brentq(f, 0.7 * T0, 1.3 * T0) - T0)
        except ValueError:
            out.append(np.nan)
    out = np.array(out)
    ok = np.isfinite(out) & (R * 9e-5 < 3e-3) & (Ts < 305)
    if ok.sum() < 20:
        return lambda t: t
    c = np.polyfit(Ts[ok], out[ok], 4)
    if np.max(abs(np.polyval(c, Ts[ok]))) < 1.0:
        return lambda t: t
    return lambda t: t + np.polyval(c, t)


def reduce(folder):
    d = Diode(folder)
    tmap = sample_t_map(d, profile_from_cv(d))
    for rec in d.scan + d.bias + d.cap:
        rec["T_K"] = tmap(rec["T_K"])
    d.scan = [r for r in d.scan if d.ok(r)]
    d.bias = [r for r in d.bias if d.ok(r)]
    d.cap = [r for r in d.cap if d.ok(r)]
    level, NT, donor = 0.45, 0.0, False
    P = profile_from_cv(d)
    # kind of centre: does the single-rate reading rise with reverse bias?
    plain = Profile(P.W, P.N, P.phi, level, NT, False)
    rates = {}
    for rec in d.bias:
        fr = fit_rec(d, plain, rec, float(rec["reverse_bias_v"][0]))
        if not fr:
            continue
        rates.setdefault(round(fr["T"] / 2.5), {})[float(rec["reverse_bias_v"][0])] = fr["e0"]
    rise = np.median([math.log(v[-10.0] / v[-2.0]) for v in rates.values() if -10.0 in v and -2.0 in v])
    donor = rise > 0.15
    eb = 0.0
    inside = []
    for rec in d.scan:
        fr = fit_rec(d, plain, rec, d.vr)
        if fr and 1.0 < fr["e0"] < 6000 and fr["amp"] > 0.03:
            inside.append(rec)
    for _ in range(3):
        P = profile_from_cv(d, level, NT, donor)
        pts = []
        for rec in inside:
            fr = fit_rec(d, P, rec, d.vr)
            if fr and 1.5 < fr["e0"] < 4000 and fr["amp"] > 0.03:
                pts.append(fr)
        T = np.array([p["T"] for p in pts])
        y = np.log(np.array([p["e0"] for p in pts]) / T ** 2)
        b, a = np.polyfit(1 / (kB * T), y, 1)
        ea = -b
        sig_a = math.exp(a) / (Vth(300) * Nc(300) / 300 ** 2)
        # capture series
        caps = {}
        for rec in d.cap:
            Tc = round(float(rec["T_K"].mean()) / 2.5) * 2.5
            fr = fit_rec(d, P, rec, d.vr)
            if fr:
                caps.setdefault(Tc, []).append((float(rec["pulse_width_s"][0]), fr["s"]))
        sig = {}
        for Tc, lst in caps.items():
            tp, s = np.array(lst).T
            frac = s / np.median(s[-3:])
            if frac[0] > 0.85:
                continue
            f = lambda tp, lc, a: a * (1 - np.exp(-10 ** lc * tp))
            p, cov = curve_fit(f, tp, s, p0=[math.log10(1 / tp[np.argmin(abs(frac - 0.6))]), np.median(s[-3:])])
            sig[Tc] = (10 ** p[0] / (Vth(Tc) * P.nfree()), math.sqrt(cov[0, 0]) * math.log(10))
        eb = 0.0
        if len(sig) == 2:
            (t1, (s1, u1)), (t2, (s2, u2)) = sorted(sig.items())
            est = kB * math.log(s2 / s1) / (1 / t1 - 1 / t2)
            if est > 3 * kB * math.hypot(u1, u2) / (1 / t1 - 1 / t2):
                eb = est
        level = ea - eb
        NT = float(np.median([p["s"] for p in pts if 5 < p["e0"] < 1500])) * 1e12
    return {"level_ev": round(level, 5), "density_cm3": float(f"{NT:.5e}"),
            "sigma_300k_cm2": float(f"{sig_a * math.exp(-eb / (kB * 300)):.5e}")}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce(sys.argv[1])) + "\n")
