#!/usr/bin/env python3
"""Reduce one DLTS diode folder to the electron trap's zero-field level, density and 300 K capture cross-section.

    python3 trap_reduce.py <diode folder> <answer.json>

Order of work:
  1. every C_pF / G_uS pair is turned into the junction capacitance: the meter sees the junction in series with the
     contact and substrate resistance, so C = C_m (1 + (G_m / (w C_m))^2);
  2. the 300 K C-V curve gives the donor profile N_d(x) = N0 (1 + g x) and the built-in voltage, through the
     depletion relation with the trap's lambda layer;
  3. each transient is read as 1/C^2 against time. In the depletion approximation W^2 - W_inf^2 is proportional to
     int x f(x, t) dx over the emitting window, so a centre with one rate gives an exact exponential in 1/C^2. A
     centre charged after emission (Poole-Frenkel) has a rate that follows the local field; then the shape is the
     window integral of exp(-e0 exp(gamma(x)) t). Which of the two the centre is follows from the bias series:
     only the field-dependent description makes the three reverse biases agree;
  4. ln(e0 / T^2) against 1/T gives the activation energy and the prefactor sigma_a;
  5. the capture series at two temperatures gives sigma_n(T) wherever the filling is resolved, and with it the
     capture barrier E_b; the level is E_a - E_b and sigma(300 K) = sigma_a exp(-E_b / k 300);
  6. N_T from the amplitude of 1/C^2 at t = 0 through the exact charge balance of the window and the layer beyond.
Steps 2-6 are repeated three times, since lambda needs the level and the C-V relation needs N_T.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

Q = 1.602176634e-19
KB = 8.617333262e-5
EPS = 11.7 * 8.8541878128e-14
T_SKIP = 6.0e-5           # meter recovery after the pulse


def nc(T):
    return 2.86e19 * (T / 300.0) ** 1.5


def vth(T):
    return 2.30e7 * np.sqrt(T / 300.0)


def read_table(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    data = np.array([[float(v) for v in r] for r in rows[1:] if r], dtype=float)
    return {h: data[:, i] for i, h in enumerate(head)}


def groups(table, keys):
    """Split a long table into records that share the values of `keys` (consecutive rows)."""
    out = []
    n = len(table["t_s"])
    start = 0
    for i in range(1, n + 1):
        if i == n or table["t_s"][i] < table["t_s"][i - 1]:
            sl = slice(start, i)
            rec = {k: table[k][sl] for k in table}
            out.append(rec)
            start = i
    return out


class Bench:
    def __init__(self, folder):
        folder = Path(folder)
        self.sheet = json.loads((folder / "diode.json").read_text())
        self.area = float(self.sheet["area_cm2"])
        self.w = 2.0 * math.pi * float(self.sheet["meter"]["frequency_hz"])
        self.vr = float(self.sheet["scan"]["reverse_bias_v"])
        self.vp = float(self.sheet["scan"]["pulse_v"])
        self.period = float(self.sheet["scan"]["period_s"])
        self.cv = read_table(folder / "cv_300k.csv")
        self.scan = groups(read_table(folder / "scan.csv"), ("T_K",))
        self.bias = groups(read_table(folder / "bias.csv"), ("T_K", "reverse_bias_v"))
        self.capture = groups(read_table(folder / "capture.csv"), ("T_K", "pulse_width_s"))
        self.checks = read_table(folder / "checks.csv")

    def sample_temperatures(self):
        """Sample temperature at each scan set point from the two-current forward check.

        V(100 uA) - V(10 uA) = n (kT/q) ln 10 + R_s * 90 uA, with R_s the series resistance the meter sees
        (Re of 1/(G + jwC) late in the transient). The ideality n is fixed at the 300 K set point, where the
        sample is at the stage temperature whatever its mounting, since the stage then matches the room."""
        stage = self.checks["T_K"]
        rs = []
        for rec in self.scan:
            y = rec["G_uS"][-20:] * 1e-6 + 1j * self.w * rec["C_pF"][-20:] * 1e-12
            rs.append(float(np.median((1.0 / y).real)))
        rs = np.interp(stage, [float(np.mean(r["T_K"])) for r in self.scan], rs)
        phi_b = first_cv(self).phi_b
        i_sat = self.area * 112.0 * stage ** 2 * np.exp(-phi_b / (KB * stage))
        ratio = np.log((1e-4 / i_sat + 1.0) / (1e-5 / i_sat + 1.0))
        u = (self.checks["VF_100uA_v"] - self.checks["VF_10uA_v"] - rs * 9e-5) / (KB * ratio)
        good = (rs * 9e-5 < 3e-3) & (stage < 305.0)          # the ohmic part must stay a small correction
        if good.sum() < 20:
            return lambda t: t
        # n T_sample with T_sample = T_stage + d0 (1 - (T_stage / 300)^4): the radiative load of the 300 K window
        # on a sample that is not clamped to the copper; linear in n and n d0
        g = 1.0 - (stage / 300.0) ** 4
        A = np.column_stack([stage, g])[good]
        coef, res, *_ = np.linalg.lstsq(A, u[good], rcond=None)
        resid = u[good] - A @ coef
        s2 = float(resid @ resid) / max(good.sum() - 2, 1)
        cov = np.linalg.inv(A.T @ A) * s2
        n, nd0 = coef
        d0 = nd0 / n
        d0_unc = math.sqrt(cov[1, 1]) / n
        if abs(d0) < 4.0 * d0_unc or abs(d0) < 0.3:          # mounted on the copper: the stage reading stands
            return lambda t: t
        return lambda t: t + d0 * (1.0 - (np.asarray(t) / 300.0) ** 4)

    def junction_c(self, c_pf, g_us):
        cm = np.asarray(c_pf) * 1e-12
        gm = np.asarray(g_us) * 1e-6
        return cm * (1.0 + (gm / (self.w * cm)) ** 2)


class Layer:
    """Donor profile, built-in voltage and the trap, with the depletion relations of the bench guide.

    The donor density is piecewise linear between knots (depth in cm) and constant beyond the end knots."""

    GRID = np.linspace(0.0, 6e-3, 30001)

    def __init__(self, kx, kn, phi_b, level, nt, donor):
        self.kx, self.kn = np.asarray(kx, float), np.asarray(kn, float)
        self.phi_b, self.level, self.nt, self.donor = phi_b, level, nt, donor
        x = self.GRID
        n = np.interp(x, self.kx, self.kn)
        dx = np.diff(x)
        self._i0 = np.concatenate([[0.0], np.cumsum(0.5 * (n[1:] + n[:-1]) * dx)])
        xn = x * n
        self._i1 = np.concatenate([[0.0], np.cumsum(0.5 * (xn[1:] + xn[:-1]) * dx)])

    def nd(self, x):
        return np.interp(x, self.kx, self.kn)

    def n_free(self, x):
        return self.nd(x) - (0.0 if self.donor else self.nt)

    def ecf(self, T, x=3e-4):
        return KB * T * math.log(nc(T) / self.n_free(x))

    def vbi(self, T):
        return self.phi_b - self.ecf(T)

    def lam(self, T, x):
        de = self.level - self.ecf(T, x)
        return math.sqrt(2.0 * EPS * max(de, 1e-4) / (Q * self.nd(x)))

    def ixd(self, a, b):
        """int_a^b x N_d(x) dx."""
        return np.interp(b, self.GRID, self._i1) - np.interp(a, self.GRID, self._i1)

    def idd(self, a, b):
        return np.interp(b, self.GRID, self._i0) - np.interp(a, self.GRID, self._i0)

    def potential(self, T, w):
        """Steady state: traps empty below the lambda point, filled beyond."""
        x1 = max(w - self.lam(T, w), 0.0)
        add = self.nt if self.donor else 0.0
        off = 0.0 if self.donor else self.nt
        return Q / EPS * (self.ixd(0.0, w) + add * x1 * x1 / 2.0 - off * (w * w - x1 * x1) / 2.0)

    def width(self, T, v):
        rhs = self.vbi(T) - v - KB * T
        lo, hi = 1e-7, 5e-3
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            if self.potential(T, mid) > rhs:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    def window(self, T, w_inf, v_pulse):
        wp = self.width(T, v_pulse)
        x1 = max(wp - self.lam(T, wp), 0.0)
        x2 = max(w_inf - self.lam(T, w_inf), x1 + 1e-7)
        return x1, x2

    def field(self, T, x, x2, w_inf):
        add = self.nt if self.donor else 0.0
        off = 0.0 if self.donor else self.nt
        xa = np.minimum(x, x2)
        inner = self.idd(xa, x2) + add * (x2 - xa)
        tail = self.idd(x2, w_inf) - off * (w_inf - x2)
        return Q / EPS * (inner + tail)

    def pf_gamma(self, T, x, x2, w_inf):
        if not self.donor:
            return np.zeros_like(x)
        F = self.field(T, x, x2, w_inf)
        return np.sqrt(Q * np.maximum(F, 0.0) / (math.pi * EPS)) / (KB * T)

    def with_trap(self, level, nt, donor):
        return Layer(self.kx, self.kn, self.phi_b, level, nt, donor)


# ---------------------------------------------------------------------------------------------------------------------
def _two_param(S, y, sig):
    """Weighted least squares of y = a + b s for every row s of S; returns chi2, a, b arrays."""
    w = 1.0 / sig ** 2
    sw, sy = w.sum(), (w * y).sum()
    ss = S @ w
    sss = (S * S) @ w
    ssy = S @ (w * y)
    det = sw * sss - ss * ss
    det = np.where(np.abs(det) < 1e-300, 1e-300, det)
    b = (sw * ssy - ss * sy) / det
    a = (sy - b * ss) / sw
    r = (y[None, :] - a[:, None] - b[:, None] * S) / sig[None, :]
    return (r * r).sum(axis=1), a, b


def fit_shape(t, y, sig, shape_of_rate, lo=-1.0, hi=5.5):
    """y = a + b * shape(e0, t); best e0 on a log grid, refined; returns (e0, a, b, chi2)."""
    grid = np.linspace(lo, hi, 66)
    chi, _, _ = _two_param(np.array([shape_of_rate(10.0 ** g, t) for g in grid]), y, sig)
    k = int(np.argmin(chi))
    a, b = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]

    def one(le):
        return _two_param(shape_of_rate(10.0 ** le, t)[None, :], y, sig)
    res = minimize_scalar(lambda le: float(one(le)[0][0]), bounds=(a, b), method="bounded",
                          options={"xatol": 1e-6})
    chi2, ca, cb = one(res.x)
    return 10.0 ** res.x, float(ca[0]), float(cb[0]), float(chi2[0])


def shape_maker(layer, T, w_inf, v_pulse, n=240):
    x1, x2 = layer.window(T, w_inf, v_pulse)
    if not layer.donor:
        return (lambda e0, t: np.exp(-e0 * t)), x1, x2
    x = np.linspace(x1, x2, n)
    wgt = x / np.trapezoid(x, x)
    gam = layer.pf_gamma(T, x, x2, w_inf)
    mult = np.exp(gam)

    def shape(e0, t):
        return np.trapezoid(wgt[None, :] * np.exp(-np.outer(t, e0 * mult)), x, axis=1)
    return shape, x1, x2


def usable(bench, rec):
    """A record is usable while the meter still resolves the junction (C_m > 0, series-resistance term small)."""
    cm = rec["C_pF"] * 1e-12
    gm = rec["G_uS"] * 1e-6
    return bool(np.all(cm > 0) and np.all(gm / (bench.w * np.abs(cm)) < 2.0))


def read_transient(bench, layer, rec, v_rev):
    keep = rec["t_s"] >= T_SKIP
    t = rec["t_s"][keep]
    c = bench.junction_c(rec["C_pF"][keep], rec["G_uS"][keep])
    y = 1.0 / c ** 2
    sig = 2.0 * 1.5e-15 / c ** 3              # 1/C^2 noise from 1.5 fF on C
    T = float(np.mean(rec["T_K"]))
    w_guess = EPS * bench.area / c[-1]
    shape, x1, x2 = shape_maker(layer, T, w_guess, bench.vp)
    e0, a, b, chi2 = fit_shape(t, y, sig, shape)
    w_inf = EPS * bench.area * math.sqrt(a)
    shape, x1, x2 = shape_maker(layer, T, w_inf, bench.vp)
    e0, a, b, chi2 = fit_shape(t, y, sig, shape)
    w_inf = EPS * bench.area * math.sqrt(a)
    return {"T": T, "e0": e0, "a": a, "b": b, "chi2": chi2, "n": len(t), "w_inf": w_inf, "x1": x1, "x2": x2,
            "snr": abs(b) / (np.median(sig) + 1e-30), "v": v_rev}


def density_from(layer, tr, fa=1.0):
    """Exact charge balance at t = 0 (traps of the window filled to fa)."""
    w_inf = tr["w_inf"]
    w0 = EPS * tr["scale"] * math.sqrt(tr["a"] + tr["b"])
    x1, x2 = tr["x1"], tr["x2"]
    phi = layer.ixd(x2, w0) - layer.ixd(x2, w_inf)
    s0 = fa * (x2 * x2 - x1 * x1) / 2.0
    if layer.donor:
        return phi / s0
    return phi / ((w0 * w0 - w_inf * w_inf) / 2.0 + s0)


def cv_widths(bench):
    c = bench.junction_c(bench.cv["C_pF"], bench.cv["G_uS"])
    return bench.cv["bias_v"], EPS * bench.area / c


def fit_cv(bench, layer):
    """Donor density on knots every 0.1 um across the depths the C-V curve reaches, and the barrier height."""
    v, w = cv_widths(bench)
    kx = np.arange(w.min(), w.max() + 1e-5, 1e-5)

    def make(p):
        return Layer(kx, np.exp(p[:-1]) * 1e15, p[-1], layer.level, layer.nt, layer.donor)

    def resid(p):
        L = make(p)
        pred = np.array([L.vbi(300.0) - L.potential(300.0, wi) - KB * 300.0 for wi in w])
        smooth = 0.05 * np.diff(p[:-1], 2)
        return np.concatenate([(pred - v) / 0.002, smooth / 0.01])
    p0 = np.concatenate([np.log(layer.nd(kx) / 1e15), [layer.phi_b]])
    r = least_squares(resid, p0, x_scale=np.concatenate([np.full(len(kx), 0.02), [0.01]]))
    return make(r.x)


def first_cv(bench):
    v, w = cv_widths(bench)
    slope, icpt = np.polyfit(v, w * w, 1)
    n = -2.0 * EPS / (Q * slope)
    vbi = icpt / (-slope) + KB * 300.0
    return Layer([0.0], [n], vbi + KB * 300.0 * math.log(nc(300.0) / n), 0.45, 0.0, False)


def arrhenius(points):
    T = np.array([p["T"] for p in points])
    e = np.array([p["e0"] for p in points])
    wts = np.array([min(p["snr"], 400.0) for p in points])
    X = 1.0 / (KB * T)
    Y = np.log(e / T ** 2)
    A = np.column_stack([np.ones_like(X), -X]) * wts[:, None]
    coef, *_ = np.linalg.lstsq(A, Y * wts, rcond=None)
    gamma = vth(300.0) * nc(300.0) / 300.0 ** 2
    return float(coef[1]), float(math.exp(coef[0]) / gamma)


def capture_rates(bench, layer):
    """Capture rate at each capture temperature from the amplitude against pulse width; None if not resolved."""
    by_T = {}
    for rec in bench.capture:
        T = round(float(np.mean(rec["T_K"])) / 2.5) * 2.5          # set points sit on a 2.5 K grid
        by_T.setdefault(T, []).append(rec)
    out = {}
    for T, recs in sorted(by_T.items()):
        widths, amps = [], []
        for rec in recs:
            tr = read_transient(bench, layer, rec, bench.vr)
            widths.append(float(rec["pulse_width_s"][0]))
            amps.append(tr["b"])
        widths, amps = np.array(widths), np.array(amps)
        sat = np.median(amps[-3:])
        frac = amps / sat
        if frac[0] > 0.85:                     # already filled at the shortest pulse: not resolved
            out[T] = None
            continue
        def res(p):
            return (sat * p[1] * (1.0 - np.exp(-10.0 ** p[0] * widths)) - amps) / (0.01 * abs(sat))
        r = least_squares(res, [math.log10(1.0 / widths[np.argmin(abs(frac - 0.63))]), 1.0])
        s2 = 2.0 * r.cost / max(len(widths) - 2, 1)
        cov = np.linalg.pinv(r.jac.T @ r.jac) * s2
        out[T] = (10.0 ** r.x[0], math.sqrt(max(cov[0, 0], 0.0)) * math.log(10.0))
    return out


def reduce(folder, talk=False):
    bench = Bench(folder)
    to_sample = bench.sample_temperatures()
    for rec in bench.scan + bench.bias + bench.capture:
        rec["T_K"] = to_sample(rec["T_K"])
    bench.scan = [r for r in bench.scan if usable(bench, r)]
    bench.bias = [r for r in bench.bias if usable(bench, r)]
    bench.capture = [r for r in bench.capture if usable(bench, r)]
    layer = first_cv(bench)
    level, nt, donor, eb = 0.45, 0.0, False, 0.0
    sig_a = 1e-15
    # records whose single-rate reading lies inside the time base carry the rate; the rest are left out
    plain = layer.with_trap(level, nt, False)
    inside = []
    for rec in bench.scan:
        tr = read_transient(bench, plain, rec, bench.vr)
        if 1.0 <= tr["e0"] <= 6000.0 and tr["snr"] > 30:
            inside.append(rec)
    for it in range(3):
        layer = layer.with_trap(level, nt, donor)
        layer = fit_cv(bench, layer)
        # which kind of centre: the bias series under both descriptions
        if it == 1:
            spread = {}
            for kind in (False, True):
                L = layer.with_trap(level, nt, kind)
                trs = [read_transient(bench, L, rec, float(rec["reverse_bias_v"][0])) for rec in bench.bias]
                by_T = {}
                for tr in trs:
                    by_T.setdefault(round(tr["T"] / 2.5), []).append(math.log(tr["e0"]))
                spread[kind] = float(np.mean([np.ptp(v) for v in by_T.values()]))
            donor = spread[True] < spread[False]
            if talk:
                print("bias-series spread of ln e0: neutral-after-emission", spread[False], "charged", spread[True])
            layer = layer.with_trap(level, nt, donor)
        pts = []
        for rec in inside:
            tr = read_transient(bench, layer, rec, bench.vr)
            if 1.5 <= tr["e0"] <= 4000.0 and tr["snr"] > 30:
                pts.append(tr)
        ea, sig_a = arrhenius(pts)
        caps = capture_rates(bench, layer)
        sig_c = {T: (c[0] / (vth(T) * layer.n_free(3e-4)), c[1]) for T, c in caps.items() if c is not None}
        eb = 0.0
        if len(sig_c) >= 2:
            (t1, (s1, u1)), (t2, (s2, u2)) = sorted(sig_c.items())[:2]
            span = 1.0 / t1 - 1.0 / t2
            est = KB * math.log(s2 / s1) / span
            unc = KB * math.hypot(u1, u2) / span
            # a temperature dependence of the capture cross-section is taken only when the series resolves it
            eb = est if est > 3.0 * unc else 0.0
        level = ea - eb
        dens = []
        for tr in pts:
            if 5.0 <= tr["e0"] <= 1500.0:
                tr["scale"] = bench.area
                L = layer.with_trap(level, nt, donor)
                x1, x2 = L.window(tr["T"], tr["w_inf"], bench.vp)
                tr["x1"], tr["x2"] = x1, x2
                dens.append(density_from(L, tr))
        nt = float(np.median(dens))
        if talk:
            print(f"pass {it}: points {len(pts)} Ea {ea:.4f} sigma_a {sig_a:.3e} Eb {eb:.4f} level {level:.4f} "
                  f"NT {nt:.4e} donor {donor} N(1um) {float(layer.nd(1e-4)):.4e} N(2.5um) {float(layer.nd(2.5e-4)):.4e}")
    sigma300 = sig_a * math.exp(-eb / (KB * 300.0))
    return {"level_ev": round(level, 5), "density_cm3": float(f"{nt:.5e}"),
            "sigma_300k_cm2": float(f"{sigma300:.5e}"),
            "centre": "donor-like" if donor else "acceptor-like", "capture_barrier_ev": round(eb, 4),
            "doping_cm3_at_2um": float(f"{float(layer.nd(2e-4)):.4e}")}


def main():
    out = reduce(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
