"""Reduces one notch-type resonator measurement to its internal quality factor at -96 dBm and its resonance frequency.

    python3 resonator_reduce.py <resonator directory> <result.json>

1. Instant response: the notch model with the Kerr-shifted detuning y, a real root of
   4 y^3 - 4 y0 y^2 + y - (y0 + A) = 0 (y0 = Q_l (f / f_r' - 1)), solved in closed form for every point; where three
   roots are real the branch follows the sweep direction (upward sweep: the smallest root, the low-amplitude branch
   continued from below; downward sweep: the largest).  Q_l is fitted at each drive power (the two-level systems'
   loss saturates); |Q_c| and phi are shared; A(-96 dBm) = 10^-0.6 (Q_l(-96) / Q_l(-90))^3 A(-90 dBm).
2. Trace model: the instant response averaged over the distribution of the wandering e of f_r' (Gauss-Hermite
   quadrature on each state's Gaussian), every pull in hertz 10^-0.15 times smaller at -90 dBm, times the environment
   a (1 + b (f - f_r)) e^{i(alpha - 2 pi f tau)}; both traces are fitted together (shared resonator, Kerr A, delay and
   slope; Q_l, amplitude and phase per trace).
3. The wandering from the time stream: every sample is projected onto the curve the instant response traces at the
   probe frequency as f_r' wanders (environment and resonator from the current fit), with the amplifier noise carried
   along the curve known from the spread across it.  The stream's drive is read from its own circle, not from the
   header: the circle's radius over the distance of its off-resonance point from the origin is the notch diameter
   Q_l / |Q_c| whatever the CW gain, and the drive whose Q_l gives the nearer diameter is taken; a stream at -90 dBm
   is projected with Q_l(-90 dBm) and A(-90 dBm), and its pulls, 10^-0.15 of those at -96 dBm, are scaled back.  The projected sequence is fitted by a two-state hidden Markov
   model (one dominating two-state fluctuator: pull d, fraction p of the time up, Markov switching) on top of a
   Gaussian part of rms sigma; where the two-state model does not beat one Gaussian clearly, e is Gaussian.  The
   stream is first divided by the CW receiver gain measured on the data: the stream's mean is the averaged response
   at the probe frequency, so |mean| over the trace model's averaged value there is the gain (the header's logged
   gain is not used).  Steps 2 and 3 alternate; the last trace fit refines d and p with sigma held.
4. Sweep direction of the -90 dBm trace: the joint fit is started with the branch rule of each direction and the
   direction with the lower cost is kept (the header's direction is only the first guess); where A(-90 dBm) exceeds
   0.5 the other direction is tried again once the wandering is known, and if it fits better, stream and trace fit
   alternate again with it.
5. Q_i = 1 / (1/Q_l(-96) - cos(phi)/|Q_c|); f_r is the mean small-signal resonance frequency.
"""
import csv
import json
import os
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402
from scipy.signal import fftconvolve  # noqa: E402

KERR_RATIO = 10 ** -0.6          # A(-96 dBm) / A(-90 dBm) at equal Q_l
PULL_RATIO = 10 ** -0.15         # every pull in hertz, -90 dBm over -96 dBm
CLOCK = [time.monotonic()]
REFINE = [True]


def load(path):
    with open(path, newline="") as fh:
        rows = [[float(v) for v in r.values()] for r in csv.DictReader(fh)]
    return np.array(rows)


def detuning(y0, amp, up):
    """Real root of 4 y^3 - 4 y0 y^2 + y - (y0 + A) = 0 on the branch the sweep direction selects."""
    y0 = np.asarray(y0, float)
    if amp == 0:
        return y0.copy()
    p = 0.25 - y0 * y0 / 3                                # y = t + y0/3 gives t^3 + p t + q = 0
    q = -2 * y0 ** 3 / 27 + y0 / 12 - (y0 + amp) / 4
    disc = (q / 2) ** 2 + (p / 3) ** 3
    sq = np.sqrt(np.clip(disc, 0, None))
    t = np.cbrt(-q / 2 + sq) + np.cbrt(-q / 2 - sq)
    three = disc <= 0
    if three.any():                                       # bistable stretch: three real roots, pick the branch
        pp = -p[three] / 3
        base = np.arccos(np.clip(3 * q[three] / (2 * p[three]) * np.sqrt(1 / pp), -1, 1)) / 3
        roots = 2 * np.sqrt(pp) * np.cos(base[None, :] - 2 * np.pi * np.arange(3)[:, None] / 3)
        t[three] = roots.min(axis=0) if up else roots.max(axis=0)
    return t + y0 / 3


def instant(f, fr, ql, qc, phi, amp, up):
    y = detuning(ql * (f / fr - 1), amp, up)
    return 1 - (ql / qc) * np.exp(1j * phi) / (1 + 2j * y)


def averaged(f, fr, ql, qc, phi, amp, up, noise, high, lw96, grid=4096):
    """Instant response averaged over the distribution of e: noise = (sigma, d, p), sigma and d in linewidths of the
    -96 dBm trace at -96 dBm (every pull in hertz 10^-0.15 times smaller at -90 dBm).  A shift of f_r' by e only moves
    the instant response along frequency (to 1 part in f_r / e), so the average is the instant response convolved
    with the Gaussian on a fine grid (exact also across the jumps past bifurcation), read at f minus each state's
    offset of the dominating fluctuator."""
    sig, d, p = noise
    scale = (PULL_RATIO if high else 1.0) * lw96
    sig_hz, d_hz = sig * scale, d * scale
    states = [(-p * d_hz, 1 - p), ((1 - p) * d_hz, p)] if d_hz > 0 else [(0.0, 1.0)]
    if sig_hz < 1.0:
        return sum(w * instant(f, fr + c, ql, qc, phi, amp, up) for c, w in states)
    lo, hi = f.min() - 7 * sig_hz - d_hz, f.max() + 7 * sig_hz + d_hz
    u = np.linspace(lo, hi, grid)
    step = u[1] - u[0]
    half = int(6 * sig_hz / step) + 1
    k = np.arange(-half, half + 1) * step
    g = np.exp(-0.5 * (k / sig_hz) ** 2)
    conv = fftconvolve(instant(u, fr, ql, qc, phi, amp, up), g / g.sum(), mode="same")
    return sum(w * (np.interp(f - c, u, conv.real) + 1j * np.interp(f - c, u, conv.imag)) for c, w in states)


class Device:
    def __init__(self, folder, directions=True, flip_high=False):
        folder = Path(folder)
        self.head = json.loads((folder / "device.json").read_text())
        self.flip = flip_high
        self.traces = []
        for name, info in self.head["sweeps"].items():
            d = load(folder / name)
            high = info["drive"].startswith("-90")
            up = info["direction"] == "up" or not directions
            self.traces.append((d[:, 0], d[:, 1] + 1j * d[:, 2], high, (not up) if (high and flip_high) else up))
        self.traces.sort(key=lambda t: t[2])                 # -96 dBm first
        # residual weights from the data's own off-resonance level (fixed, so the fit cannot buy them)
        self.weight = [1 / ((0.0008 if t[2] else 0.0016) * np.abs(np.r_[t[1][:100], t[1][-100:]]).mean())
                       for t in self.traces]
        ts = self.head["time_stream"]
        d = load(folder / ts["file"])
        self.stream = d[:, 0] + 1j * d[:, 1]
        self.logged_gain = 10 ** (float(ts.get("receiver_gain_db", 0.0)) / 20)
        self.used = self.stream / self.logged_gain
        self.fp = float(ts["probe_frequency_hz"])
        self.f0 = float(np.median(np.concatenate([t[0] for t in self.traces])))


# parameters: ppm, ln Q_l(-96), ln Q_l(-90), ln |Q_c|, phi, sqrt A(-90), tau ns, b per MHz, then a, alpha per trace
def unpack(p, f0, shared_q=False):
    fr = f0 * (1 + p[0] * 1e-6)
    ql96, ql90 = np.exp(p[1]), np.exp(p[1] if shared_q else p[2])
    amp90 = p[5] * p[5]
    return fr, ql96, ql90, np.exp(p[3]), p[4], amp90 * KERR_RATIO * (ql96 / ql90) ** 3, amp90, p[6], p[7]


def env_at(p, f0, k, f):
    fr, tau, b = f0 * (1 + p[0] * 1e-6), p[6], p[7]
    a, alpha = p[8 + 2 * k], p[9 + 2 * k]
    return a * (1 + b * (f - fr) / 1e6) * np.exp(1j * (alpha - 2 * np.pi * f * tau * 1e-9))


def trace_model(p, dev, k, noise, kerr=True, shared_q=False, sub=1):
    fr, ql96, ql90, qc, phi, amp96, amp90 = unpack(p, dev.f0, shared_q)[:7]
    f, _, high, up = dev.traces[k]
    f = f[::sub]
    ql, amp = (ql90, amp90) if high else (ql96, amp96)
    return env_at(p, dev.f0, k, f) * averaged(f, fr, ql, qc, phi, amp if kerr else 0.0, up, noise, high, fr / ql96)


def start(dev):
    """Linear notch reading of the -96 dBm trace for starting values."""
    f, s, _, _ = dev.traces[0]
    edge = np.r_[np.arange(f.size // 6), np.arange(f.size - f.size // 6, f.size)]
    a0 = np.abs(s[edge]).mean()
    i0 = int(np.argmin(np.abs(s)))
    depth = 1 - np.abs(s[i0]) / a0
    below = np.flatnonzero(np.abs(s) / a0 < 1 - depth / 2)
    fwhm = max(f[below[-1]] - f[below[0]], 3 * (f[1] - f[0])) if below.size else (f[-1] - f[0]) / 8
    ql0 = f[i0] / fwhm
    qc0 = ql0 / max(depth, 0.05)
    alpha0 = np.angle(s[edge[0]])
    ppm0 = (f[i0] / dev.f0 - 1) * 1e6
    return np.r_[ppm0, np.log(ql0), np.log(ql0), np.log(qc0), 0.0, 0.0, 0.0, 0.0, a0, alpha0, a0, alpha0]


def residuals(q, dev, noise, kerr=True, shared_q=False, sub=1):
    out = []
    for k, (f, s, high, _) in enumerate(dev.traces):
        r = (trace_model(q, dev, k, noise, kerr, shared_q, sub) - s[::sub]) * dev.weight[k]
        out += [r.real, r.imag]
    return np.concatenate(out)


def fit(dev, p0, noise, kerr=True, shared_q=False, nfev=300, sub=1, free_tel=False):
    """Joint fit of both traces; with free_tel the fluctuator's pull d and occupation p are refined too (sigma held)."""
    if not free_tel:
        sol = least_squares(residuals, np.asarray(p0, float), method="lm", x_scale="jac", max_nfev=nfev, ftol=1e-11,
                            xtol=1e-11, args=(dev, noise, kerr, shared_q, sub))
        sol.noise = noise
        return sol
    tel = lambda x: (noise[0], abs(x[12]), min(max(x[13], 0.05), 0.95))  # noqa: E731
    sol = least_squares(lambda x: residuals(x[:12], dev, tel(x), kerr, shared_q, sub), np.r_[p0, noise[1], noise[2]],
                        method="lm", x_scale="jac", max_nfev=nfev, ftol=1e-11, xtol=1e-11)
    sol.noise = tel(sol.x)
    sol.x = sol.x[:12]
    return sol


def circle(z):
    """Centre and radius of the circle through the stream samples (algebraic fit refined geometrically)."""
    x, y = z.real, z.imag
    m = np.c_[x, y, np.ones_like(x)]
    sol = np.linalg.lstsq(m, x * x + y * y, rcond=None)[0]
    c = complex(sol[0] / 2, sol[1] / 2)
    rad = np.sqrt(sol[2] + abs(c) ** 2)
    for _ in range(20):                                   # geometric refinement of the algebraic circle
        r = np.abs(z - c)
        u = (z - c) / r
        jac = np.c_[-u.real, -u.imag, -np.ones_like(r)]
        step = np.linalg.lstsq(jac, -(r - rad), rcond=None)[0]
        c += complex(step[0], step[1])
        rad += step[2]
        if np.abs(step).max() < 1e-12 * max(rad, 1e-30):
            break
    return c, rad


def stream_diameter(dev, p):
    """The notch diameter Q_l / |Q_c| of the stream's own circle, free of the CW gain: radius over the distance of the
    circle's off-resonance point from the origin, that point lying in the direction arg(environment) + phi from the
    centre (current trace fit)."""
    phi = unpack(p, dev.f0)[4]
    c, rad = circle(dev.stream)
    off = c + rad * np.exp(1j * (np.angle(env_at(p, dev.f0, 0, dev.fp)) + phi))
    return 2 * rad / abs(off)


def project(dev, p, noise, gain=True, high=False):
    """Wandering of every stream sample in linewidths of the -96 dBm trace at -96 dBm (up to a constant), and the
    variance the amplifier noise gives it, from the geometry of the resonance circle the samples lie on: the circle is
    fitted to the samples themselves (so the CW gain does not enter), the off-resonance point of the circle lies in
    the direction arg(environment) + phi from its centre (current trace fit), the angle 2 beta from that point gives
    the Kerr-shifted detuning y = tan(beta) / 2, and the Kerr law gives y_0 = y - A / (1 + 4 y^2) with the A of the
    stream's drive.  A stream at -90 dBm counts its detuning in the -90 dBm linewidth and its pulls are 10^-0.15 of
    those at -96 dBm, so its wandering is multiplied by (Q_l(-96) / Q_l(-90)) / 10^-0.15."""
    fr, ql96, ql90, qc, phi, amp96, amp90 = unpack(p, dev.f0)[:7]
    amp, scale = (amp90, ql96 / ql90 / PULL_RATIO) if high else (amp96, 1.0)
    z = dev.stream
    c, rad = circle(z)
    far = np.angle(env_at(p, dev.f0, 0, dev.fp)) + phi
    beta = np.angle(np.exp(1j * (far + np.pi)) / ((z - c) / np.abs(z - c))) / 2
    yk = np.tan(beta) / 2
    y0 = yk - amp / (1 + 4 * yk * yk)
    w2 = np.mean((np.abs(z - c) - rad) ** 2)               # one quadrature of the amplifier noise
    dy = 0.5 / np.cos(beta) ** 2 * (1 - amp * (-8 * yk) / (1 + 4 * yk * yk) ** 2)
    return -y0 * scale, w2 / (2 * rad) ** 2 * dy ** 2 * scale ** 2


def project_curve(dev, p, noise, gain=True, high=False):
    """Stream samples projected onto the instant-response curve at the probe frequency of the current fit at the
    stream's drive (the stream divided by the measured or the logged CW gain): wandering in linewidths of the -96 dBm
    trace at -96 dBm (a -90 dBm stream's pulls 10^-0.15 of those) and along-curve noise variance."""
    fr, ql96, ql90, qc, phi, amp96, amp90 = unpack(p, dev.f0)[:7]
    lw = fr / ql96
    ql, amp, pull = (ql90, amp90, PULL_RATIO) if high else (ql96, amp96, 1.0)
    env = env_at(p, dev.f0, 0, dev.fp)
    if gain:
        pred = env * averaged(np.array([dev.fp]), fr, ql, qc, phi, amp, True, noise, high, lw)[0]
        dev.used = dev.stream / (np.abs(dev.stream.mean()) / np.abs(pred))
    grid = np.linspace(-4.0, 4.0, 3201)
    curve = env * instant(np.full(grid.size, dev.fp), fr + grid * lw * pull, ql, qc, phi, amp, True)
    best = np.empty(dev.used.size, int)
    for k in range(0, dev.used.size, 500):
        best[k:k + 500] = np.argmin(np.abs(dev.used[k:k + 500, None] - curve[None, :]), axis=1)
    j = np.clip(best, 1, grid.size - 2)
    d2 = lambda k: np.abs(dev.used - curve[k]) ** 2  # noqa: E731
    y1, y2, y3 = d2(j - 1), d2(j), d2(j + 1)
    den = y1 - 2 * y2 + y3
    shift = 0.5 * (y1 - y3) / np.where(np.abs(den) > 0, den, 1)
    step = grid[1] - grid[0]
    eps = grid[j] + np.clip(shift, -1, 1) * step
    slope = np.abs(curve[j + 1] - curve[j - 1]) / (2 * step)
    across = np.mean(np.min(np.c_[y1, y2, y3], axis=1))
    return eps, across / slope ** 2


def hmm(eps, w2, iters=40):
    """Two-state hidden Markov model with Gaussian emissions (means m0 < m1, variance s2 + w2_k) by Baum-Welch:
    returns (s2, m0, m1, p up, a00 + a11 - 1 the state's lag-1 correlation, log-likelihood)."""
    n = eps.size
    xs, ws = eps.tolist(), w2.tolist()
    m0, m1 = (float(v) for v in np.percentile(eps, [25, 75]))
    s2 = max(float(np.var(eps)) / 4, 1e-12)
    a00, a01, a10, a11 = 0.8, 0.2, 0.2, 0.8
    ll_old = -np.inf
    for _ in range(iters):
        b0, b1 = [], []
        for x, w in zip(xs, ws):
            v = s2 + w
            k = 1 / np.sqrt(2 * np.pi * v)
            b0.append(max(k * np.exp(-0.5 * (x - m0) ** 2 / v), 1e-300))
            b1.append(max(k * np.exp(-0.5 * (x - m1) ** 2 / v), 1e-300))
        f0, f1, cs = [0.0] * n, [0.0] * n, [0.0] * n
        u0, u1 = 0.5 * b0[0], 0.5 * b1[0]
        c = u0 + u1
        f0[0], f1[0], cs[0] = u0 / c, u1 / c, c
        for t in range(1, n):
            u0 = (f0[t - 1] * a00 + f1[t - 1] * a10) * b0[t]
            u1 = (f0[t - 1] * a01 + f1[t - 1] * a11) * b1[t]
            c = u0 + u1
            f0[t], f1[t], cs[t] = u0 / c, u1 / c, c
        g0, g1 = [0.0] * n, [0.0] * n
        r0, r1 = 1.0, 1.0
        x00 = x01 = x10 = x11 = 0.0
        g0[-1], g1[-1] = f0[-1], f1[-1]
        for t in range(n - 2, -1, -1):
            e0, e1 = b0[t + 1] * r0 / cs[t + 1], b1[t + 1] * r1 / cs[t + 1]
            x00 += f0[t] * a00 * e0
            x01 += f0[t] * a01 * e1
            x10 += f1[t] * a10 * e0
            x11 += f1[t] * a11 * e1
            r0, r1 = a00 * e0 + a01 * e1, a10 * e0 + a11 * e1
            q0, q1 = f0[t] * r0, f1[t] * r1
            g0[t], g1[t] = q0 / (q0 + q1), q1 / (q0 + q1)
        a00, a01 = x00 / (x00 + x01), x01 / (x00 + x01)
        a10, a11 = x10 / (x10 + x11), x11 / (x10 + x11)
        g = np.c_[g0, g1]
        wts = g / (s2 + w2)[:, None]
        m0, m1 = (float(v) for v in (wts * eps[:, None]).sum(axis=0) / wts.sum(axis=0))
        s2 = max(float(np.mean(g[:, 0] * (eps - m0) ** 2 + g[:, 1] * (eps - m1) ** 2 - w2)), 1e-12)
        ll = float(np.sum(np.log(cs)))
        if abs(ll - ll_old) < 1e-4:
            break
        ll_old = ll
    if m1 < m0:
        m0, m1, a00, a11, a01, a10 = m1, m0, a11, a00, a10, a01
    return s2, m0, m1, a01 / (a01 + a10), a00 + a11 - 1, ll


def stream_noise(dev, p, noise, gain=True, telegraph=True, geometry=True):
    """(sigma, d, p) of the wandering from the stream, in linewidths of the -96 dBm trace.  A dominating fluctuator
    is taken only where the two-state model fits clearly better than one Gaussian and its states persist from sample
    to sample (lag-1 correlation of the state above 0.3); a single Gaussian distorted by a wrong curve does neither."""
    eps, w2 = (project if geometry else project_curve)(dev, p, noise, gain, dev.high)
    gauss_s2 = max(float(np.var(eps) - np.mean(w2)), 0.0)
    if not telegraph:
        return (np.sqrt(gauss_s2), 0.0, 0.0)
    v = gauss_s2 + w2
    ll1 = float(np.sum(-0.5 * (eps - eps.mean()) ** 2 / v - 0.5 * np.log(2 * np.pi * v)))
    s2, m0, m1, pu, rho, ll2 = hmm(eps, w2)
    if ll2 - ll1 < 25 or rho < 0.3 or m1 - m0 < 0.05:
        return (np.sqrt(gauss_s2), 0.0, 0.0)
    return (np.sqrt(s2), m1 - m0, pu)


def reduce_device(folder, stream_drive="fit", **kw):
    """Options other than the defaults are the shortcuts tabulated by authoring/evidence/bench.py.  The stream's drive
    (stream_drive="fit"): the reduction is made with the header's drive first; where the stream circle's notch
    diameter differs from Q_l / |Q_c| at that drive by more than 3 %, it is made again with the other drive (started
    from the first result) and the drive whose diameter agrees better is kept.  "logged": the header's drive only."""
    high = json.loads((Path(folder) / "device.json").read_text())["time_stream"]["drive"].startswith("-90")
    out, mismatch = solve(folder, high, **kw)
    if stream_drive == "logged" or not kw.get("jitter", True) or not kw.get("stream", True) or mismatch <= 0.03:
        return out
    other, other_mismatch = solve(folder, not high, warm=solve.state, **kw)
    return other if other_mismatch < mismatch else out


def solve(folder, high, jitter=True, kerr=True, stream=True, directions="fit", cw_gain="measured", telegraph=True,
          shared_q=False, geometry=True, warm=None):
    """Reduction with the stream taken at -90 dBm (high) or at -96 dBm: (result, relative difference between the
    stream circle's diameter and Q_l / |Q_c| at that drive).  warm: (parameters, flip, wandering) of an earlier run
    to start from instead of the multistart."""
    CLOCK[0] = time.monotonic()
    options = (False, True) if directions == "fit" and warm is None else (False,)
    best, dev = None, None
    for flip in options:
        if warm is not None:
            dev = Device(folder, directions is not False, warm[1])
            dev.high = high
            break
        trial = Device(folder, directions is not False, flip)
        trial.high = high
        p = start(trial)
        if shared_q:
            p[2] = p[1]
        for amp0 in (0.05, 0.4, 1.0, 2.0):
            q = p.copy()
            q[5] = np.sqrt(amp0)
            sol = fit(trial, q, (0.0, 0.0, 0.0), kerr, shared_q, nfev=200, sub=4)
            if best is None or sol.cost < best.cost:
                best, dev = sol, trial
    p = best.x if warm is None else warm[0]
    noise = (0.0, 0.0, 0.0) if warm is None else warm[2]
    if jitter and stream:
        def settle(dev, p, noise, rounds, first):
            for rnd in range(rounds):
                new = stream_noise(dev, p, noise, cw_gain == "measured", telegraph and (rnd > 0 or not first), geometry)
                sol = fit(dev, p, new, kerr, shared_q, sub=2)
                moved = abs(new[0] - noise[0]) + abs(new[1] - noise[1]) + abs(new[2] - noise[2])
                p, noise = sol.x, new
                if (moved < 2e-3 and rnd >= 2) or time.monotonic() - CLOCK[0] > 30:
                    break
            return p, noise
        p, noise = settle(dev, p, noise, 6, warm is None)
        if directions == "fit" and unpack(p, dev.f0)[6] > 0.5:
            other = Device(folder, True, not dev.flip)
            other.high = high
            alt = fit(other, p, noise, kerr, shared_q, sub=2)
            mine = 0.5 * np.sum(residuals(p, dev, noise, kerr, shared_q, 2) ** 2)
            if alt.cost < mine:                     # the other direction: the wandering is estimated again with it
                dev = other
                p, noise = settle(dev, alt.x, noise, 4, False)
        sol = fit(dev, p, noise, kerr, shared_q, free_tel=noise[1] > 0 and REFINE[0])
        p = sol.x
    elif jitter:
        x0 = np.r_[p, 0.1]

        p = least_squares(lambda x: residuals(x[:12], dev, (abs(x[12]), 0.0, 0.0), kerr, shared_q), x0,
                          method="lm", x_scale="jac", max_nfev=300).x[:12]
    fr, ql96, ql90, qc, phi = unpack(p, dev.f0, shared_q)[:5]
    qi = 1 / (1 / ql96 - np.cos(phi) / qc)
    mismatch = abs(stream_diameter(dev, p) / ((ql90 if high else ql96) / qc) - 1)
    solve.state = (p, dev.flip, noise)
    return {"internal_quality_factor": float(qi), "resonance_frequency_ghz": float(fr / 1e9)}, mismatch


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_device(sys.argv[1])))
