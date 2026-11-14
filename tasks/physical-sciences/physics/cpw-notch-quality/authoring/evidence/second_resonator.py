"""Second, independent reduction of one resonator (same call as the reference, different numerics).

    python second_resonator.py <resonator directory> <result.json>

- Parameters: f_r, ln Q_i(-96 dBm), ln of the ratio Q_i(-90 dBm) / Q_i(-96 dBm), ln |Q_c|, phi, sqrt(A at -90 dBm),
  delay, slope, and an amplitude and phase per trace (Q_l at each power follows from Q_i and Q_c, and
  A(-96 dBm) = 10^-0.6 (Q_l(-96) / Q_l(-90))^3 A(-90 dBm)), fitted by Levenberg-Marquardt least squares.
- Kerr detuning: Newton iterations from the linear detuning where the cubic's discriminant shows a single real root;
  where it shows three, or Newton has not reached a root, the roots are the eigenvalues of the companion matrix and
  the sweep direction picks the outer real one on its side (up: the smallest, down: the largest).
- Average over the wandering: a shift of f_r only moves the instant response along frequency, so the response is
  tabulated once per trace at half the trace's frequency step (as cell averages: the cell holding the jump of a
  bistable trace mixes the two branches by the parts of the cell on either side of the fold), convolved directly
  (no FFT) with a Gaussian kernel truncated at +-5 sigma with trapezoid end weights, and read by linear
  interpolation at f minus each state's offset of the dominating fluctuator (weights 1 - p and p); a trial step of
  a fit that makes the wandering wider than the sweep has it capped at the sweep's span.
- Wandering from the time stream: an algebraic (Kasa) circle fit to the gain-corrected stream samples, whose
  off-resonance point is the circle point nearest to the fitted environment at the probe frequency; a sample's angle
  2 beta from that point gives the Kerr-shifted detuning y = tan(beta) / 2, the Kerr law gives
  y_0 = y - A / (1 + 4 y^2), and the radial scatter about the circle gives the amplifier noise to subtract.  The
  wandering sequence counts as a dominating fluctuator plus a Gaussian part when its lag-1 autocorrelation exceeds
  0.15 (independent Gaussian draws have none); then an equal-width two-Gaussian mixture fitted by expectation-
  maximisation (the samples' time order unused) gives sigma, the pull d and the fraction up p; otherwise sigma is the
  rms of the wandering less the amplifier noise.
- Order (every fit but the last on every third trace point): the -96 dBm trace alone from three Kerr strengths, then
  three rounds of wandering estimate (the first Gaussian only) and -96 dBm fit; the -90 dBm trace joins (equal Q_i,
  its amplitude and phase from the -96 dBm trace) from two Kerr strengths (and both directions where it can be
  bistable) and the lowest misfit is kept; wandering
  estimate and joint fit alternate (up to three rounds); the last joint fit, on all points, frees d and p with sigma
  held.
- CW receiver gain, model-free: the own points of the trace at the stream's drive around the probe frequency (a
  quadratic in f through the nine nearest, real and imaginary parts) give the averaged response there, and
  |stream mean| over it is the gain; the stream is divided by it before anything else (the header's logged gain is
  not used).
- Drive of the stream (the header's is not trusted): the whole reduction is run with the stream at -96 dBm and at
  -90 dBm; at -90 dBm the stream's detuning is read with Q_l(-90 dBm) and A(-90 dBm), and its pulls are carried to -96 dBm linewidths by (Q_l(-96) / Q_l(-90)) 10^0.15 (until the -90 dBm
  trace has joined, Q_l(-90) = Q_l(-96), and such a stream gives only a Gaussian wandering); the run whose last joint
  fit leaves the smaller cost is kept (a run that breaks down counts as infinite cost).
- Sweep direction of the -90 dBm trace: where the -96 dBm trace's Kerr strength carried to -90 dBm exceeds 0.5 (a
  trace that can be bistable), the -90 dBm trace joins with each direction's branch rule and the lowest misfit
  decides (the header's direction is not trusted).
"""
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402



def read(path):
    with open(path, newline="") as fh:
        return np.array([[float(v) for v in r.values()] for r in csv.DictReader(fh)])


def solve_y(y0, amp, up):
    if amp == 0:
        return y0
    # discriminant of 4 y^3 - 4 y0 y^2 + y - (y0 + A)
    a3, a2, a1, a0 = 4.0, -4 * y0, 1.0, -(y0 + amp)
    disc = (18 * a3 * a2 * a1 * a0 - 4 * a2 ** 3 * a0 + a2 ** 2 * a1 ** 2 - 4 * a3 * a1 ** 3 - 27 * a3 ** 2 * a0 ** 2)
    y = y0 + amp / (1 + 4 * y0 ** 2)
    for _ in range(60):                                                   # Newton on the single-root points
        g = 4 * y ** 3 - 4 * y0 * y ** 2 + y - (y0 + amp)
        dg = 12 * y ** 2 - 8 * y0 * y + 1
        step = g / np.where(np.abs(dg) > 1e-12, dg, 1e-12)
        y = y - step
        if np.max(np.abs(step)) < 1e-12:
            break
    g = 4 * y ** 3 - 4 * y0 * y ** 2 + y - (y0 + amp)
    near = (disc > 0) | ~(np.abs(g) < 1e-9)   # three real roots, or Newton stalled near the other branch's fold
    if near.any():
        m = np.zeros((near.sum(), 3, 3))
        m[:, 0, 0] = y0[near]
        m[:, 0, 1] = -0.25
        m[:, 0, 2] = (y0[near] + amp) / 4
        m[:, 1, 0] = m[:, 2, 1] = 1.0
        ev = np.linalg.eigvals(m)
        real = np.abs(ev.imag) <= 1e-5 * (1 + np.abs(ev.real))   # a double root comes out with a small imaginary part
        lone = real.sum(axis=1) < 2                           # one real root: the eigenvalue nearest the real axis
        pick = np.where(real, ev.real, np.inf if up else -np.inf)
        pick = pick.min(axis=1) if up else pick.max(axis=1)
        alone = ev.real[np.arange(ev.shape[0]), np.argmin(np.abs(ev.imag), axis=1)]
        y = y.copy()
        y[near] = np.where(lone, alone, pick)
    return y


def response(f, fr, ql, qc, phi, amp, up):
    y = solve_y(ql * (f / fr - 1), amp, up)
    return 1 - (ql / qc) * np.exp(1j * phi) / (1 + 2j * y)


def folds(amp):
    """The two fold points (y, y0) of y0 = y - A / (1 + 4 y^2) past bifurcation: the real negative roots of
    16 y^4 + 8 y^2 + 8 A y + 1 = 0, where dy0/dy vanishes."""
    r = np.roots([16.0, 0.0, 8.0, 8.0 * amp, 1.0])
    y = np.sort(r[np.abs(r.imag) < 1e-9].real)
    return [(v, v - amp / (1 + 4 * v * v)) for v in y]


def cell_table(u, fr, ql, qc, phi, amp, up):
    """The instant response on a uniform grid u as cell averages: the cell holding the jump (where the branch the
    sweep follows ends) gets the two branches in proportion to the parts of the cell on either side of it, so the
    table moves continuously with the jump."""
    tab = response(u, fr, ql, qc, phi, amp, up)
    if amp <= 4 * np.sqrt(3) / 9:
        return tab
    fl = folds(amp)
    if len(fl) < 2:
        return tab
    (ya, y0max), (yb, y0min) = fl[0], fl[1]
    y0j, yend = (y0max, ya) if up else (y0min, yb)
    uj = fr * (1 + y0j / ql)
    step = u[1] - u[0]
    i = int(np.floor((uj - u[0]) / step + 0.5))
    if not 0 <= i < u.size:
        return tab
    left = (uj - (u[i] - step / 2)) / step                    # part of the cell below the jump frequency
    y0i = ql * (u[i] / fr - 1)
    roots = np.sort(np.roots([4.0, -4 * y0i, 1.0, -(y0i + amp)]))
    real = roots[np.abs(roots.imag) < 1e-9].real
    before = (real.min() if up else real.max()) if (y0i < y0max if up else y0i > y0min) else yend
    after = real.max() if up else real.min()
    val = lambda y: 1 - (ql / qc) * np.exp(1j * phi) / (1 + 2j * y)  # noqa: E731
    lo, hi = (val(before), val(after)) if up else (val(after), val(before))
    tab[i] = left * lo + (1 - left) * hi
    return tab


class Run:
    def __init__(self, folder):
        folder = Path(folder)
        self.meta = json.loads((folder / "device.json").read_text())
        self.tr = []
        for name, info in self.meta["sweeps"].items():
            d = read(folder / name)
            self.tr.append([d[:, 0], d[:, 1] + 1j * d[:, 2], info["drive"].startswith("-90"), info["direction"] == "up"])
        self.k96 = next(i for i, t in enumerate(self.tr) if not t[2])
        d = read(folder / self.meta["time_stream"]["file"])
        self.z = d[:, 0] + 1j * d[:, 1]
        self.fp = float(self.meta["time_stream"]["probe_frequency_hz"])
        self.fref = float(self.tr[0][0].mean())
        self.w = [1 / ((0.0008 if t[2] else 0.0016) * np.abs(t[1][:50]).mean()) for t in self.tr]
        self.high = False

    def par(self, p):
        fr = self.fref * (1 + p[0] * 1e-6)
        qi96, qc, phi = np.exp(p[1]), np.exp(p[3]), p[4]
        qi90 = qi96 * np.exp(p[2])
        ql96 = 1 / (1 / qi96 + np.cos(phi) / qc)
        ql90 = 1 / (1 / qi90 + np.cos(phi) / qc)
        amp90 = p[5] ** 2
        return fr, ql96, ql90, qc, phi, amp90 * 10 ** -0.6 * (ql96 / ql90) ** 3, amp90, p[6], p[7]

    def model(self, p, k, wand, thin=1):
        fr, ql96, ql90, qc, phi, amp96, amp90, tau, b = self.par(p)
        f, _, high, up = self.tr[k]
        f = f[::thin]
        ql, amp = (ql90, amp90) if high else (ql96, amp96)
        env = p[8 + 2 * k] * (1 + b * (f - fr) / 1e6) * np.exp(1j * (p[9 + 2 * k] - 2 * np.pi * f * tau * 1e-9))
        sig, d, pu = wand
        hz = fr / ql96 * (10 ** -0.15 if high else 1.0)                  # wandering in hertz per -96 dBm linewidth
        span = f.max() - f.min()                       # a trial step that widens the line past the sweep is capped there
        s_hz, d_hz = min(sig * hz, span), min(d * hz, span)
        states = [(-pu * d_hz, 1 - pu), ((1 - pu) * d_hz, pu)] if d > 0 else [(0.0, 1.0)]
        if s_hz < 1.0:
            return env * sum(w * response(f, fr + c, ql, qc, phi, amp, up) for c, w in states)
        df = self.tr[k][0][1] - self.tr[k][0][0]
        step = df / 2
        half = int(np.ceil(5 * s_hz / step))
        reach = (half + 2) * step + d_hz
        table_f = np.arange(f.min() - reach, f.max() + reach + step, step)
        table = cell_table(table_f, fr, ql, qc, phi, amp, up)        # a shift of f_r only moves the response in f
        kern = np.exp(-0.5 * (np.arange(-half, half + 1) * step / s_hz) ** 2)
        kern[[0, -1]] *= 0.5                                          # trapezoid end weights
        kern /= kern.sum()
        smooth = np.convolve(table.real, kern, mode="same") + 1j * np.convolve(table.imag, kern, mode="same")
        return env * sum(w * (np.interp(f - c, table_f, smooth.real) + 1j * np.interp(f - c, table_f, smooth.imag))
                         for c, w in states)

    def fit(self, p, wand, which, nfev=300, free=False, thin=1):
        """Least squares over the traces in `which`; only the parameters those traces depend on move (with the
        -96 dBm trace alone: not the Q_i ratio nor the -90 dBm amplitude and phase); `free` adds d and p."""
        p = np.asarray(p, float)
        idx = [0, 1, 3, 4, 5, 6, 7, 8 + 2 * which[0], 9 + 2 * which[0]] if len(which) == 1 else list(range(12))

        def full(q):
            x = p.copy()
            x[idx] = q[:len(idx)]
            return x

        def res(q):
            wd = (wand[0], abs(q[-2]), min(max(q[-1], 0.05), 0.95)) if free else wand
            x = full(q)
            out = []
            for k in which:
                r = (self.model(x, k, wd, thin) - self.tr[k][1][::thin]) * self.w[k]
                out += [r.real, r.imag]
            return np.concatenate(out)
        x0 = np.r_[p[idx], wand[1], wand[2]] if free else p[idx]
        sol = least_squares(res, x0, method="lm", x_scale="jac", max_nfev=nfev, ftol=1e-10, xtol=1e-10)
        sol.wand = (wand[0], abs(sol.x[-2]), min(max(sol.x[-1], 0.05), 0.95)) if free else wand
        sol.x = full(sol.x)
        return sol

    def wandering(self, p):
        """Wandering of each stream sample in -96 dBm linewidths (up to a constant) and the amplifier noise variance
        along it: an algebraic circle fit to the gain-corrected samples; the circle's off-resonance point is the one
        nearest to the fitted environment at the probe frequency; a sample's angle 2 beta from it gives the detuning
        y = tan(beta) / 2, and the Kerr law y_0 = y - A / (1 + 4 y^2)."""
        fr, ql96, ql90, qc, phi, amp96, amp90, tau, b = self.par(p)
        k = self.k96                                    # the cables are the same for every record
        amp, scale = (amp90, ql96 / ql90 * 10 ** 0.15) if self.high else (amp96, 1.0)
        env = p[8 + 2 * k] * (1 + b * (self.fp - fr) / 1e6) * np.exp(1j * (p[9 + 2 * k] - 2 * np.pi * self.fp * tau * 1e-9))
        x, y = self.z.real, self.z.imag
        cx, cy, c0 = np.linalg.lstsq(np.column_stack([x, y, np.ones_like(x)]), x * x + y * y, rcond=None)[0]
        c = (cx + 1j * cy) / 2
        rad = np.sqrt(c0 + abs(c) ** 2)
        far = c + rad * (env - c) / abs(env - c)
        rot = (self.z - c) / (far - c)
        beta = -np.angle(-rot) / 2                             # z - c = -(far - c) exp(-2 i beta)
        yk = np.tan(beta) / 2
        y0 = yk - amp / (1 + 4 * yk ** 2)
        radial = np.abs(self.z - c) - rad
        dy = 1 / (4 * np.cos(beta) ** 2 * rad)
        return -y0 * scale, np.mean(radial ** 2) * dy ** 2 * scale ** 2

    def fluctuator(self, p, allow=True):
        ek, w2 = self.wandering(p)
        x = ek - ek.mean()
        rho = float(np.mean(x[1:] * x[:-1]) / np.mean(x * x))
        if rho < 0.15 or not allow:
            return (float(np.sqrt(max(np.var(ek) - np.mean(w2), 0.0))), 0.0, 0.0)
        m0, m1 = np.percentile(ek, [20, 80])
        s2, pu = max(np.var(ek) / 4, 1e-6), 0.5
        for _ in range(300):                                   # equal-width two-Gaussian mixture, EM
            v = s2 + w2
            l0 = (1 - pu) * np.exp(-0.5 * (ek - m0) ** 2 / v)
            l1 = pu * np.exp(-0.5 * (ek - m1) ** 2 / v)
            g1 = l1 / (l0 + l1)
            g0 = 1 - g1
            n0, n1 = (g0 / v).sum(), (g1 / v).sum()
            new0, new1 = (g0 * ek / v).sum() / n0, (g1 * ek / v).sum() / n1
            s2n = max(float(np.mean(g0 * (ek - new0) ** 2 + g1 * (ek - new1) ** 2 - w2)), 1e-8)
            pun = float(g1.mean())
            done = abs(new0 - m0) + abs(new1 - m1) + abs(s2n - s2) + abs(pun - pu) < 1e-9
            m0, m1, s2, pu = new0, new1, s2n, pun
            if done:
                break
        if m1 < m0:
            m0, m1, pu = m1, m0, 1 - pu
        return (float(np.sqrt(s2)), float(m1 - m0), float(pu))


def run(folder):
    """The reduction with the stream at each drive; the one whose last joint fit has the smaller cost (a run whose fit
    breaks down, a wandering that turns non-finite, counts as infinite cost)."""
    runs = []
    for high in (False, True):
        try:
            runs.append(reduce_at(folder, high))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            runs.append((np.inf, None))
    return min(runs, key=lambda t: t[0])[1]


def reduce_at(folder, high):
    r = Run(folder)
    r.high = high
    f, s, _, _ = r.tr[1 - r.k96 if high else r.k96]
    near = np.argsort(np.abs(f - r.fp))[:9]
    x = (f[near] - r.fp) / 1e3
    at_fp = np.polyval(np.polyfit(x, s[near].real, 2), 0.0) + 1j * np.polyval(np.polyfit(x, s[near].imag, 2), 0.0)
    r.z = r.z / (np.abs(r.z.mean()) / np.abs(at_fp))
    f, s, _, _ = r.tr[r.k96]
    lvl = np.abs(np.r_[s[:50], s[-50:]]).mean()
    i0 = int(np.argmin(np.abs(s)))
    depth = 1 - abs(s[i0]) / lvl
    half = np.flatnonzero(np.abs(s) / lvl < 1 - depth / 2)
    ql0 = f[i0] / max(f[half[-1]] - f[half[0]], 3 * (f[1] - f[0]))
    qc0 = ql0 / max(depth, 0.05)
    qi0 = 1 / max(1 / ql0 - 1 / qc0, 0.2 / ql0)
    alpha0 = float(np.angle(s[:20].mean()))
    none = (0.0, 0.0, 0.0)
    best = None
    for a0 in (0.2, 0.6, 1.2):                                  # the -96 dBm trace alone, three Kerr strengths
        p0 = np.r_[(f[i0] / r.fref - 1) * 1e6, np.log(qi0), 0.0, np.log(qc0), 0.0, np.sqrt(a0), 0.0, 0.0,
                   lvl, alpha0, lvl, alpha0]
        sol = r.fit(p0, none, [r.k96], nfev=100, thin=3)
        if best is None or sol.cost < best.cost:
            best = sol
    p, wand = best.x, none
    for rnd in range(3):                    # the wandering on the -96 dBm trace alone: Gaussian first, then the rest
        wand = r.fluctuator(p, allow=rnd > 0 and not high)  # a -90 dBm stream gives only its rms until that trace joins
        p = r.fit(p, wand, [r.k96], nfev=100, thin=3).x
    k90 = 1 - r.k96
    joint = []
    logged = r.tr[k90][3]
    # the -90 dBm trace joins from two Kerr strengths, and where the -96 dBm trace's Kerr strength carried to -90 dBm
    # exceeds 0.5 (a trace that can be bistable) with each direction's branch rule (the header is not trusted)
    for up in ((logged, not logged) if p[5] ** 2 * 10 ** 0.6 > 0.5 else (logged,)):
        r.tr[k90][3] = up
        for scale in (1.0, 0.4):
            q = p.copy()
            q[2] = 0.0
            q[8 + 2 * k90], q[9 + 2 * k90] = p[8 + 2 * r.k96], p[9 + 2 * r.k96]
            q[5] = p[5] * np.sqrt(scale * 10 ** 0.6)
            sol = r.fit(q, wand, [0, 1], nfev=120, thin=3)
            joint.append((sol.cost, up, sol.x))
    _, r.tr[k90][3], p = min(joint, key=lambda t: t[0])
    for rnd in range(3):
        new = r.fluctuator(p)
        moved = sum(abs(a - b) for a, b in zip(new, wand))
        p, wand = r.fit(p, new, [0, 1], nfev=120, thin=3).x, new
        if moved < 3e-3:
            break
    sol = r.fit(p, wand, [0, 1], free=wand[1] > 0)
    fr, ql96, _, qc, phi = r.par(sol.x)[:5]
    return sol.cost, {"internal_quality_factor": float(1 / (1 / ql96 - np.cos(phi) / qc)),
                      "resonance_frequency_ghz": float(fr / 1e9)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(run(sys.argv[1])))
