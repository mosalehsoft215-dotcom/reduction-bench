"""Linear d-q model of a round-rotor generator in a sudden short-circuit test (authoring side).

Per unit throughout; inductances are reactances at rated frequency, the rotor is referred to the stator on the
reciprocal (Lad-base) system, and currents use the internal sign convention (the recorded line current is the
negative of the internal phase current).

    state  i = [i_d, i_q, i_fd, i_1d, i_2d, i_1q, i_2q]
    (1/wb) L di/dt = -(R + w S L) i + V,        V = [0, 0, v_f, 0, 0, 0, 0]

w is the per-unit speed, wb the rated angular frequency. The series impedance between the machine terminals and
the short (bus duct, shorting assembly, anything else on the day) is part of the stator rows of L and R; the field
loop carries the winding plus everything the exciter puts in series with it. Phase quantities follow
x_a = x_d cos(th) - x_q sin(th) with th = th0 + w wb t; phases b and c lag by 120 and 240 degrees.

Sequence: open circuit until the second pole closes (t2); a line-to-line short between the first two poles until
the third closes (t3), integrated with RK4 on a fine step; a three-phase short afterwards, solved exactly through
the eigenvectors of the constant-coefficient system.

Before t2 the rotor is either settled (field current constant, dampers dead) or, with pre = (t_step, deficit), still
answering a step of field voltage made at t_step: the field current at t2 is then short of its final value by the
fraction deficit, and the field voltage is the one that holds the final value. E0 is the open-circuit voltage at t2.
"""
import numpy as np

SHIFT = np.array([0.0, -2.0 * np.pi / 3.0, 2.0 * np.pi / 3.0])
ROT = np.zeros((7, 7))
ROT[0, 1], ROT[1, 0] = -1.0, 1.0


def circuit(p, x_series=0.0, r_series=0.0):
    """Inductance and resistance matrices. p holds Ll Lad Laq Lf1d Lfd R fd L1d R1d L2d R2d L1q R1q L2q R2q Ra;
    Lfd and Rfd here are the whole field loop."""
    ls = p["Ll"] + x_series
    ad, aq, c = p["Lad"], p["Laq"], p["Lf1d"]
    m = ad + c
    L = np.array([
        [ls + ad, 0.0, ad, ad, ad, 0.0, 0.0],
        [0.0, ls + aq, 0.0, 0.0, 0.0, aq, aq],
        [ad, 0.0, m + p["Lfd"], m, m, 0.0, 0.0],
        [ad, 0.0, m, m + p["L1d"], m, 0.0, 0.0],
        [ad, 0.0, m, m, m + p["L2d"], 0.0, 0.0],
        [0.0, aq, 0.0, 0.0, 0.0, aq + p["L1q"], aq],
        [0.0, aq, 0.0, 0.0, 0.0, aq, aq + p["L2q"]]])
    R = np.diag([p["Ra"] + r_series, p["Ra"] + r_series, p["Rfd"], p["R1d"], p["R2d"], p["R1q"], p["R2q"]])
    return L, R


def _park(th):
    ang = th + SHIFT
    P = (2.0 / 3.0) * np.vstack([np.cos(ang), -np.sin(ang)])
    dP = (2.0 / 3.0) * np.vstack([-np.sin(ang), -np.cos(ang)])
    Pinv = np.column_stack([np.cos(ang), -np.sin(ang)])
    dPinv = np.column_stack([-np.sin(ang), -np.cos(ang)])
    return P, dP, Pinv, dPinv


class ShortCircuit:
    """One test: machine p, series impedance (x_s, r_s) to the short, open-circuit voltage E0 (phase peak, pu),
    rotor angle th0 at t = 0, speed w (pu), rated angular frequency wb, pole times t2 <= t3, the two phases
    (indices) that close first."""

    def __init__(self, p, x_s, r_s, E0, th0, w, wb, t2, t3, pair, dt_ll=2e-6, pre=None):
        self.p, self.w, self.wb, self.th0 = p, w, wb, th0
        self.t2, self.t3 = t2, t3
        self.L, self.R = circuit(p, x_s, r_s)
        self.Lm, self.Rm = circuit(p, 0.0, 0.0)
        self._open_circuit(E0, pre)
        self.V = np.zeros(7)
        self.V[2] = p["Rfd"] * self.ifd_final
        C = np.zeros(3)
        C[pair[0]], C[pair[1]] = 1.0, -1.0
        self.C = C
        self._line_to_line(dt_ll)
        self._three_phase()

    def _open_circuit(self, E0, pre):
        """Rotor currents before t2. Settled: ifd constant. pre = (t_step, deficit): exact response of the
        open-circuit d-axis rotor (field, damper 1, damper 2) to a field-voltage step at t_step."""
        lad = self.p["Lad"]
        self.pre = pre
        if pre is None:
            self.ifd_final = E0 / (self.w * lad)
            self.x_pre = np.array([self.ifd_final, 0.0, 0.0])
            return
        t_step, deficit = pre
        Lr, Rr = self.L[2:5, 2:5], self.R[2:5, 2:5]
        lam, vec = np.linalg.eig(-self.wb * np.linalg.solve(Lr, Rr))
        self.pre_lam, self.pre_vec = lam, vec
        unit = np.linalg.solve(vec, np.array([-1.0, 0.0, 0.0], complex))       # from rest at r = 0, final 1
        g = (vec @ (unit * np.exp(lam * (self.t2 - t_step)))).real
        shortfall = -g[0]                                   # field current still missing at t2 when r = 0
        r = 1.0 - deficit / shortfall                       # starting level that leaves exactly that deficit
        self.pre_coef = (1.0 - r) * unit                    # deviation from the final state, at t_step
        dev = (vec @ (self.pre_coef * np.exp(lam * (self.t2 - t_step)))).real
        x2 = np.array([1.0, 0.0, 0.0]) + dev
        self.ifd_final = E0 / (self.w * lad * x2.sum())
        self.pre_start = r
        self.x_pre = self.ifd_final * x2

    def _pre_state(self, t):
        """Rotor currents and derivatives (3, n) for t < t2."""
        n = t.size
        if self.pre is None:
            return np.tile(self.x_pre[:, None], (1, n)), np.zeros((3, n))
        t_step = self.pre[0]
        tau = np.maximum(t - t_step, 0.0)
        E = self.pre_coef[:, None] * np.exp(self.pre_lam[:, None] * tau[None, :])
        x = self.ifd_final * (np.array([1.0, 0.0, 0.0])[:, None] + (self.pre_vec @ E).real)
        d = self.ifd_final * (self.pre_vec @ (self.pre_lam[:, None] * E)).real
        before = t < t_step
        x[:, before] = np.array([self.pre_start * self.ifd_final, 0.0, 0.0])[:, None]
        d[:, before] = 0.0
        return x, d

    def _pre_field_antiderivative(self, t):
        """A function whose differences are integrals of i_fd over t < t2."""
        if self.pre is None:
            return self.ifd_final * t
        t_step = self.pre[0]
        tau = t - t_step
        after = np.maximum(tau, 0.0)
        grow = (np.exp(self.pre_lam[:, None] * after[None, :]) - 1.0) / self.pre_lam[:, None]
        G = self.ifd_final * (after + (self.pre_vec[0] @ (self.pre_coef[:, None] * grow)).real)
        return np.where(tau < 0.0, self.pre_start * self.ifd_final * tau, G)

    def theta(self, t):
        return self.th0 + self.w * self.wb * t

    def _G(self, t):
        P, dP, Pinv, dPinv = _park(self.theta(t))
        G = np.zeros((7, 6)); G[:2, 0] = P @ self.C; G[2:, 1:] = np.eye(5)
        dG = np.zeros((7, 6)); dG[:2, 0] = dP @ self.C
        H = np.zeros((6, 7)); H[0, :2] = self.C @ Pinv; H[1:, 2:] = np.eye(5)
        return G, dG, H

    def _rhs(self, t, x):
        G, dG, H = self._G(t)
        lhs = H @ self.L @ G / self.wb
        f = H @ self.V - H @ (self.R @ G @ x + self.w * self.L @ dG @ x + self.w * ROT @ self.L @ G @ x)
        return np.linalg.solve(lhs, f)

    def _line_to_line(self, dt):
        span = self.t3 - self.t2
        x = np.concatenate([[0.0], self.x_pre, [0.0, 0.0]])
        if span <= 0.0:
            self.ll_t = np.array([self.t2]); self.ll_x = x[None, :]; self.ll_dx = np.zeros((1, 6))
            self.x3 = x
            return
        n = max(2, int(np.ceil(span / dt)))
        h = span / n
        ts = self.t2 + h * np.arange(n + 1)
        xs = np.empty((n + 1, 6)); ds = np.empty((n + 1, 6))
        for k in range(n):
            t = ts[k]
            k1 = self._rhs(t, x)
            xs[k], ds[k] = x, k1
            k2 = self._rhs(t + h / 2, x + h / 2 * k1)
            k3 = self._rhs(t + h / 2, x + h / 2 * k2)
            k4 = self._rhs(t + h, x + h * k3)
            x = x + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        xs[n], ds[n] = x, self._rhs(ts[n], x)
        self.ll_t, self.ll_x, self.ll_dx = ts, xs, ds
        self.x3 = x

    def _three_phase(self):
        G, _, _ = self._G(self.t3)
        i3 = G @ self.x3
        Linv = np.linalg.inv(self.L)
        A = -self.wb * Linv @ (self.R + self.w * ROT @ self.L)
        b = self.wb * Linv @ self.V
        self.yss = -np.linalg.solve(A, b)
        lam, vec = np.linalg.eig(A)
        self.lam, self.vec = lam, vec
        self.coef = np.linalg.solve(vec, (i3 - self.yss).astype(complex))

    def state(self, t):
        """dq currents and their time derivatives (7, n) at the times t."""
        t = np.asarray(t, float)
        I = np.zeros((7, t.size)); D = np.zeros((7, t.size))
        m = t < self.t2
        if m.any():
            I[2:5, m], D[2:5, m] = self._pre_state(t[m])
        m = (t >= self.t2) & (t < self.t3)
        if m.any():
            idx = np.nonzero(m)[0]
            tg, xg, dg = self.ll_t, self.ll_x, self.ll_dx
            j = np.clip(np.searchsorted(tg, t[idx]) - 1, 0, len(tg) - 2)
            h = (tg[j + 1] - tg[j])[:, None]
            s = (t[idx] - tg[j])[:, None] / h
            x = ((2 * s ** 3 - 3 * s ** 2 + 1) * xg[j] + (s ** 3 - 2 * s ** 2 + s) * h * dg[j]
                 + (-2 * s ** 3 + 3 * s ** 2) * xg[j + 1] + (s ** 3 - s ** 2) * h * dg[j + 1])
            dx = ((6 * s ** 2 - 6 * s) * xg[j] + (3 * s ** 2 - 4 * s + 1) * h * dg[j]
                  + (-6 * s ** 2 + 6 * s) * xg[j + 1] + (3 * s ** 2 - 2 * s) * h * dg[j + 1]) / h
            for col, k in enumerate(idx):
                G, dG, _ = self._G(t[k])
                I[:, k] = G @ x[col]
                D[:, k] = G @ dx[col] + self.w * self.wb * (dG @ x[col])
        m = t >= self.t3
        if m.any():
            E = self.coef[:, None] * np.exp(self.lam[:, None] * (t[m] - self.t3)[None, :])
            I[:, m] = (self.vec @ E).real + self.yss[:, None]
            D[:, m] = (self.vec @ (self.lam[:, None] * E)).real
        return I, D

    def phase_currents(self, t):
        """Internal phase currents (3, n)."""
        I, _ = self.state(t)
        ang = self.theta(np.asarray(t, float))[None, :] + SHIFT[:, None]
        return I[0] * np.cos(ang) - I[1] * np.sin(ang)

    def terminal_voltages(self, t):
        """Phase voltages at the machine terminals (3, n), upstream of the series impedance."""
        I, D = self.state(t)
        psi = self.Lm @ I
        vdq = self.Rm[:2, :2] @ I[:2] + (self.Lm[:2] @ D) / self.wb + self.w * (ROT @ psi)[:2]
        ang = self.theta(np.asarray(t, float))[None, :] + SHIFT[:, None]
        return vdq[0] * np.cos(ang) - vdq[1] * np.sin(ang)

    def field_integral(self, t, t_start):
        """Integral of i_fd from t_start (< t2) to each t (seconds x pu)."""
        t = np.asarray(t, float)
        G = self._pre_field_antiderivative
        out = G(np.minimum(t, self.t2)) - G(np.full(t.shape, float(t_start)))
        if self.t3 > self.t2:
            tg, fg = self.ll_t, self.ll_x[:, 1]
            cum = np.concatenate([[0.0], np.cumsum(0.5 * (fg[1:] + fg[:-1]) * np.diff(tg))])
            out = out + np.interp(np.clip(t, self.t2, self.t3), tg, cum)
        m = t > self.t3
        if m.any():
            dt = t[m] - self.t3
            grow = (np.exp(self.lam[:, None] * dt[None, :]) - 1.0) / self.lam[:, None]   # every mode decays
            out[m] = out[m] + (self.vec[2] @ (self.coef[:, None] * grow)).real + self.yss[2] * dt
        return out

    def airgap_peak(self, t_end=1.0, step=5e-5):
        """Largest magnitude of the air-gap (magnetising) flux after the fault, pu."""
        t = np.arange(self.t2, self.t2 + t_end, step)
        I, _ = self.state(t)
        ud = self.p["Lad"] * (I[0] + I[2] + I[3] + I[4])
        uq = self.p["Laq"] * (I[1] + I[5] + I[6])
        return float(np.max(np.hypot(ud, uq)))


def operational(p):
    """Operational parameters of the bare machine (d and q), for plausibility checks and design data."""
    L, R = circuit(p)
    out = {}
    for axis, st, rot in (("d", 0, [2, 3, 4]), ("q", 1, [5, 6])):
        Lrr, Rrr = L[np.ix_(rot, rot)], R[np.ix_(rot, rot)]
        open_t = np.sort(1.0 / np.linalg.eigvals(np.linalg.solve(Lrr, Rrr)).real)[::-1]
        Ls = L[st, st]
        shorted = Lrr - np.outer(L[rot, st], L[st, rot]) / Ls
        short_t = np.sort(1.0 / np.linalg.eigvals(np.linalg.solve(shorted, Rrr)).real)[::-1]
        out["X" + axis] = Ls
        x = Ls
        for k in range(len(rot)):
            x = x * short_t[k] / open_t[k]
            out["X" + axis + "'" * (k + 1)] = x
            out["T" + axis + "0" + "'" * (k + 1)] = open_t[k]
            out["T" + axis + "'" * (k + 1)] = short_t[k]
    return out
