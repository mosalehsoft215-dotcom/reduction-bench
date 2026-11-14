"""Physical model of an Au / n-Si Schottky diode holding one electron trap, as measured on a DLTS bench.

Depletion approximation, abrupt metal side, depth x from the metal (cm). Units: cm, V, F, s, eV.

Charge in the epi layer (per q):  N+(x, t) = N_d(x) + D N_T - N_T f(x, t)
  D = 1 for a donor-like centre (positive when empty), 0 for an acceptor-like one (negative when filled).
Bias condition:                     (q / eps) * int_0^W x N+(x, t) dx = V_bi(T) - V - kT/q
Lambda edge:                        lambda = sqrt(2 eps (E_F - E_T) / (q N_d(edge))),  E_c - E_F = kT ln(N_c / n)
Occupancy at reverse bias:          f = 0 for x < W_p - lambda_p (never filled), f = 1 beyond W_inf - lambda_inf,
                                    f = f_a exp(-e_n(x) t) in between, f_a the repetitive-pulsing steady state.
Rates:                              e_n = sigma_n(T) v_th N_c exp(-(dE - dE_PF(x)) / kT),  c_n = sigma_n v_th n,
                                    sigma_n(T) = sigma_inf exp(-E_b / kT),
                                    dE_PF = sqrt(q F(x) / (pi eps)) for a donor-like centre, 0 otherwise.
Meter (1 MHz):                      R_s(T) in series with (C parallel G_leak), reported as parallel C_m and G_m.
"""
from dataclasses import dataclass, field

import numpy as np

Q = 1.602176634e-19
K_EV = 8.617333262e-5
EPS = 11.7 * 8.8541878128e-14


def n_c(T):
    return 2.86e19 * (T / 300.0) ** 1.5


def v_th(T):
    return 2.30e7 * np.sqrt(T / 300.0)


@dataclass
class Trap:
    dE: float              # E_c - E_T at zero field, eV
    NT: float              # cm^-3
    sigma_inf: float       # cm^2
    Eb: float = 0.0        # capture barrier, eV
    donor: bool = False    # True: positive when empty (Poole-Frenkel applies)

    def sigma(self, T):
        return self.sigma_inf * np.exp(-self.Eb / (K_EV * T))


@dataclass
class Diode:
    area: float                 # cm^2
    N0: float                   # donor density at the metal, cm^-3
    grad: float                 # relative increase of N_d per cm of depth
    phi_b: float                # barrier height, eV
    trap: Trap
    rs300: float = 15.0         # ohm
    rs_ea: float = 0.0          # eV, freeze-out activation of R_s
    g_leak: float = 2e-9        # S
    freq: float = 1.0e6         # Hz
    extra: dict = field(default_factory=dict)

    # ---- doping ---------------------------------------------------------------------------------------------------
    # N_d(x) = N0 (1 + grad x) + step * N0 * s((x - x_step) / w_step), s a smooth 0 -> 1 ramp (two-layer epi)
    def __post_init__(self):
        x = np.linspace(0.0, 6e-3, 60001)
        n = self.nd(x)
        self._x = x
        self._i0 = np.concatenate([[0.0], np.cumsum(0.5 * (n[1:] + n[:-1]) * np.diff(x))])
        xn = x * n
        self._i1 = np.concatenate([[0.0], np.cumsum(0.5 * (xn[1:] + xn[:-1]) * np.diff(x))])

    def nd(self, x):
        base = self.N0 * (1.0 + self.grad * np.asarray(x, dtype=float))
        step = self.extra.get("step", 0.0)
        if step:
            z = (np.asarray(x, dtype=float) - self.extra["x_step"]) / self.extra["w_step"]
            base = base + step * self.N0 * 0.5 * (1.0 + np.tanh(z))
        return base

    def _I0(self, x):
        return np.interp(x, self._x, self._i0)

    def _I1(self, x):
        return np.interp(x, self._x, self._i1)

    def bulk_n(self, x):
        """Free-electron density in neutral material at depth x (acceptor-like traps hold N_T electrons)."""
        return self.nd(x) - (0.0 if self.trap.donor else self.trap.NT)

    def ef_depth(self, T, x):
        """E_c - E_F in neutral material at depth x, eV."""
        return K_EV * T * np.log(n_c(T) / self.bulk_n(x))

    def vbi(self, T, x):
        return self.phi_b - self.ef_depth(T, x)

    def lam(self, T, x):
        de = self.trap.dE - self.ef_depth(T, x)
        return np.sqrt(2.0 * EPS * max(de, 1e-4) / (Q * self.nd(x)))

    # ---- piecewise integrals ----------------------------------------------------------------------------------------
    def _int_x_nb(self, a, b):
        """int_a^b x * (N_d(x) - (1-D) N_T) dx, the always-filled region beyond the window."""
        off = (0.0 if self.trap.donor else self.trap.NT)
        return (self._I1(b) - self._I1(a)) - off * (b ** 2 - a ** 2) / 2.0

    def _int_x_nempty(self, a, b):
        """int_a^b x * (N_d(x) + D N_T) dx, traps empty."""
        add = (self.trap.NT if self.trap.donor else 0.0)
        return (self._I1(b) - self._I1(a)) + add * (b ** 2 - a ** 2) / 2.0

    # ---- steady widths --------------------------------------------------------------------------------------------
    def _steady_width(self, T, V, filled_window_start=None):
        """Width with traps empty below the lambda point and filled beyond it (self-consistent lambda).

        filled_window_start: if given, traps are also filled from this depth on (used for the pulse state,
        where the window filled by earlier pulses is irrelevant: the pulse state is computed with the traps
        beyond W_p - lambda_p filled and the rest empty)."""
        rhs = self.vbi(T, 3e-4) - V - K_EV * T
        lo, hi = 1e-7, 5e-3
        for _ in range(80):
            w = 0.5 * (lo + hi)
            x1 = max(w - self.lam(T, w), 0.0)
            pot = Q / EPS * (self._int_x_nempty(0.0, x1) + self._int_x_nb(x1, w))
            if pot > rhs:
                hi = w
            else:
                lo = w
        return 0.5 * (lo + hi)

    def widths(self, T, Vr, Vp):
        w_inf = self._steady_width(T, Vr)
        w_p = self._steady_width(T, Vp)
        x1 = max(w_p - self.lam(T, w_p), 0.0)
        x2 = max(w_inf - self.lam(T, w_inf), x1)
        return w_inf, w_p, x1, x2

    def field_profile(self, T, Vr, x):
        """Steady reverse-bias field (V/cm) at depths x inside the window."""
        w_inf, _, _, x2 = self.widths(T, Vr, 0.0)
        add = (self.trap.NT if self.trap.donor else 0.0)
        off = (0.0 if self.trap.donor else self.trap.NT)
        xa = np.minimum(x, x2)
        inner = (self._I0(x2) - self._I0(xa)) + add * (x2 - xa)
        tail = (self._I0(w_inf) - self._I0(x2)) - off * (w_inf - x2)
        return Q / EPS * (inner + tail)

    # ---- rates -------------------------------------------------------------------------------------------------------
    def emission(self, T, F=None):
        tr = self.trap
        e0 = tr.sigma(T) * v_th(T) * n_c(T) * np.exp(-tr.dE / (K_EV * T))
        if F is None or not tr.donor:
            return e0 if F is None else e0 * np.ones_like(F)
        dpf = np.sqrt(Q * np.maximum(F, 0.0) / (np.pi * EPS))
        return e0 * np.exp(dpf / (K_EV * T))

    def capture(self, T, x=3e-4):
        return self.trap.sigma(T) * v_th(T) * self.bulk_n(x)

    # ---- transient ---------------------------------------------------------------------------------------------------
    def transient(self, T, Vr, Vp, tp, t, t_rep, nx=700):
        """Junction capacitance (F) at times t after the end of a filling pulse."""
        w_inf, w_p, x1, x2 = self.widths(T, Vr, Vp)
        tr = self.trap
        D = 1.0 if tr.donor else 0.0
        if x2 <= x1:
            return np.full_like(t, EPS * self.area / w_inf)
        x = np.linspace(x1, x2, nx)
        F = self.field_profile(T, Vr, x)
        e = self.emission(T, F)
        p = np.exp(-self.capture(T) * tp)
        qd = np.exp(-e * t_rep)
        fa = (1.0 - p) / (1.0 - p * qd)
        f = fa[None, :] * np.exp(-np.outer(t, e))                      # (nt, nx)
        nplus = self.nd(x)[None, :] + D * tr.NT - tr.NT * f
        integ = np.trapezoid(x[None, :] * nplus, x, axis=1)
        base = Q / EPS * (self._int_x_nempty(0.0, x1) + integ)
        rhs = self.vbi(T, 3e-4) - Vr - K_EV * T
        # beyond x2: always-filled profile; solve (Q/EPS) int_{x2}^{W} x N_b dx = rhs - base for W >= x2
        need = (rhs - base) * EPS / Q
        lo = np.full_like(t, x2)
        hi = np.full_like(t, 5e-3)
        for _ in range(70):
            w = 0.5 * (lo + hi)
            val = self._int_x_nb(x2, w)
            big = val > need
            hi = np.where(big, w, hi)
            lo = np.where(big, lo, w)
        w = 0.5 * (lo + hi)
        return EPS * self.area / w

    def steady_cv(self, T, volts):
        return np.array([EPS * self.area / self._steady_width(T, v) for v in volts])

    # ---- sample temperature and forward check -------------------------------------------------------------------------
    def sample_t(self, t_stage):
        """A sample not clamped to the copper sits away from the stage by the radiative load of the 300 K window."""
        d0 = self.extra.get("dT0", 0.0)
        return t_stage + d0 * (1.0 - (t_stage / 300.0) ** 4)

    def forward_v(self, T, current):
        """Thermionic emission over the barrier, ideality n, plus the drop across R_s."""
        n = self.extra.get("ideality", 1.03)
        i_sat = self.area * 112.0 * T ** 2 * np.exp(-self.phi_b / (K_EV * T))
        return n * K_EV * T * np.log(current / i_sat + 1.0) + current * self.rs(T)

    # ---- meter -------------------------------------------------------------------------------------------------------
    def rs(self, T):
        return self.rs300 * np.exp(self.rs_ea / K_EV * (1.0 / T - 1.0 / 300.0))

    def meter(self, C, T):
        w = 2.0 * np.pi * self.freq
        yd = self.g_leak + 1j * w * C
        y = 1.0 / (self.rs(T) + 1.0 / yd)
        return y.imag / w, y.real
