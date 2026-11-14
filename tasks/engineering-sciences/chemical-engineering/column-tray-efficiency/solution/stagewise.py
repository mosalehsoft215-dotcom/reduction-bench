"""Steady-state binary tray column: NRTL liquid, ideal vapour, Murphree vapour efficiency on every tray, full energy
balance with heat loss per tray, subcooled reflux, total or partial condenser, kettle reboiler as an equilibrium stage,
optional liquid side draw. Molar units inside (mol/s, J/mol, Pa, K).

Stage numbering: trays 1..N from the top, reboiler N+1. Condenser handled apart.
"""
import numpy as np
from scipy.optimize import root

R_GAS = 8.314462618
T0 = 298.15


class Mixture:
    """Binary mixture, component 1 the lighter one. Data dictionary as written in the plant sheet."""

    def __init__(self, d):
        self.d = d
        self.M = np.array([d["light"]["molar_mass_g_mol"], d["heavy"]["molar_mass_g_mol"]]) / 1000.0
        self.ant = np.array([[c["antoine"]["A"], c["antoine"]["B"], c["antoine"]["C"]] for c in (d["light"], d["heavy"])])
        self.cpl = np.array([d["light"]["cp_liquid_J_molK"], d["heavy"]["cp_liquid_J_molK"]])
        self.cpv = np.array([d["light"]["cp_vapour_J_molK"], d["heavy"]["cp_vapour_J_molK"]])
        self.lam = np.array([d["light"]["dH_vap_nbp_J_mol"], d["heavy"]["dH_vap_nbp_J_mol"]])
        n = d["nrtl"]
        self.a12, self.b12, self.a21, self.b21, self.alpha = n["a12"], n["b12_K"], n["a21"], n["b21_K"], n["alpha"]
        self.Tb = np.array([self._tsat(i, 101325.0) for i in range(2)])

    def psat(self, T):
        """Pa; Antoine in log10(mmHg) with t in C. T may be an array; returns shape (..., 2)."""
        t = np.asarray(T)[..., None] - 273.15
        A, B, C = self.ant[:, 0], self.ant[:, 1], self.ant[:, 2]
        return 10.0 ** (A - B / (C + t)) * 133.322368

    def _tsat(self, i, P):
        A, B, C = self.ant[i]
        return B / (A - np.log10(P / 133.322368)) - C + 273.15

    def gamma(self, x1, T):
        x1 = np.asarray(x1)
        x2 = 1.0 - x1
        t12 = self.a12 + self.b12 / T
        t21 = self.a21 + self.b21 / T
        G12 = np.exp(-self.alpha * t12)
        G21 = np.exp(-self.alpha * t21)
        l1 = x2 ** 2 * (t21 * (G21 / (x1 + x2 * G21)) ** 2 + t12 * G12 / (x2 + x1 * G12) ** 2)
        l2 = x1 ** 2 * (t12 * (G12 / (x2 + x1 * G12)) ** 2 + t21 * G21 / (x1 + x2 * G21) ** 2)
        return np.stack([np.exp(l1), np.exp(l2)], axis=-1)

    def bubble(self, x1, P):
        """Bubble temperature and equilibrium vapour mole fraction of the light component."""
        x1 = np.asarray(x1, dtype=float)
        P = np.broadcast_to(np.asarray(P, dtype=float), x1.shape)
        T = x1 * self._tsat(0, P) + (1 - x1) * self._tsat(1, P)
        X = np.stack([x1, 1 - x1], axis=-1)
        for _ in range(30):
            k = X * self.gamma(x1, T) * self.psat(T)
            s = k.sum(-1)
            f = np.log(s / P)
            k2 = X * self.gamma(x1, T + 0.01) * self.psat(T + 0.01)
            df = (np.log(k2.sum(-1) / P) - f) / 0.01
            step = f / df
            T = T - step
            if np.all(np.abs(step) < 1e-9):
                break
        k = X * self.gamma(x1, T) * self.psat(T)
        return T, k[..., 0] / k.sum(-1)

    def dew_liquid(self, y1, P):
        """Liquid in equilibrium with vapour y1 at P (partial condenser): returns x1, T."""
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            _, ys = self.bubble(mid, P)
            if ys < y1:
                lo = mid
            else:
                hi = mid
        x = 0.5 * (lo + hi)
        T, _ = self.bubble(x, P)
        return x, T

    def hL(self, x1, T):
        return x1 * self.cpl[0] * (T - T0) + (1 - x1) * self.cpl[1] * (T - T0)

    def hV(self, y1, T):
        hv = self.cpl * (self.Tb - T0) + self.lam + self.cpv * (np.asarray(T)[..., None] - self.Tb)
        return y1 * hv[..., 0] + (1 - y1) * hv[..., 1]

    def molar_mass(self, x1):
        return x1 * self.M[0] + (1 - x1) * self.M[1]


def solve_column(mix, col, op, E_rect, E_strip, guess=None):
    """col: N, feed_tray, side_tray (or None), condenser 'total'|'partial', dP_tray_Pa, dP_condenser_Pa, dP_reboiler_Pa,
    heat_loss_W_K (array per tray, UA), op: F, zF (mol frac), hF (J/mol), L0 (reflux mol/s), TR (reflux T, total only),
    D (mol/s), S (side draw mol/s), P_top (Pa, top tray vapour space), T_amb.
    E_rect applies to trays above the feed tray (1..f-1), E_strip to f..N."""
    N, f = col["N"], col["feed_tray"]
    side = col.get("side_tray")
    P = op["P_top"] + col["dP_tray_Pa"] * np.arange(N)  # tray pressures, tray 1 = top
    P_reb = P[-1] + col["dP_reboiler_Pa"]
    Pst = np.append(P, P_reb)
    E = np.where(np.arange(1, N + 1) < f, E_rect, E_strip)
    UA = np.asarray(col["heat_loss_W_K"], dtype=float)
    F, zF, hF, L0, D = op["F"], op["zF"], op["hF"], op["L0"], op["D"]
    S = op.get("S", 0.0) if side else 0.0
    V1 = L0 + D
    Fj = np.zeros(N + 1)
    Fj[f - 1] = F
    Sj = np.zeros(N + 1)
    if side:
        Sj[side - 1] = S
    # water entering the kettle from outside the column (a leaking steam tube), mol/s; carries no light component
    Wj = np.zeros(N + 1)
    Wj[N] = op.get("W_reb", 0.0)

    partial = col["condenser"] == "partial"
    Pd = op["P_top"] - col["dP_condenser_Pa"]

    def profile(z):
        x = np.clip(z[: N + 1], 1e-10, 1 - 1e-10)
        V = np.concatenate([[V1], z[N + 1:2 * N + 1], [0.0]])  # V[j] leaving stage j+1 (index 0 = tray 1), V[N+1] = 0
        T, ys = mix.bubble(x, Pst)
        y = np.empty(N + 1)
        y[N] = ys[N]
        for j in range(N - 1, -1, -1):
            y[j] = y[j + 1] + E[j] * (ys[j] - y[j + 1])
        return x, V, T, y

    def reflux_state(z, y1):
        if not partial:
            return y1, op["TR"], y1, None
        # partial condenser: the reflux composition x0 is the last unknown; the drum balance
        # V1 y1 = L0 x0 + D yD, with yD in equilibrium with x0 at drum pressure, is its equation
        x0 = float(np.clip(z[-1], 1e-12, 1 - 1e-12))
        TR, yD = mix.bubble(x0, Pd)
        return x0, float(TR), float(yD), (L0 * x0 + D * float(yD) - V1 * y1) / F

    def resid(z):
        x, V, T, y = profile(z)
        x0, TR, _, drum = reflux_state(z, y[0])
        Lin = np.empty(N + 1)
        L = np.empty(N + 1)
        Lin[0] = L0
        for j in range(N + 1):
            L[j] = Lin[j] + V[j + 1] - V[j] + Fj[j] - Sj[j] + Wj[j]
            if j < N:
                Lin[j + 1] = L[j]
        xin = np.concatenate([[x0], x[:-1]])
        Tin = np.concatenate([[TR], T[:-1]])
        yin = np.concatenate([y[1:], [0.0]])
        Tvin = np.concatenate([T[1:], [T[-1]]])
        comp = Lin * xin + V[1:] * yin + Fj * zF - L * x - V[:-1] * y - Sj * x
        hin = mix.hL(xin, Tin)
        en = (Lin * hin + V[1:] * mix.hV(yin, Tvin) + Fj * hF - L * mix.hL(x, T) - V[:-1] * mix.hV(y, T) - Sj * mix.hL(x, T))
        en_tr = en[:N] - UA * (T[:N] - op["T_amb"])
        out = np.concatenate([comp / F, en_tr / (F * 4e4)])
        return out if drum is None else np.append(out, drum)

    size = 2 * N + 1 + (1 if partial else 0)
    starts = [] if guess is None or len(guess) != size else [np.asarray(guess, dtype=float)]
    for top, bot in ((0.95, 0.02), (0.99, 0.005), (0.8, 0.05), (0.999, 0.001)):
        z0 = np.concatenate([np.linspace(top, bot, N + 1), np.full(N, V1)])
        starts.append(np.append(z0, top) if partial else z0)
    for start in starts:
        sol = root(resid, start, method="hybr", options={"xtol": 1e-12, "maxfev": 20000})
        if not np.all(np.isfinite(sol.x)) or np.max(np.abs(resid(sol.x))) > 1e-9:
            sol = root(resid, sol.x if np.all(np.isfinite(sol.x)) else start, method="lm",
                       options={"xtol": 1e-13, "ftol": 1e-13, "maxiter": 20000})
        if np.all(np.isfinite(sol.x)) and np.max(np.abs(resid(sol.x))) < 1e-9:
            break
    x, V, T, y = profile(sol.x)
    x0, TR, yD, _ = reflux_state(sol.x, y[0])
    Lin = np.empty(N + 1)
    L = np.empty(N + 1)
    Lin[0] = L0
    for j in range(N + 1):
        L[j] = Lin[j] + V[j + 1] - V[j] + Fj[j] - Sj[j] + Wj[j]
        if j < N:
            Lin[j + 1] = L[j]
    # reboiler duty from its energy balance
    B = L[N]
    Q = B * mix.hL(x[N], T[N]) + V[N] * mix.hV(y[N], T[N]) - L[N - 1] * mix.hL(x[N - 1], T[N - 1])
    res = np.max(np.abs(resid(sol.x)))
    return dict(x=x, y=y, T=T, V=V, L=L, xD=yD if col["condenser"] == "partial" else y[0], x0=x0, TR=TR, xB=x[N],
                B=B, Q_reb=Q, P=Pst, ok=bool(res < 1e-9), resid=res, z=sol.x)
