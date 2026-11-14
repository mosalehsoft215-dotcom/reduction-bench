"""Physics of a constant-load creep test as simulated for this task.

The specimen between the extensometer ridges is cut into axial slices (parallel length, the two circular fillets
and the bits of shoulder up to the ridges). Every slice creeps by its own true strain e:

    de/dt = A * s^n * exp(-Q / (R T)) * (1 + beta * exp(-e / e_p))        (per second)
    s     = F * exp(e) / a_hot                                             (true stress, constant volume)
    a_hot = a_cold * (1 + alpha_m(T) * (T - 20 C))^2                       (section at the slice temperature)

with T the slice temperature from the furnace profile at that moment. The extensometer reads the change in the
hot, crept, elastically stretched length between the ridges since it was zeroed at temperature before loading.
"""
import numpy as np
from scipy.integrate import solve_ivp

R_GAS = 8.314462618
T_REF_C = 20.0


def table(xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    return lambda v: np.interp(v, xs, ys)


class Slices:
    """Axial slices of one specimen between the ridges (cold dimensions, x = 0 at the middle of the gauge)."""

    def __init__(self, d, Lp, rho, D, ridge_spacing, n_par=160, n_fil=70, n_sh=6):
        rg, rs = d / 2, D / 2
        if rho < rs - rg:
            raise ValueError("fillet radius smaller than the step")
        xf = np.sqrt(rho ** 2 - (rho - (rs - rg)) ** 2)
        sh = (ridge_spacing - Lp) / 2 - xf
        if sh <= 0:
            raise ValueError("ridges inside the fillet")
        e_par = np.linspace(0, Lp / 2, n_par // 2 + 1)
        e_fil = np.linspace(0, xf, n_fil + 1)
        e_sh = np.linspace(0, sh, n_sh + 1)
        u_edges = np.concatenate([e_par, Lp / 2 + e_fil[1:], Lp / 2 + xf + e_sh[1:]])
        u_mid = 0.5 * (u_edges[1:] + u_edges[:-1])
        du = np.diff(u_edges)
        w = u_mid - Lp / 2
        r = np.where(w <= 0, rg, np.where(w < xf, rg + rho - np.sqrt(np.clip(rho ** 2 - w ** 2, 0, None)), rs))
        self.x = np.concatenate([-u_mid[::-1], u_mid])
        self.dx = np.concatenate([du[::-1], du])
        self.area = np.pi * np.concatenate([r[::-1], r]) ** 2
        self.Lp, self.d, self.span = Lp, d, ridge_spacing


class Furnace:
    """Temperature along the specimen: T = T_c(t) + c (1 - (x/H)^2) - c + k x / H, with the curvature c and the tilt k
    growing linearly with the set temperature above the lowest one of the programme (base)."""

    def __init__(self, c0, c1, k0, k1, H=25.0, offset=0.0, base=600.0):
        self.c0, self.c1, self.k0, self.k1, self.H, self.offset, self.base = c0, c1, k0, k1, H, offset, base

    def T_C(self, x, centre_C, set_C):
        c = self.c0 + self.c1 * (set_C - self.base)
        k = self.k0 + self.k1 * (set_C - self.base)
        s = x / self.H
        return centre_C + self.offset + c * (1 - s ** 2) - c + k * s


def setpoint_history(t, set_C, excursions=(), outages=(), room_C=22.0, tau_cool=1.6, tau_heat=0.9):
    """Temperature at the reference point (gauge middle without the profile curvature) at times t (h)."""
    t = np.asarray(t, float)
    T = np.full_like(t, set_C)
    for t1, t2, dT in excursions:
        on = (t >= t1) & (t < t2)
        T[on] += dT * (1 - np.exp(-(t[on] - t1) / 0.4))
        off = t >= t2
        T[off] += dT * (1 - np.exp(-(t2 - t1) / 0.4)) * np.exp(-(t[off] - t2) / 0.4)
    for t1, t2 in outages:
        on = (t >= t1) & (t < t2)
        T[on] = room_C + (set_C - room_C) * np.exp(-(t[on] - t1) / tau_cool)
        low = room_C + (set_C - room_C) * np.exp(-(t2 - t1) / tau_cool)
        after = t >= t2
        T[after] = set_C - (set_C - low) * np.exp(-(t[after] - t2) / tau_heat)
    return T


def run_test(sl, mat, F_N, set_C, furnace, t_end, excursions=(), outages=()):
    """Integrate one test; returns a function giving the extension (mm) and slice temperatures at any time set."""
    alpha, E_GPa = table(*mat["alpha"]), table(*mat["E"])
    lnA, n, Q, beta, ep = mat["lnA"], mat["n"], mat["Q"], mat["beta"], mat["ep"]
    x = sl.x

    def temps(t):
        Tc = setpoint_history(np.atleast_1d(t), set_C, excursions, outages)
        return furnace.T_C(x[None, :], Tc[:, None], set_C)

    def rhs(t, e):
        T = temps(t)[0]
        k = 1 + alpha(T) * (T - T_REF_C)
        a_hot = sl.area * k ** 2
        s = F_N * np.exp(e) / a_hot
        rate = np.exp(lnA + n * np.log(s) - Q / (R_GAS * (T + 273.15)))
        return rate * (1 + beta * np.exp(-e / ep)) * 3600.0

    marks = sorted({0.0, t_end, *[v for a, b, _ in excursions for v in (a, b)], *[v for a, b in outages for v in (a, b)]})
    marks = [m for m in marks if 0 <= m <= t_end]
    pieces, e0 = [], np.zeros_like(x)
    for a, b in zip(marks[:-1], marks[1:]):
        sol = solve_ivp(rhs, (a, b), e0, method="LSODA", rtol=1e-8, atol=1e-12, dense_output=True,
                        max_step=max((b - a) / 50, 0.05))
        pieces.append((a, b, sol.sol))
        e0 = sol.y[:, -1]

    T_zero = temps(0.0)[0]
    L_zero = np.sum(sl.dx * (1 + alpha(T_zero) * (T_zero - T_REF_C)))

    def state(times):
        times = np.atleast_1d(np.asarray(times, float))
        ext = np.empty_like(times)
        Tall = temps(times)
        for j, tj in enumerate(times):
            for a, b, f in pieces:
                if a <= tj <= b:
                    e = f(tj)
                    break
            T = Tall[j]
            k = 1 + alpha(T) * (T - T_REF_C)
            s = F_N * np.exp(e) / (sl.area * k ** 2)
            ext[j] = np.sum(sl.dx * k * np.exp(e) * (1 + s / (E_GPa(T) * 1e3))) - L_zero
        return ext, Tall

    return state


def steady_rate(mat, stress_MPa, T_C):
    return float(np.exp(mat["lnA"] + mat["n"] * np.log(stress_MPa) - mat["Q"] / (R_GAS * (T_C + 273.15))))
