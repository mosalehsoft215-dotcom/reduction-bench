"""Planar cyclic voltammetry by implicit finite volumes.

Not used for the data (every sweep is at a microdisk, see
../provenance/microdisk.py). calibrate.py uses it for the shortcut of
modelling the microdisk as a planar electrode, and as the planar half of
the "planar + steady state" shortcut.

One-electron couple O + e- = R at a planar electrode:
  * Butler-Volmer kinetics (E0, k0, alpha)
  * optional follow-up reaction of R in solution, R -> Z (rate constant kc),
    which may run back, Z -> R (kb); Z is not electroactive in the window
  * semi-infinite linear diffusion of O and R; Z diffuses like R
  * the cell: uncompensated resistance Ru in series with the interface,
    double-layer capacitance Cdl in parallel with the faradaic path, so
      i = -F A r + Cdl dE_int/dt,   E_int = E_app - i Ru
    with r the net reduction rate (mol cm^-2 s^-1)

Currents follow the IUPAC convention (anodic positive). Units: cm, s,
mol/cm^3, V, A, ohm, F.

Time stepping is backward Euler on an expanding grid; at every step the
surface flux, the interfacial potential and the capacitive current are
solved together, so iR drop and charging are exact within the step.

With kb > 0 and D_Z = D_R, the pair R, Z is carried as S = R + Z, which
only diffuses, and W = kc R - kb Z, which diffuses and decays at kc + kb;
both are then single-species problems like O.
"""

import numpy as np
from scipy.linalg import cho_solve_banded, cholesky_banded
from scipy.optimize import brentq

FARADAY = 96485.33212
GAS = 8.314462618


def expanding_grid(h0, gamma, length):
    x = [0.0]
    h = h0
    while x[-1] < length:
        x.append(x[-1] + h)
        h *= gamma
    return np.array(x)


class _Species:
    """Diffusion (+ first-order decay) of one species; the far node is held
    at the bulk value. The system for the free nodes is symmetric positive
    definite, so it is factored once."""

    def __init__(self, x, diff, decay, dt, bulk):
        h = np.diff(x)
        vol = np.empty(x.size)
        vol[0] = h[0] / 2
        vol[1:-1] = (h[:-1] + h[1:]) / 2
        vol[-1] = h[-1] / 2
        cu = diff / h
        m = x.size - 1                       # free nodes 0 .. m-1
        diag = vol[:m] / dt + decay * vol[:m]
        diag[:m] += cu[:m]
        diag[1:m] += cu[:m - 1]
        off = -cu[:m - 1]
        ab = np.zeros((2, m))
        ab[0, 1:] = off                      # upper form for cholesky_banded
        ab[1, :] = diag
        self.chol = cholesky_banded(ab)
        self.vol = vol[:m]
        self.dt = dt
        self.edge = cu[m - 1] * bulk         # far-field Dirichlet term
        self.c = np.full(m, bulk, float)
        unit = np.zeros(m)
        unit[0] = 1.0
        self.green = cho_solve_banded((self.chol, False), unit)

    def base(self):
        rhs = self.vol / self.dt * self.c
        rhs[-1] += self.edge
        return cho_solve_banded((self.chol, False), rhs)


def simulate(e_program, dt, *, e0, k0, alpha, d_ox, d_red, c_ox, c_red, area,
             ru=0.0, cdl=0.0, kc=0.0, kb=0.0, temp=298.15, h0=None, gamma=1.02):
    """Current (A) for the applied potentials e_program[n] at t_n = n dt.
    The bulk holds O at c_ox and R at c_red with no Z; with a reversible
    follow-up step (kb > 0) the bulk must hold no R."""
    if kb > 0 and c_red > 0:
        raise ValueError("a reversible follow-up step needs c_red = 0 in the bulk")
    f = FARADAY / (GAS * temp)
    t_total = dt * (len(e_program) - 1)
    dmax = max(d_ox, d_red)
    length = 6.0 * np.sqrt(dmax * t_total) + 1e-4
    if h0 is None:
        h0 = 0.02 * np.sqrt(min(d_ox, d_red) * dt)
    x = expanding_grid(h0, gamma, length)
    ox = _Species(x, d_ox, 0.0, dt, c_ox)
    if kb > 0:
        s_sp = _Species(x, d_red, 0.0, dt, 0.0)          # S = R + Z
        w_sp = _Species(x, d_red, kc + kb, dt, 0.0)      # W = kc R - kb Z
        g_ox = ox.green[0]
        g_red = (kb * s_sp.green[0] + kc * w_sp.green[0]) / (kc + kb)
    else:
        red = _Species(x, d_red, kc, dt, c_red)
        g_ox, g_red = ox.green[0], red.green[0]

    k_ir = ru * cdl / dt
    far = ru * FARADAY * area
    e_int = e_program[0]
    current = np.empty(len(e_program))
    current[0] = 0.0
    for n in range(1, len(e_program)):
        b_ox = ox.base()
        a_ox = b_ox[0]
        if kb > 0:
            b_s, b_w = s_sp.base(), w_sp.base()
            a_red = (kb * b_s[0] + b_w[0]) / (kc + kb)
        else:
            b_red = red.base()
            a_red = b_red[0]
        e_prev = e_int
        e_app = e_program[n]

        def e_of(r):
            return (e_app + far * r + k_ir * e_prev) / (1.0 + k_ir)

        def g(r):
            eta = f * (e_of(r) - e0)
            co = a_ox - r * g_ox
            cr = a_red + r * g_red
            return r - k0 * (co * np.exp(-alpha * eta) - cr * np.exp((1 - alpha) * eta))

        r_hi = a_ox / g_ox
        r_lo = -a_red / g_red
        if r_hi - r_lo <= 0:
            r = 0.0
        else:
            r = brentq(g, r_lo, r_hi, xtol=1e-22, rtol=1e-13, maxiter=200)
        e_int = e_of(r)
        ox.c = b_ox - r * ox.green
        if kb > 0:
            s_sp.c = b_s + r * s_sp.green
            w_sp.c = b_w + kc * r * w_sp.green
        else:
            red.c = b_red + r * red.green
        current[n] = -FARADAY * area * r + cdl * (e_int - e_prev) / dt
    return current


def triangle(e_start, e_switch, v, de):
    """Linear sweep e_start -> e_switch -> e_start sampled every de volts;
    returns (potentials, dt)."""
    n_half = int(round(abs(e_switch - e_start) / de))
    leg = np.linspace(e_start, e_switch, n_half + 1)
    prog = np.concatenate([leg, leg[-2::-1]])
    return prog, de / v

