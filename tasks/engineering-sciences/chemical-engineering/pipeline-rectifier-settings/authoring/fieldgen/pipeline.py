"""Impressed-current cathodic protection of a coated pipeline as a nonlinear leaky transmission line.

Nodes every dx metres; electrically insulated ends. Unknown: pipe metal potential V (vs remote earth).
Ground beds are point current sources at the surface: soil potential beside the pipe
    Vs(x) = sum_k I_k rho_k / (2 pi sqrt((x - x_bed_k)^2 + offset_k^2)).
Pipe-to-adjacent-soil voltage U = V - Vs; drive = -U.
Pickup per node through coating defects: J = G (drive - eta), G = g(x) dx,
interface (steel/soil) law: J = K (exp(eta / b) - 1), K = k(x) dx.
A shorted casing at a crossing drains Gc * drive at its node through the short (bare, grounded casing steel
with no polarization of the pipe there).
Potentials vs CSE: E_on = E_nat + U; instant-off (IR-free, polarized) E_off = E_nat - eta.
Rectifier k drains I_k from the pipe at its drain node. Pipe longitudinal resistance r per metre.
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def eta_from_drive(drive, G, K, b, iters=70):
    drive = np.maximum(drive, 0.0)
    lo = np.zeros_like(drive)
    hi = drive.copy()
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        f = G * (drive - mid) - K * np.expm1(np.minimum(mid / b, 60.0))
        lo = np.where(f > 0, mid, lo)
        hi = np.where(f > 0, hi, mid)
    return 0.5 * (lo + hi)


class Pipeline:
    def __init__(self, length_m, dx, r_per_m, g, k, b, e_nat, beds, points=()):
        self.x = np.arange(0.0, length_m + 1e-9, dx)
        self.n = len(self.x)
        self.dx = dx
        self.r = r_per_m
        self.b = b
        self.G = (g(self.x) if callable(g) else np.full(self.n, g)) * dx
        self.K = (k(self.x) if callable(k) else np.full(self.n, k)) * dx
        self.Gc = np.zeros(self.n)
        for xp, Gp in points:
            self.Gc[int(round(xp / dx))] += Gp
        self.e_nat = e_nat(self.x) if callable(e_nat) else np.full(self.n, e_nat)
        self.beds = beds
        cond = 1.0 / (r_per_m * dx)
        main = np.zeros(self.n)
        main[:-1] += cond
        main[1:] += cond
        off = -cond * np.ones(self.n - 1)
        self.lap = diags([off, main, off], [-1, 0, 1], format="csc")

    def soil(self, currents):
        vs = np.zeros(self.n)
        for bed, I in zip(self.beds, currents):
            d = np.sqrt((self.x - bed["x_bed"]) ** 2 + bed["offset"] ** 2)
            vs += bed["rho"] * I / (2 * np.pi * d)
        return vs

    def solve(self, currents, tol=1e-11, max_iter=100):
        vs = self.soil(currents)
        inj = np.zeros(self.n)
        for bed, I in zip(self.beds, currents):
            inj[int(round(bed["x_drain"] / self.dx))] -= I
        V = -0.3 - vs
        for it in range(max_iter):
            drive = vs - V
            eta = eta_from_drive(drive, self.G, self.K, self.b)
            J = self.G * (np.maximum(drive, 0) - eta) + self.Gc * drive
            res = self.lap @ V - (J + inj)
            fprime = self.K / self.b * np.exp(np.minimum(eta / self.b, 60.0))
            dJdrive = self.G * fprime / (self.G + fprime) * (drive > 0) + self.Gc
            jac = self.lap + diags(dJdrive, 0, format="csc")
            step = spsolve(jac, -res)
            V = V + step
            if np.max(np.abs(step)) < tol:
                break
        drive = vs - V
        eta = eta_from_drive(drive, self.G, self.K, self.b)
        J = self.G * (np.maximum(drive, 0) - eta) + self.Gc * drive
        flow = np.cumsum(J + inj)  # current along the pipe between node i and i+1, toward increasing x
        return dict(V=V, Vs=vs, e_on=self.e_nat + (V - vs), e_off=self.e_nat - eta, eta=eta, J=J,
                    pipe_current=np.concatenate([flow[:-1], [0.0]]), iterations=it + 1)


def margins(p, currents, low=-0.85, high=-1.15):
    """(protection margin, overprotection margin) in volts; both >= 0 when the design is acceptable."""
    eo = p.solve(currents)["e_off"]
    return low - eo.max(), eo.min() - high
