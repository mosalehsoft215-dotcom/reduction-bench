"""Physics of the separator used for the simulated trial: hollow fibres, solved for all runs at once.

Shell (outside the fibres): gas from the feed end z = 0 to the retentate end z = L in plug flow at constant pressure;
its temperature goes linearly from the inlet to the outlet thermocouple value. Bores: closed at one end; the permeate
flows to the open end (the feed end z = 0 in counter-current modules, the retentate end z = L in co-current ones),
through the bores inside that potting (no permeation there) to the nozzle, where the pressure is measured; laminar
Hagen-Poiseuille flow in the potted bores and along the active length, at the local gas temperature (the potted bores
at the temperature of the gas at their end). Flux of component i through the outer skin: Q_i(T) times the
Peng-Robinson fugacity difference against the mixed permeate present in the bore at the same position, with
Q_i(T) = Q_i(35 C) exp(-E_i / R (1/T - 1/T_ref)). The bore pressure profile is iterated with damping.

Broken fibres (leak): a fibre broken at z_b no longer separates. Shell gas of the composition present at z_b enters
its open bore there and flows along the bore and through the potted bores of the open end into the permeate header:
laminar, compressible, isothermal in each slice at the local temperature,
n = N_b pi d^4 (p_shell^2 - p_nozzle^2) / (256 mu R integral(T dz)) over that path. The shell loses that gas at z_b.
"""
import numpy as np

R = 8.314462618
COMP = ("CO2", "CH4", "C2H6", "N2")
TC = np.array([304.13, 190.56, 305.32, 126.19])
PC = np.array([73.773e5, 45.992e5, 48.722e5, 33.958e5])
OMEGA = np.array([0.2239, 0.0115, 0.0995, 0.0372])
KIJ = np.array([[0.0, 0.0919, 0.1322, -0.0170],
                [0.0919, 0.0, 0.0026, 0.0311],
                [0.1322, 0.0026, 0.0, 0.0515],
                [-0.0170, 0.0311, 0.0515, 0.0]])
GPU = 3.348e-10                    # mol / (m2 s Pa)
T_REF = 308.15
MU = np.array([1.53e-5, 1.13e-5, 0.95e-5, 1.80e-5])       # pure-gas viscosities near 35 C, Pa s
MW = np.array([44.01, 16.04, 30.07, 28.01])
_KAPPA = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA ** 2
_B = 0.07780 * R * TC / PC


def fugacity(y, p, t):
    """fugacities (Pa) of the components of a gas mixture; y (..., 4); p and t of shape y.shape[:-1]"""
    p = np.asarray(p, float)[..., None]
    t = np.asarray(t, float)[..., None]
    a = 0.45724 * R ** 2 * TC ** 2 / PC * (1 + _KAPPA * (1 - np.sqrt(t / TC))) ** 2
    aij = np.sqrt(a[..., :, None] * a[..., None, :]) * (1 - KIJ)
    sum_a = np.einsum("...ij,...j->...i", aij, y)
    am = (y * sum_a).sum(-1, keepdims=True)
    bm = (y * _B).sum(-1, keepdims=True)
    A = am * p / (R * t) ** 2
    Bm = bm * p / (R * t)
    z = np.ones_like(A)
    for _ in range(12):                 # vapour root of the cubic
        f = z ** 3 - (1 - Bm) * z ** 2 + (A - 3 * Bm ** 2 - 2 * Bm) * z - (A * Bm - Bm ** 2 - Bm ** 3)
        z = z - f / (3 * z ** 2 - 2 * (1 - Bm) * z + (A - 3 * Bm ** 2 - 2 * Bm))
    s2 = np.sqrt(2.0)
    ln_phi = (_B / bm * (z - 1) - np.log(z - Bm)
              - A / (2 * s2 * Bm) * (2 * sum_a / am - _B / bm) * np.log((z + (1 + s2) * Bm) / (z + (1 - s2) * Bm)))
    return y * p * np.exp(ln_phi)


def wilke(y):
    """viscosity of a gas mixture (Pa s) by Wilke's rule, y (..., 4)"""
    phi = (1 + np.sqrt(MU[:, None] / MU[None, :]) * (MW[None, :] / MW[:, None]) ** 0.25) ** 2 / np.sqrt(8 * (1 + MW[:, None] / MW[None, :]))
    return (y * MU / np.einsum("...j,ij->...i", y, phi)).sum(-1)


def simulate(runs, geom, q_ref_gpu, e_kj, n=100, sweeps=200, tol=1e-12, outlet="feed", leak=None):
    """retentate (shell outflow), permeate (bores) and leak molar flows (mol/s, arrays (runs, 4)).
    runs: feed (runs, 4) mol/s, p_feed, p_perm (Pa), t_feed, t_ret (K), arrays (runs,).
    geom: fibres (open fibres), od, id, active, potting (m). Permeate viscosity from Wilke's rule.
    outlet: "feed" (counter-current) or "ret" (co-current), the end where the bores are open.
    leak: None or dict(z = distance of the break from the feed-end face of the active length (m), count = fibres)."""
    feed = np.asarray(runs["feed"], float)
    nr = feed.shape[0]
    length = geom["active"]
    zf = np.linspace(0.0, 1.0, n + 1)
    kb = None
    if leak is not None:
        sb = leak["z"] / length
        zf = np.unique(np.round(np.concatenate([zf, [sb]]), 12))
        kb = int(np.argmin(np.abs(zf - sb)))
    nz = zf.size
    dz = np.diff(zf) * length
    fibres = geom["fibres"] - (leak["count"] if leak is not None else 0)
    da = fibres * np.pi * geom["od"]
    tf, tr = runs["t_feed"], runs["t_ret"]
    t = tf[:, None] + (tr - tf)[:, None] * zf
    q = np.asarray(q_ref_gpu) * GPU * np.exp(-np.asarray(e_kj) * 1e3 / R * (1 / t[..., None] - 1 / T_REF))
    pf = np.repeat(runs["p_feed"][:, None], nz, 1)
    yb = np.repeat((feed / feed.sum(1, keepdims=True))[:, None, :], nz, 1)
    pb = np.repeat(runs["p_perm"][:, None], nz, 1)
    g = np.zeros((nr, nz, 4))
    hp = 256 * R / (np.pi * geom["id"] ** 4 * fibres)
    pot = geom.get("potting", 0.0)
    lk = np.zeros((nr, 4))
    if leak is not None:        # integral of T over the open bore from the break to the nozzle
        if outlet == "feed":
            tint = length * (sb * tf + (tr - tf) * sb ** 2 / 2) + pot * tf
        else:
            tint = length * ((1 - sb) * tf + (tr - tf) * (1 - sb ** 2) / 2) + pot * tr
        cl = leak["count"] * np.pi * geom["id"] ** 4 * (runs["p_feed"] ** 2 - runs["p_perm"] ** 2) / (256 * R * tint)
    for sweep in range(sweeps):
        fb = fugacity(yb, pb, t)
        flow = np.zeros((nr, nz, 4))
        flow[:, 0] = feed
        flux = np.zeros((nr, nz, 4))
        for k in range(nz):
            x = flow[:, k] / flow[:, k].sum(1, keepdims=True)
            flux[:, k] = q[:, k] * (fugacity(x, pf[:, k], t[:, k]) - fb[:, k]) * da
            if k == kb:         # shell gas of the local composition leaves through the broken bores
                lk = (cl / wilke(x))[:, None] * x
                flow[:, k] = flow[:, k] - lk
            if k < nz - 1:
                pred = flow[:, k] - flux[:, k] * dz[k]
                xp = pred / pred.sum(1, keepdims=True)
                f2 = q[:, k + 1] * (fugacity(xp, pf[:, k + 1], t[:, k + 1]) - fb[:, k + 1]) * da
                flow[:, k + 1] = flow[:, k] - 0.5 * (flux[:, k] + f2) * dz[k]
        seg = 0.5 * (flux[:, 1:] + flux[:, :-1]) * dz[None, :, None]
        g_new = np.zeros_like(g)
        if outlet == "feed":    # bores closed at z = L, permeate collected towards z = 0
            g_new[:, :-1] = np.cumsum(seg[:, ::-1], axis=1)[:, ::-1]
        else:                   # bores closed at z = 0, permeate collected towards z = L
            g_new[:, 1:] = np.cumsum(seg, axis=1)
        tot = g_new.sum(-1, keepdims=True)
        y_new = np.where(tot > 0, g_new / np.where(tot > 0, tot, 1.0), flux / flux.sum(-1, keepdims=True))
        gm = 0.5 * (tot[:, 1:, 0] + tot[:, :-1, 0])
        tm = 0.5 * (t[:, 1:] + t[:, :-1])
        mu = wilke(y_new)
        mm = 0.5 * (mu[:, 1:] + mu[:, :-1])
        inc = hp * mm * gm * tm * dz[None, :]
        # the whole permeate of a fibre passes the potted bores of the open end, at the gas temperature of that end
        if outlet == "feed":
            p_edge2 = runs["p_perm"] ** 2 + hp * mu[:, 0] * tot[:, 0, 0] * tf * pot
            p2 = p_edge2[:, None] + np.concatenate([np.zeros((nr, 1)), np.cumsum(inc, 1)], 1)
        else:
            p_edge2 = runs["p_perm"] ** 2 + hp * mu[:, -1] * tot[:, -1, 0] * tr * pot
            p2 = p_edge2[:, None] + np.concatenate([np.cumsum(inc[:, ::-1], 1)[:, ::-1], np.zeros((nr, 1))], 1)
        p_new = np.sqrt(p2)
        change = np.abs(g_new - g).max() / g_new.max()
        g = g_new
        yb = 0.5 * yb + 0.5 * y_new
        pb = 0.5 * pb + 0.5 * p_new
        if change < tol:
            break
    return flow[:, -1], (g[:, 0] if outlet == "feed" else g[:, -1]), lk


CP = np.array([[19.80, 7.344e-2, -5.602e-5, 1.715e-8],     # ideal-gas heat capacity, J/(mol K), T in K
               [19.25, 5.213e-2, 1.197e-5, -1.132e-8],
               [5.409, 1.781e-1, -6.938e-5, 8.713e-9],
               [31.15, -1.357e-2, 2.680e-5, -1.168e-8]])


def enthalpy(y, p, t):
    """molar enthalpy (J/mol) of a gas mixture relative to the ideal gas at 298.15 K: ideal part + Peng-Robinson
    departure; y (4,), scalars p (Pa), t (K)"""
    t0 = 298.15
    ig = sum(y[i] * (CP[i, 0] * (t - t0) + CP[i, 1] / 2 * (t ** 2 - t0 ** 2) + CP[i, 2] / 3 * (t ** 3 - t0 ** 3)
                     + CP[i, 3] / 4 * (t ** 4 - t0 ** 4)) for i in range(4))
    sq = np.sqrt(t / TC)
    alpha = (1 + _KAPPA * (1 - sq)) ** 2
    dalpha = -_KAPPA * (1 + _KAPPA * (1 - sq)) / np.sqrt(t * TC)
    ac = 0.45724 * R ** 2 * TC ** 2 / PC
    a, da = ac * alpha, ac * dalpha
    root = np.sqrt(np.outer(a, a))
    am = float(y @ (root * (1 - KIJ)) @ y)
    dam = float(y @ ((np.outer(da, a) + np.outer(a, da)) / (2 * root) * (1 - KIJ)) @ y)
    bm = float(y @ _B)
    A, Bm = am * p / (R * t) ** 2, bm * p / (R * t)
    z = 1.0
    for _ in range(30):
        f = z ** 3 - (1 - Bm) * z ** 2 + (A - 3 * Bm ** 2 - 2 * Bm) * z - (A * Bm - Bm ** 2 - Bm ** 3)
        z -= f / (3 * z ** 2 - 2 * (1 - Bm) * z + (A - 3 * Bm ** 2 - 2 * Bm))
    s2 = np.sqrt(2.0)
    dep = R * t * (z - 1) + (t * dam - am) / (2 * s2 * bm) * np.log((z + (1 + s2) * Bm) / (z + (1 - s2) * Bm))
    return ig + dep


def throttle(y, p1, t1, p2):
    """temperature after an isenthalpic expansion from (p1, t1) to p2 (bisection)"""
    h = enthalpy(y, p1, t1)
    lo, hi = t1 - 60.0, t1 + 5.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if enthalpy(y, p2, mid) < h else (lo, mid)
    return 0.5 * (lo + hi)
