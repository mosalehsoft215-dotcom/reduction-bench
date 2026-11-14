"""Axial load of a roof tie (strand backstay or solid bar) before and after its retensioning, from two ambient
vibration records.

    python3 tie_loads.py <tie folder> <output file>

1. Backstays: spectrum of each record by Welch (Hann, sixth-length segments); peaks with a prominence of at least
   12 dB, each located by a Lorentzian on a sloping floor fitted to the thirteen spectral lines around it.  Peaks of
   the structure and stray peaks stay in: the ladder match charges an unmatched peak a fixed cap.
2. In-plane model of the tie (Hermite beam elements on the unit length): bending stiffness EI, axial load T, mass per
   metre, deflection held at both ends, end slopes free (backstays), held, or restrained by a rotational spring (a
   bar end on a bedded nut, c = k_end Lv / EI), the sag-extensibility stiffness of a flat-sagged inclined tie, and a
   factor on the mass per metre.
3. Bars (EI = E pi d^4 / 64, vertical, no sag, vibrating over the clear length Lv between the block faces).  The three
   pickups give the cross-spectral matrix (Welch, Hann, eighth-length segments); peaks of its trace located by the
   Lorentzian, and at each peak the operating shape, the principal eigenvector of the matrix summed over three lines,
   turned real.  Frequencies alone cannot tell a stiffer bar on softer beds from a slacker bar on harder ones; the
   shape near each block shows how far the end turns.  The dimensionless ladder Omega_k(xi; c_lo, c_hi), xi = Lv
   sqrt(T / EI), is scanned for each restraint pair of a coarse grid (0.3, 2, 12, held at each end); the six best
   are polished by least squares over (T_pre, T_post, log c_lo, log c_hi) on matched (peak, mode) pairs of both
   records, residuals the log frequencies and the shape components, pairs re-made over five rounds; the lowest
   total wins.  A bar with plain nuts comes out with both c at the upper bound (held).
4. Backstays (EI of the strand bundle, E A from the steel area, sheet length, ends free to rotate).  A scan of T per
   record with the sheet mass, then a least-squares fit of (T_pre / kappa, T_post / kappa) at fixed kappa, kappa
   being the ratio of the vibrating mass per metre to the sheet value.  The ladders only fix r = T / (kappa m); on a
   short, taut backstay the first symmetric mode does not move with the weight at all, so kappa comes from the shim
   change instead: the shims draw that much strand through the ground anchor, and the unstressed strand length
   between the anchorages of a flat-sagged elastic tie on its chord, referred to one temperature, is
   S0 = (L + q^2 L^3 / (24 T^2)) / ((1 + T / EA)(1 + 12e-6 (t - 20))), q = kappa m g cos(theta), t the steel
   temperature of the record (a backstay warmed by 8 K between its records grows by 1 mm per 10 m, as much as a
   few hundred kN of load change on a taut one).  With T = kappa m r the sag term does not depend on kappa, so S0(pre) - S0(post) = shim
   is solved for kappa (bracketing root), the ladders are refitted at that kappa (the sag term of a slack backstay
   moves the first symmetric mode) and the two steps alternate four times.
"""
import csv
import json
import os
import sys
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")            # small matrices: one BLAS thread is the fast setting on two CPUs

import numpy as np  # noqa: E402
from scipy.linalg import eigh
from scipy.optimize import curve_fit, least_squares
from scipy.signal import csd, find_peaks, welch

G = 9.80665
CAP = 0.02
STRAND_D = 0.0157
KMAX = 12
ALPHA = 12e-6                                               # thermal expansion of steel per kelvin
SAG_FLOOR = 0.5
BED_GRID = (0.3, 2.0, 12.0, np.inf)                          # k_end Lv / EI of the start grid (inf: held)
NSTART = 6
SIG_F, SIG_U = 0.0015, 0.01                                # peak frequency (relative) and shape component scatter
XI_COARSE = np.geomspace(3.0, 120.0, 30)
XI_FINE = np.geomspace(3.0, 120.0, 260)
_UNIT = {}


def unit_mesh(nel):
    """bending, geometric and mass matrices and the uniform-load vector on the unit length (dofs w, dw/dx)"""
    if nel not in _UNIT:
        h = 1.0 / nel
        nd = 2 * (nel + 1)
        Kb, Kg, M, c = np.zeros((nd, nd)), np.zeros((nd, nd)), np.zeros((nd, nd)), np.zeros(nd)
        hh = h * h
        kb = np.array([[12, 6 * h, -12, 6 * h], [6 * h, 4 * hh, -6 * h, 2 * hh],
                       [-12, -6 * h, 12, -6 * h], [6 * h, 2 * hh, -6 * h, 4 * hh]]) / h ** 3
        kg = np.array([[36, 3 * h, -36, 3 * h], [3 * h, 4 * hh, -3 * h, -hh],
                       [-36, -3 * h, 36, -3 * h], [3 * h, -hh, -3 * h, 4 * hh]]) / (30 * h)
        me = np.array([[156, 22 * h, 54, -13 * h], [22 * h, 4 * hh, 13 * h, -3 * hh],
                       [54, 13 * h, 156, -22 * h], [-13 * h, -3 * hh, -22 * h, 4 * hh]]) * h / 420
        for e in range(nel):
            s = slice(2 * e, 2 * e + 4)
            Kb[s, s] += kb
            Kg[s, s] += kg
            M[s, s] += me
            c[s] += np.array([h / 2, hh / 12, h / 2, -hh / 12])
        _UNIT[nel] = (Kb, Kg, M, c)
    return _UNIT[nel]


def station_vector(nel, s):
    """Hermite interpolation row for the deflection at station s (0-1) of the unit length"""
    h = 1.0 / nel
    e = min(int(s * nel), nel - 1)
    t = s * nel - e
    row = np.zeros(2 * (nel + 1))
    row[2 * e:2 * e + 4] = (1 - 3 * t * t + 2 * t ** 3, h * (t - 2 * t * t + t ** 3), 3 * t * t - 2 * t ** 3, h * (t ** 3 - t * t))
    return row


def omega(nel, xi2, held, eps=0.0, s=0.5, sag=0.0, kmax=KMAX):
    """dimensionless circular frequencies of  w'''' - xi2 w'' + sag (int w) = Omega^2 (1 + eps delta(x - s)) w"""
    Kb, Kg, M, c = unit_mesh(nel)
    K = Kb + xi2 * Kg
    nd = K.shape[0]
    if sag:
        K = K + sag * np.outer(c, c)
    if eps:
        row = station_vector(nel, s)
        M = M + eps * np.outer(row, row)
    drop = {True: (0, 1, nd - 2, nd - 1), "top": (0, nd - 2, nd - 1), False: (0, nd - 2)}[held]
    keep = [j for j in range(nd) if j not in drop]
    w = eigh(K[np.ix_(keep, keep)], M[np.ix_(keep, keep)], eigvals_only=True, subset_by_index=[0, kmax - 1])
    return np.sqrt(np.maximum(w, 0.0))


def stay_ladder(T, L, EI, EA, m, theta, kappa=1.0, nel=100, kmax=26):
    mk = m * kappa
    q = mk * G * np.cos(np.radians(theta))
    sag = q * L * L / (8 * T)
    Le = L * (1 + 8 * (sag / L) ** 2)
    stiff = EA / Le * (q / T) ** 2 * L ** 5 / EI           # rank-one term in units of EI / L^3 on the unit length
    return omega(nel, L * L * T / EI, False, sag=stiff, kmax=kmax) / (2 * np.pi) * np.sqrt(EI / (mk * L ** 4))


def lorentz(x, A, f0, g, c0, c1):
    return A / ((x - f0) ** 2 + g * g) + c0 + c1 * (x - f0)


def peaks(acc, fs, segments):
    f, p = welch(acc - acc.mean(), fs=fs, window="hann", nperseg=len(acc) // segments)
    return peaks_of(f, p)


def peaks_of(f, p):
    db = 10 * np.log10(p + 1e-30)
    idx, _ = find_peaks(db, prominence=12.0)
    df = f[1] - f[0]
    out = []
    for i in idx:
        if f[i] < 0.12 or i < 7 or i + 7 >= len(f):
            continue
        x, y = f[i - 6:i + 7], p[i - 6:i + 7]
        try:
            prm = curve_fit(lorentz, x, y, p0=(p[i] * df * df, f[i], df, 0.0, 0.0), maxfev=6000)[0]
            out.append(prm[1] if abs(prm[1] - f[i]) < 2 * df else f[i])
        except RuntimeError:
            out.append(f[i])
    return np.array(out)


def cost(obs, model):
    band = model[(model > obs.min() * 0.97) & (model < obs.max() * 1.03)]
    d = np.abs(np.log(obs[:, None] / model[None, :])).min(axis=1)
    hit = np.abs(np.log(band[:, None] / obs[None, :])).min(axis=1) if len(band) else np.zeros(0)
    return np.sum(np.minimum(d, CAP) ** 2) + 0.35 * CAP ** 2 * np.sum(hit > CAP)


def costs(obs, ladders):
    """cost() for many model ladders at once (rows of ladders)"""
    d = np.abs(np.log(obs[None, :, None] / ladders[:, None, :])).min(axis=2)
    band = (ladders > obs.min() * 0.97) & (ladders < obs.max() * 1.03)
    hit = np.abs(np.log(ladders[:, :, None] / obs[None, None, :])).min(axis=2) > CAP
    return np.sum(np.minimum(d, CAP) ** 2, axis=1) + 0.35 * CAP ** 2 * np.sum(band & hit, axis=1)


def pairs(obs, model, tol):
    k = np.abs(np.log(obs[:, None] / model[None, :])).argmin(axis=1)
    ok = np.abs(np.log(obs / model[k])) < tol
    return obs[ok], k[ok]


def polish(x0, scale, bounds, ladder_of, obs_sets, res, extra=None, rounds=4, tol=0.012, wide=False):
    """least squares on matched (peak, mode) pairs, residuals in Hz over the line spacing; pairs re-made each round.
    wide: two extra opening rounds with a looser match (4 %, then 2.5 %) and relative residuals, so that a start a few
    per cent off still pairs every peak with its own mode before the Hz-weighted rounds"""
    x = np.clip(np.array(x0, float), bounds[0], bounds[1])
    plan = ([(0.04, True), (0.025, True)] if wide else []) + [(tol, False)] * rounds
    for tl, rel in plan:
        sets = [pairs(obs, ladder_of(x, i), tl) for i, obs in enumerate(obs_sets)]

        def fun(z, rel=rel):
            if rel:
                r = [np.log(o / ladder_of(z, i)[k]) / 1e-3 for i, (o, k) in enumerate(sets)]
            else:
                r = [(o - ladder_of(z, i)[k]) / res for i, (o, k) in enumerate(sets)]
            if extra is not None:
                r.append(extra(z))
            return np.concatenate(r)
        x = least_squares(fun, x0=x, x_scale=scale, bounds=bounds, diff_step=1e-4).x
    return x


def load(folder):
    sheet = json.loads((folder / "tie_sheet.json").read_text())
    recs = {}
    for tag, name in sheet["files"].items():
        with open(folder / name) as fh:
            rows = list(csv.reader(fh))
        cols = np.array(rows[1:], float)
        recs[tag] = cols[:, 0] if cols.shape[1] == 1 else cols
    return sheet, recs


def bar_shapes(nel, xi2, beds, stations, kmax=KMAX):
    """dimensionless circular frequencies of  w_xxxx - xi2 w_xx = Omega^2 w  on the unit length, deflection held at both
    ends, end slopes restrained by beds = (c_lo, c_hi) = k_end Lv / EI (inf: held); and the deflection of each mode at
    the stations, as unit vectors"""
    Kb, Kg, M, _ = unit_mesh(nel)
    K = Kb + xi2 * Kg
    nd = K.shape[0]
    drop = [0, nd - 2]
    for dof, c in zip((1, nd - 1), beds):
        if np.isinf(c):
            drop.append(dof)
        else:
            K[dof, dof] += c
    keep = [j for j in range(nd) if j not in drop]
    w, v = eigh(K[np.ix_(keep, keep)], M[np.ix_(keep, keep)], subset_by_index=[0, kmax - 1])
    full = np.zeros((nd, kmax))
    full[keep] = v
    at = np.array([station_vector(nel, x) for x in stations]) @ full
    return np.sqrt(np.maximum(w, 0.0)), at / np.linalg.norm(at, axis=0)


def bar_peaks(acc, fs, segments):
    """peaks of the summed auto-spectra of the three pickups, each located by a Lorentzian; at each peak the operating
    shape: principal eigenvector of the 3 x 3 cross-spectral matrix summed over the three lines around the peak,
    turned real (phase of its largest entry), as a unit vector"""
    x = acc - acc.mean(axis=0)
    nper = len(x) // segments
    f, _ = welch(x[:, 0], fs=fs, window="hann", nperseg=nper)
    G = np.zeros((len(f), 3, 3), complex)
    for a in range(3):
        for b in range(a, 3):
            G[:, a, b] = csd(x[:, a], x[:, b], fs=fs, window="hann", nperseg=nper)[1]
            G[:, b, a] = np.conj(G[:, a, b])
    tr = np.real(np.einsum("kii->k", G))
    freqs, shapes = [], []
    for fk in peaks_of(f, tr):
        i = int(np.argmin(np.abs(f - fk)))
        val, vec = np.linalg.eigh(G[i - 1:i + 2].sum(axis=0))
        u = vec[:, -1]
        u = np.real(u * np.exp(-1j * np.angle(u[np.argmax(np.abs(u))])))
        freqs.append(fk)
        shapes.append(u / np.linalg.norm(u))
    return np.array(freqs), np.array(shapes).reshape(-1, 3)


def solve_bar(sheet, rb, ra, fs, shapes=True, beds=True, segments=8):
    """T before and after with the two end restraints shared by both records: a scan of xi per record for each
    restraint pair of a coarse grid, then least squares on matched (peak, mode) pairs over (T_pre, T_post, log c_lo,
    log c_hi), residuals the log peak frequencies in units of SIG_F and the three operating-shape
    components in units of SIG_U (their scatter between records), pairs re-made each round; the restraint pair with the lowest total wins.
    shapes=False drops the shape residuals (frequencies only); beds=False holds both ends (plain nuts)"""
    m = sheet["unit_mass_kg_per_m"]
    Lv = sheet["clear_length_m"]
    d = sheet["make_up"]["bar_diameter_mm"] / 1000
    EI = sheet["young_modulus_mpa"] * 1e6 * np.pi * d ** 4 / 64
    hts = sheet["pickups"]["metres_above_lower_block_face"]
    st = [hts[k] / Lv for k in ("a_lo_mg", "a_mid_mg", "a_hi_mg")]
    obs = [bar_peaks(r, fs, segments) for r in (rb, ra)]
    scale = np.sqrt(EI / (m * Lv ** 4)) / (2 * np.pi)
    grid = [(c1, c2) for c1 in BED_GRID for c2 in BED_GRID] if beds else [(np.inf, np.inf)]
    starts = []
    for c1, c2 in grid:
        table = np.array([bar_shapes(36, xi * xi, (c1, c2), st)[0] for xi in XI_COARSE])
        fine = np.column_stack([np.interp(np.log(XI_FINE), np.log(XI_COARSE), table[:, k]) for k in range(KMAX)])
        Ts, tot = [], 0.0
        for fo, _ in obs:
            cs = costs(fo, fine * scale)
            j = int(np.argmin(cs))
            tot += cs[j]
            Ts.append((XI_FINE[j] / Lv) ** 2 * EI)
        starts.append((tot, c1, c2, Ts))
    starts.sort(key=lambda b: b[0])
    picks = starts[:NSTART]
    winner = None
    for _, c1, c2, Ts in picks:
        z0 = [Ts[0], Ts[1]] + ([np.log(min(c1, 1e4)), np.log(min(c2, 1e4))] if beds else [])
        lo = [1e3, 1e3] + ([np.log(0.05), np.log(0.05)] if beds else [])
        hi = [1e8, 1e8] + ([np.log(1e4), np.log(1e4)] if beds else [])

        def model(z, i):
            b = (np.exp(z[2]), np.exp(z[3])) if beds else (np.inf, np.inf)
            om, at = bar_shapes(60, Lv * Lv * max(z[i], 1.0) / EI, b, st)
            return om * scale, at
        z = np.clip(np.array(z0, float), lo, hi)
        for tl in (0.04, 0.025, 0.012, 0.012, 0.012):
            sets = []
            for i, (fo, uo) in enumerate(obs):
                fm, _ = model(z, i)
                k = np.abs(np.log(fo[:, None] / fm[None, :])).argmin(axis=1)
                ok = np.abs(np.log(fo / fm[k])) < tl
                sets.append((fo[ok], uo[ok], k[ok]))

            def fun(y):
                r = []
                for i, (fo, uo, k) in enumerate(sets):
                    fm, am = model(y, i)
                    r.append(np.log(fo / fm[k]) / SIG_F)
                    if shapes:
                        v = am[:, k].T
                        sgn = np.sign(np.sum(v * uo, axis=1, keepdims=True) + 1e-12)
                        r.append(((uo - sgn * v) / SIG_U).ravel())
                return np.concatenate(r)
            sol = least_squares(fun, x0=z, x_scale=[z[0], z[1]] + [1.0] * (len(z) - 2), bounds=(lo, hi), diff_step=1e-4)
            z = sol.x
        total = 2 * sol.cost + 25.0 * sum(len(o[0]) - len(t[0]) for o, t in zip(obs, sets))
        if winner is None or total < winner[0]:
            winner = (total, z)
    return winner[1][0], winner[1][1]


def solve_stay(sheet, pb, pa, res, mass_factor):
    m = sheet["unit_mass_kg_per_m"]
    L = sheet["working_point_length_m"]
    E = sheet["young_modulus_mpa"] * 1e6
    EI = sheet["make_up"]["strand_count"] * E * np.pi * STRAND_D ** 4 / 64
    EA = sheet["steel_area_mm2"] * 1e-6 * E
    th = sheet["chord_slope_deg"]
    first = []
    for obs in (pb, pa):
        spacing = np.median(np.diff(np.sort(obs)))
        T0 = 4 * m * L * L * spacing ** 2
        for lo, hi, num in ((0.15, 6.0, 160), (0.96, 1.04, 41)):
            grid = np.geomspace(lo * T0, hi * T0, num)
            T0 = grid[int(np.argmin([cost(obs, stay_ladder(T, L, EI, EA, m, th, nel=70)) for T in grid]))]
        first.append(T0)
    lad = lambda z, i: stay_ladder(z[i] * z[2], L, EI, EA, m, th, z[2])     # z = (T_pre / kappa, T_post / kappa, kappa)
    q = m * G * np.cos(np.radians(th))
    lam2 = max((q * L / T) ** 2 * EA / T for T in first)    # Irvine's parameter with the sheet mass
    free = mass_factor and lam2 > SAG_FLOOR
    winner = None
    for kappa in ((1.0, 1.05, 1.10, 1.15, 1.20) if free else (1.0,)):
        lo = [1e3, 1e3, 0.85 if free else 1.0 - 1e-9]
        hi = [1e9, 1e9, 1.40 if free else 1.0 + 1e-9]
        x = polish([first[0], first[1], kappa], [first[0], first[1], 0.05], (lo, hi), lad, (pb, pa), res,
                   extra=lambda z: np.array([(z[2] - 1.0) / 0.25]))
        total = cost(pb, lad(x, 0)) + cost(pa, lad(x, 1)) + 2e-4 * ((x[2] - 1.0) / 0.25) ** 2
        if winner is None or total < winner[0]:
            winner = (total, x)
    x = winner[1]
    return x[0] * x[2], x[1] * x[2]


def shim_kappa(sheet, rb, ra, m, L, EA, th, thermal=True):
    """kappa from the shim change, given r = T / (kappa m) in each record"""
    from scipy.optimize import brentq
    shim = sheet["retensioning"]["shim_change_mm"] / 1e3
    sag = (G * np.cos(np.radians(th))) ** 2 * L ** 3 / 24                 # q^2 / T^2 = (g cos(theta) / r)^2
    temp = sheet["steel_temperature_c"]
    hb, ha = ((1 + ALPHA * (temp["pre"] - 20.0)), (1 + ALPHA * (temp["post"] - 20.0))) if thermal else (1.0, 1.0)

    def gap(k):
        tb, ta = k * m * rb, k * m * ra
        return (L + sag / rb ** 2) / ((1 + tb / EA) * hb) - (L + sag / ra ** 2) / ((1 + ta / EA) * ha) - shim
    try:
        return brentq(gap, 0.6, 2.0, xtol=1e-7)
    except ValueError:
        return 1.0


def solve_stay_shim(sheet, pb, pa, res, floor=False, thermal=True):
    m = sheet["unit_mass_kg_per_m"]
    L = sheet["working_point_length_m"]
    E = sheet["young_modulus_mpa"] * 1e6
    EI = sheet["make_up"]["strand_count"] * E * np.pi * STRAND_D ** 4 / 64
    EA = sheet["steel_area_mm2"] * 1e-6 * E
    th = sheet["chord_slope_deg"]
    first = []
    for obs in (pb, pa):
        spacing = np.median(np.diff(np.sort(obs)))
        T0 = 4 * m * L * L * spacing ** 2
        for lo, hi, num in ((0.15, 6.0, 160), (0.96, 1.04, 41)):
            grid = np.geomspace(lo * T0, hi * T0, num)
            T0 = grid[int(np.argmin([cost(obs, stay_ladder(T, L, EI, EA, m, th, nel=70)) for T in grid]))]
        first.append(T0)
    kappa, z = 1.0, list(first)
    for _ in range(4):
        lad = lambda x, i, k=kappa: stay_ladder(x[i] * k, L, EI, EA, m, th, k)          # x = (T_pre / kappa, T_post / kappa)
        z = polish(z, [z[0], z[1]], ([1e3, 1e3], [1e9, 1e9]), lad, (pb, pa), res)
        kappa = shim_kappa(sheet, z[0] / m, z[1] / m, m, L, EA, th, thermal)
        kappa = max(kappa, 1.0) if floor else kappa
    lad = lambda x, i, k=kappa: stay_ladder(x[i] * k, L, EI, EA, m, th, k)
    z = polish(z, [z[0], z[1]], ([1e3, 1e3], [1e9, 1e9]), lad, (pb, pa), res)
    return z[0] * kappa, z[1] * kappa


def solve(folder, shapes=True, beds=True, mass_factor=True, shim=True, floor=False, thermal=True):
    """the reference is solve(folder); the keywords switch parts off for the route tables of authoring/evidence:
    shapes=False fits a bar on its frequencies alone, beds=False holds both bar ends; mass_factor / shim (backstays:
    filling, shim change), floor keeps the mass factor of the shim route at one or above (a filling only adds mass),
    thermal=False leaves the steel temperature out of the shim balance"""
    sheet, recs = load(folder)
    fs = sheet["samples_per_second"]
    if "bar_diameter_mm" in sheet["make_up"]:
        Tb, Ta = solve_bar(sheet, recs["pre"], recs["post"], fs, shapes, beds)
        return Tb / 1e3, Ta / 1e3
    segments = 6
    pb, pa = peaks(recs["pre"], fs, segments), peaks(recs["post"], fs, segments)
    res = segments * fs / len(recs["pre"])                  # spectral line spacing of the Welch segments
    if shim and mass_factor:
        Tb, Ta = solve_stay_shim(sheet, pb, pa, res, floor, thermal)
    else:
        Tb, Ta = solve_stay(sheet, pb, pa, res, mass_factor)
    return Tb / 1e3, Ta / 1e3


def main():
    tb, ta = solve(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps({"load_pre_kn": round(tb, 3), "load_post_kn": round(ta, 3)}) + "\n")


if __name__ == "__main__":
    main()
