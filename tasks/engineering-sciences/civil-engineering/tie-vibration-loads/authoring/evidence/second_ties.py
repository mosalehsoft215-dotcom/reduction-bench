"""Independent evaluation of one roof tie (does not import the reference, no beam elements, no tables).

Differences from the reference:
  - spectra by Welch with the default window (fifth-length segments), peaks of 12 dB prominence located by a parabola
    through the logarithm of the five lines around the top;
  - bars: Rayleigh-Ritz with 24 trial functions, each a combination of two shifted Legendre polynomials that
    vanishes at both ends; the end restraints add c_lo phi'(0) phi'(0)^T + c_hi phi'(1) phi'(1)^T (c = k_end Lv / EI,
    up to 1e4 for a held end); operating shapes at the three pickups from the H1 transfer of each pickup to the
    liveliest one at the top line of each peak (real part, unit vector);
  - backstays: sine series (exact for a string with small bending), the sag-extensibility term as the rank-one matrix
    of the integrals of the sines;
  - the unknowns come from a scan (bars: both loads for each restraint pair of a 4 x 4 grid) of a capped log-distance
    misfit plus, for bars, the shape misfit of the matched modes, then Nelder-Mead on the same misfit in (T_pre,
    T_post, log c_lo, log c_hi) from the two best scan points.  Backstays: both ratios T / (kappa m) by Nelder-Mead at a fixed mass factor
    kappa; kappa from the shim change (unstressed strand length between the anchorages of a flat-sagged elastic
    tie at 20 C, (L + q^2 L^3 / (24 T^2)) / ((1 + T / EA)(1 + 12e-6 (t - 20)))) by bisection; three rounds of the two.

    python second_ties.py <tie folder> <output file>
"""
import csv
import json
import os
import sys
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")            # small matrices: one BLAS thread is the fast setting on two CPUs

import numpy as np  # noqa: E402
from numpy.polynomial import legendre
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.signal import csd, find_peaks, welch

G = 9.80665
NB = 24
GX, GW = legendre.leggauss(120)
GX, GW = 0.5 * (GX + 1), 0.5 * GW


def spectrum_peaks(acc, fs):
    """peaks of the spectrum (several pickups: the sum of their auto-spectra)"""
    f, p = welch((acc - np.mean(acc, axis=0)).T, fs=fs, nperseg=len(acc) // 5)
    p = p.sum(axis=0) if p.ndim > 1 else p
    lg = np.log(p + 1e-300)
    idx, _ = find_peaks(10 * np.log10(p + 1e-30), prominence=12.0, distance=3)
    out = []
    for i in idx:
        if f[i] > 0.12 and 2 <= i < len(f) - 2:
            c = np.polyfit(np.arange(-2, 3), lg[i - 2:i + 3], 2)
            shift = -c[1] / (2 * c[0]) if c[0] < 0 else 0.0
            out.append(f[i] + np.clip(shift, -1, 1) * (f[1] - f[0]))
    return np.array(out)


def ritz_basis():
    """trial functions on 0..1 built from Legendre polynomials P_j(2x - 1): P_j - P_(j+2) vanishes at both ends"""
    t = 2 * GX - 1
    coefs = []
    for j in range(NB):
        c = np.zeros(NB + 3)
        c[j], c[j + 2] = 1.0, -1.0
        coefs.append(c / np.sqrt(np.sum(GW * legendre.legval(t, c) ** 2)))
    v0 = np.array([legendre.legval(t, c) for c in coefs])
    v1 = np.array([2 * legendre.legval(t, legendre.legder(c, 1)) for c in coefs])
    v2 = np.array([4 * legendre.legval(t, legendre.legder(c, 2)) for c in coefs])
    s0 = np.array([2 * legendre.legval(-1.0, legendre.legder(c, 1)) for c in coefs])
    s1 = np.array([2 * legendre.legval(1.0, legendre.legder(c, 1)) for c in coefs])
    return np.array(coefs).T, (v2 * GW) @ v2.T, (v1 * GW) @ v1.T, (v0 * GW) @ v0.T, np.outer(s0, s0), np.outer(s1, s1)


BASIS = ritz_basis()


def bar_modes(xi2, c_lo, c_hi, rows, kmax=10):
    """dimensionless frequencies (Hz on the unit bar with EI = m = 1) and unit shape vectors at the stations, rows being
    the trial functions evaluated there"""
    coefs, Kb, Kg, M, E0, E1 = BASIS
    w, v = eigh(Kb + xi2 * Kg + c_lo * E0 + c_hi * E1, M, subset_by_index=[0, kmax - 1])
    at = rows @ v
    return np.sqrt(w) / (2 * np.pi), at / np.linalg.norm(at, axis=0)


def stay_freqs(T, L, EI, EA, m, theta, kappa, n=40, kmax=26):
    mk = m * kappa
    j = np.arange(1, n + 1)
    q = mk * G * np.cos(np.radians(theta))
    sag = q * L * L / (8 * T)
    Le = L * (1 + 8 * (sag / L) ** 2)
    c = (1 - np.cos(j * np.pi)) / (j * np.pi) * L           # integral of sin(j pi x / L) over the length
    K = np.diag(EI * (j * np.pi / L) ** 4 * L / 2 + T * (j * np.pi / L) ** 2 * L / 2) + EA / Le * (q / T) ** 2 * np.outer(c, c)
    w = eigh(K, eigvals_only=True) / (mk * L / 2)
    return np.sqrt(w[:kmax]) / (2 * np.pi)


def misfit(obs, model, cap):
    if not np.all(np.isfinite(model)):
        return 1e9
    d = np.abs(np.log(obs[:, None] / model[None, :])).min(axis=1)
    band = model[(model > obs.min() * 0.97) & (model < obs.max() * 1.03)]
    lonely = np.abs(np.log(band[:, None] / obs[None, :])).min(axis=1) > cap if len(band) else np.zeros(0)
    return np.sum(np.minimum(d, cap) ** 2) + 0.35 * cap * cap * np.sum(lonely)


def misfit_rows(obs, table, cap):
    """misfit of every row of table (one ladder per row) in one pass"""
    lt = np.log(table)
    d = np.abs(np.log(obs)[None, :, None] - lt[:, None, :]).min(axis=2)
    inside = (table > obs.min() * 0.97) & (table < obs.max() * 1.03)
    lonely = np.abs(lt[:, :, None] - np.log(obs)[None, None, :]).min(axis=2) > cap
    return np.sum(np.minimum(d, cap) ** 2, axis=1) + 0.35 * cap * cap * np.sum(inside & lonely, axis=1)


def load(folder):
    sheet = json.loads((folder / "tie_sheet.json").read_text())
    recs = {}
    for tag, name in sheet["files"].items():
        with open(folder / name) as fh:
            rows = list(csv.DictReader(fh))
        recs[tag] = np.array([[float(r[k]) for k in rows[0]] for r in rows])
        if recs[tag].shape[1] == 1:
            recs[tag] = recs[tag][:, 0]
    return sheet, recs


def shapes_at(acc, fs, freqs):
    """H1 transfer of each pickup to the liveliest one at the top line of each peak, real part, unit vector"""
    x = acc - acc.mean(axis=0)
    nper = len(x) // 5
    f, auto = welch(x.T, fs=fs, nperseg=nper)
    out = []
    for fk in freqs:
        i = int(np.argmin(np.abs(f - fk)))
        ref = int(np.argmax(auto[:, i]))
        h = np.array([csd(x[:, ref], x[:, j], fs=fs, nperseg=nper)[1][i] for j in range(3)]) / auto[ref, i]
        u = np.real(h)
        out.append(u / np.linalg.norm(u))
    return np.array(out).reshape(-1, 3)


def bar_misfit(obs, shp, model, at, cap):
    d = np.abs(np.log(obs[:, None] / model[None, :]))
    k = d.argmin(axis=1)
    dk = d[np.arange(len(obs)), k]
    J = np.sum(np.minimum(dk, cap) ** 2) / 0.0015 ** 2
    for i in np.flatnonzero(dk < cap):
        v = at[:, k[i]]
        J += min(np.sum((shp[i] - np.sign(v @ shp[i]) * v) ** 2), 0.05) / 0.01 ** 2
    return J


def bar(sheet, recs, fs):
    m, Lv = sheet["unit_mass_kg_per_m"], sheet["clear_length_m"]
    d = sheet["make_up"]["bar_diameter_mm"] / 1000
    EI = sheet["young_modulus_mpa"] * 1e6 * np.pi * d ** 4 / 64
    h = sheet["pickups"]["metres_above_lower_block_face"]
    st = [h["a_lo_mg"] / Lv, h["a_mid_mg"] / Lv, h["a_hi_mg"] / Lv]
    obs = []
    for tag in ("pre", "post"):
        pk = spectrum_peaks(recs[tag], fs)
        obs.append((pk, shapes_at(recs[tag], fs, pk)))
    scale = np.sqrt(EI / (m * Lv ** 4))
    xis = np.geomspace(3.0, 120.0, 60)
    ranked = []
    st = legendre.legval(2 * np.asarray(st) - 1, BASIS[0]).T
    for c1 in (0.5, 3.0, 20.0, 1e4):
        for c2 in (0.5, 3.0, 20.0, 1e4):
            table = [bar_modes(x * x, c1, c2, st) for x in xis]
            tot, Ts = 0.0, []
            for pk, shp in obs:
                js = [bar_misfit(pk, shp, fr * scale, at, 0.02) for fr, at in table]
                j = int(np.argmin(js))
                tot += js[j]
                Ts.append((xis[j] / Lv) ** 2 * EI)
            ranked.append((tot, c1, c2, Ts))
    ranked.sort(key=lambda r: r[0])
    best = None
    for tot, c1, c2, Ts in ranked[:2]:
        def J(z, cap):
            if not (z[0] > 0 and z[1] > 0 and -3.0 < z[2] < 9.3 and -3.0 < z[3] < 9.3):
                return 1e12
            out = 0.0
            for T, (pk, shp) in zip(z[:2], obs):
                fr, at = bar_modes(Lv * Lv * T / EI, np.exp(z[2]), np.exp(z[3]), st)
                out += bar_misfit(pk, shp, fr * scale, at, cap)
            return out
        z = np.array([Ts[0], Ts[1], np.log(c1), np.log(c2)])
        for cap in (0.02, 0.008):
            z = minimize(J, z, args=(cap,), method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-9, "maxiter": 450}).x
        score = J(z, 0.02)
        if best is None or score < best[0]:
            best = (score, z)
    return best[1][0], best[1][1]


def stay(sheet, pb, pa):
    m, L = sheet["unit_mass_kg_per_m"], sheet["working_point_length_m"]
    E = sheet["young_modulus_mpa"] * 1e6
    EI = sheet["make_up"]["strand_count"] * E * np.pi * 0.0157 ** 4 / 64
    EA = sheet["steel_area_mm2"] * 1e-6 * E
    th = sheet["chord_slope_deg"]
    first = []
    for obs in (pb, pa):
        T0 = 4 * m * L * L * np.median(np.diff(np.sort(obs))) ** 2
        for lo, hi, num in ((0.2, 5.0, 240), (0.97, 1.03, 61)):
            grid = np.geomspace(lo * T0, hi * T0, num)
            T0 = grid[int(np.argmin([misfit(obs, stay_freqs(T, L, EI, EA, m, th, 1.0), 0.02) for T in grid]))]
        first.append(T0)
    shim = sheet["retensioning"]["shim_change_mm"] / 1e3
    gc = G * np.cos(np.radians(th))

    temp = sheet["steel_temperature_c"]

    def strand_length(T, mk, t):                            # unstressed, referred to 20 C (steel 12e-6 per K)
        return (L + (mk * gc) ** 2 * L ** 3 / (24 * T * T)) / ((1 + T / EA) * (1 + 12e-6 * (t - 20.0)))

    kappa, z = 1.0, np.array(first) / m                     # z = T / (kappa m) per record
    for _ in range(3):
        def J(y, cap, k=kappa):
            return misfit(pb, stay_freqs(y[0] * k * m, L, EI, EA, m, th, k), cap) + misfit(pa, stay_freqs(y[1] * k * m, L, EI, EA, m, th, k), cap)
        for cap in (0.02, 0.008, 0.004):
            z = minimize(J, z, args=(cap,), method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-13, "maxiter": 3000}).x
        lo, hi = 0.6, 2.0
        f = lambda k: strand_length(k * m * z[0], k * m, temp["pre"]) - strand_length(k * m * z[1], k * m, temp["post"]) - shim
        if f(lo) * f(hi) > 0:
            break
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
        kappa = 0.5 * (lo + hi)
    return z[0] * kappa * m, z[1] * kappa * m


def solve(folder):
    sheet, recs = load(folder)
    fs = sheet["samples_per_second"]
    if "bar_diameter_mm" in sheet["make_up"]:
        tb, ta = bar(sheet, recs, fs)
    else:
        tb, ta = stay(sheet, spectrum_peaks(recs["pre"], fs), spectrum_peaks(recs["post"], fs))
    return tb / 1e3, ta / 1e3


def main(folder, out):
    tb, ta = solve(Path(folder))
    Path(out).write_text(json.dumps({"load_pre_kn": round(tb, 3), "load_post_kn": round(ta, 3)}) + chr(10))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
