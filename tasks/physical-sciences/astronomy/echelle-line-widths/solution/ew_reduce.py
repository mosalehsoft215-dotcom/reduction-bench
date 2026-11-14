"""Reduces one exposure of echelle order 74 (target star plus its telluric standard) to the star's radial velocity and
the equivalent width of the K I 7665 A feature.

    python3 ew_reduce.py <exposure directory> <result.json>

Forward model of the notes on a 0.004 A grid (five samples per 0.02 A pixel): counts = cubic continuum x
[stellar x O2 x interstellar transmission] convolved with a Gaussian instrumental profile of FWHM 7650 A / R, binned
to the pixels, plus a scattered-light constant; the pipeline scale is the true wavelength plus an offset; each O2
line is a Lorentzian of HWHM 0.018 A cut off 4 A from its centre.
1. Standard (a short exposure, airmass constant at the header value): its offset, scattered light, cubic continuum,
   resolving power R, the airmass-1 peak optical depth of every listed O2 line for the standard's air and the
   interstellar K I doublet (peak depth, velocity started from the deepest pixel near 7698.96 A clear of O2 lines,
   Doppler parameter) are fitted by least squares; the line depths enter through an analytic Jacobian. R comes from
   the shapes of the O2 lines, whose own width is fixed, so it need not be the 60 000 of the evenly filled slit.
2. Target: the O2 transmission is the time mean of exp(-X(t) tau) with X moving linearly from the header's value at
   shutter opening to a fitted value at shutter closing; tau are the standard's depths times the ratio of the lower
   levels' Boltzmann populations at the two air temperatures, (1 + 250 u) exp(-c2 E'' u) with u = 1/T_target -
   1/T_standard fitted (the depths' trend with E'' gives it); fitted with them are the target's offset, scattered
   light, cubic continuum, resolving power, radial velocity, the common Doppler width and the central optical depth
   of each listed stellar line. The closing airmass comes from how the O2 lines of the target, weak and saturated,
   differ from those of the standard.
3. Moonlight and gas, each kept when it lowers chi-square by more than 400: scattered moonlight (the solar spectrum
   of the atlas at a small fitted velocity under the same O2 band, and over 2 % of the light); an interstellar K I
   doublet in the target's own light (started at the lowest three-pixel mean of data / model - 1 along the
   7698.96 A member's track, 0.15 A clear of O2 lines and 0.3 A of the star's lines, and deeper than 0.03).  Then a
   fringe sinusoid from a periodogram of data / model - 1 over periods of 7-18 A, kept when its amplitude exceeds
   1.5 % (no chi-square test).
4. Equivalent width: the fitted stellar transmission alone (no interstellar gas), convolved with the target's
   instrumental profile, integrated as 1 - F over +-0.9 A around 7664.91 (1 + v/c).
"""
import json
import math
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy import sparse  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

C_KMS = 299792.458
PIX, FINE = 0.02, 0.004
OVER = 5
GAMMA = 0.018
R_SLIT = 60000.0
SIG_PER_R = 7650.0 / 2.3548 / FINE          # Gaussian sigma in fine samples = SIG_PER_R / R
T0 = time.monotonic()


def read_order(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1]


def blur(y, power):
    return gaussian_filter1d(y, SIG_PER_R / power, mode="nearest", truncate=6.0)


def bin_pix(y):
    return y.reshape(-1, OVER).mean(1)


def fine_grid(lam_obs, offset):
    return ((lam_obs - offset)[:, None] + (np.arange(OVER) - 2) * FINE).ravel()


def cubic(lam_f, c):
    x = (lam_f - 7650.0) / 60.0
    return c[0] * (1 + x * (c[1] + x * (c[2] + x * c[3])))


def o2_basis(lam_f, lines):
    cols, rows, vals = [], [], []
    for j, c in enumerate(lines):
        i0, i1 = np.searchsorted(lam_f, [c - 4.0, c + 4.0])
        rows.append(np.arange(i0, i1))
        cols.append(np.full(i1 - i0, j))
        vals.append(1 / (1 + ((lam_f[i0:i1] - c) / GAMMA) ** 2))
    return sparse.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(lam_f.size, len(lines)))


def mean_transmission(tau, x_open, x_close):
    """Time mean of exp(-X tau) for X moving linearly from x_open to x_close over the exposure."""
    step = (x_close - x_open) * tau
    safe = np.where(np.abs(step) < 1e-9, 1.0, step)
    return np.exp(-x_open * tau) * np.where(np.abs(step) < 1e-9, 1.0 - 0.5 * step, -np.expm1(-step) / safe)


def envelope_cubic(lam, cnt):
    keep = np.ones(lam.size, bool)
    for _ in range(8):
        p = np.polyfit(lam[keep] - 7650.0, cnt[keep], 3)
        keep = cnt > np.polyval(p, lam - 7650.0) * 0.995 if keep.sum() > 400 else keep
    p = np.polyfit(lam[keep] - 7650.0, cnt[keep], 3)
    return [p[3], p[2] * 60.0 / p[3], p[1] * 3600.0 / p[3], p[0] * 216000.0 / p[3]]


def line_windows(lam_f, lines):
    """Per O2 line: a window of fine samples (multiple of OVER at both ends) round its core and its profile there."""
    out = []
    for c in lines:
        mid = int(np.searchsorted(lam_f, c))
        a = max(0, (mid - 700) // OVER * OVER)
        b = min(lam_f.size, (mid + 700) // OVER * OVER)
        out.append((a, b, 1 / (1 + ((lam_f[a:b] - c) / GAMMA) ** 2)))
    return out


K_DOUBLET = ((7664.91, 1.0), (7698.96, 0.499))      # interstellar K I: rest air wavelength, share of the peak depth


def tau_ism(lf, t0, v, b):
    out = np.zeros_like(lf)
    for c, share in K_DOUBLET:
        centre = c * (1 + v / C_KMS)
        width = centre * abs(b) / C_KMS
        i0, i1 = np.searchsorted(lf, [centre - 6 * width, centre + 6 * width])
        out[i0:i1] += share * t0 * np.exp(-((lf[i0:i1] - centre) / width) ** 2)
    return out


def fit_standard(lam, cnt, airmass, lines, free_power=True, ism=True):
    n_l = len(lines)
    c0 = envelope_cubic(lam, cnt)
    norm = cnt / bin_pix(cubic(fine_grid(lam, 0.0), c0))
    d0 = [-np.log(np.clip(norm[np.argmin(np.abs(lam - c))], 0.02, 1.0)) / airmass * 2.0 + 0.05 for c in lines]
    w = 1 / np.sqrt(np.clip(cnt, 25, None))
    cache = {}

    def basis_at(off):
        if off not in cache:
            if len(cache) > 3:
                cache.clear()
            lf = fine_grid(lam, off)
            cache[off] = (lf, o2_basis(lf, lines), line_windows(lf, lines))
        return cache[off]

    def unpack(p):
        return p[0], p[1], p[2:6], p[6], p[7:7 + n_l], (p[7 + n_l:] if ism else (0.0, 0.0, 2.0))

    def resid(p):
        off, floor, c, power, d, g_ism = unpack(p)
        lf, basis, _ = basis_at(off)
        light = np.exp(-airmass * (basis @ d) - tau_ism(lf, *g_ism)) * cubic(lf, c)
        return (bin_pix(blur(light, power)) + floor - cnt) * w

    def jac(p):
        off, floor, c, power, d, g_ism = unpack(p)
        lf, basis, wins = basis_at(off)
        g = -airmass * np.exp(-airmass * (basis @ d) - tau_ism(lf, *g_ism)) * cubic(lf, c)
        j = np.zeros((lam.size, p.size))
        for k, (a, b, prof) in enumerate(wins):
            j[a // OVER:b // OVER, 7 + k] = bin_pix(blur(prof * g[a:b], power)) * w[a // OVER:b // OVER]
        j[:, 1] = w
        base = resid(p)
        for k in (0, 2, 3, 4, 5, 6) + ((7 + n_l, 8 + n_l, 9 + n_l) if ism else ()):
            q = p.copy()
            h = 1e-5 if k == 0 else 1e-6 * max(abs(p[k]), 1.0)
            q[k] += h
            j[:, k] = (resid(q) - base) / h
        return j

    # interstellar K I start: velocity from the deepest residual near 7698.96 A, where no O2 line sits
    o2_near = [l for l in lines if abs(l - 7698.96) < 2]
    near = np.abs(lam - 7698.96) < 1.2
    if o2_near:
        near &= np.all(np.abs(lam[:, None] - np.array(o2_near)[None, :]) > 0.12, axis=1)
    v0 = (lam[near][np.argmin(norm[near])] / 7698.96 - 1) * C_KMS if near.any() else 0.0
    p0 = np.r_[0.0, max(np.percentile(cnt, 0.5), 0.0), c0, R_SLIT, d0, 0.3, v0, 2.0]
    lo = np.r_[-0.3, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 3.0e4, np.zeros(n_l), 0.0, -80.0, 0.5]
    hi = np.r_[0.3, np.inf, np.inf, np.inf, np.inf, np.inf, 1.6e5, np.full(n_l, np.inf), 5.0, 80.0, 8.0]
    if not free_power:
        lo[6], hi[6] = R_SLIT - 1e-6, R_SLIT + 1e-6
    if not ism:                                   # the standard taken as featureless: no interstellar parameters
        p0, lo, hi = p0[:7 + n_l], lo[:7 + n_l], hi[:7 + n_l]
    sol = least_squares(resid, p0, jac=jac, bounds=(lo, hi), x_scale="jac", max_nfev=200)
    return sol.x[:7 + n_l]


def stellar_tau(lf, rv, b, depths, rest):
    s = 1 + rv / C_KMS
    tau = np.zeros_like(lf)
    for c, t in zip(rest, depths):
        i0, i1 = np.searchsorted(lf, [c * s - 6 * b, c * s + 6 * b])
        tau[i0:i1] += t * np.exp(-((lf[i0:i1] - c * s) / b) ** 2)
    return tau


SUN_DEPTHS = [1.95, 0.46, 0.42, 0.31, 0.55, 0.27, 0.48, 0.36, 0.22]      # replaced by lines.json (solar atlas)
SUN_DOP = 0.072
C2, T_AIR = 1.4388, 250.0                    # second radiation constant (cm K); the air's temperature, about 250 K


def warmth(energies, uu):
    """Ratio of the target's airmass-1 O2 depths to the standard's for uu = 1e4 (1/T_target - 1/T_standard): the
    lower levels' Boltzmann populations, (T_s / T_t) exp(-c2 E'' (1/T_t - 1/T_s)), with T_s / T_t = 1 + T_s u."""
    u = uu * 1e-4
    return (1 + T_AIR * u) * np.exp(-C2 * energies * u)


def ism_track(lam, ratio, off, lines, rest, rv):
    """Velocity of the deepest residual along the 7698.96 A member's track, clear of O2 and stellar lines."""
    o2a, st = np.asarray(lines), np.asarray(rest) * (1 + rv / C_KMS)
    best = (0.0, None)
    for v in np.arange(-60.0, 60.01, 0.5):
        at = 7698.96 * (1 + v / C_KMS)
        if np.min(np.abs(o2a - at)) < 0.15 or np.min(np.abs(st - at)) < 0.3:
            continue
        i = int(np.argmin(np.abs(lam - off - at)))
        depth = float(np.mean(ratio[max(i - 1, 0):i + 2]))
        if depth < best[0]:
            best = (depth, v)
    return best[1]


def fit_target(lam, cnt, x_open, lines, o2_depths, rest, free_power=True, free_close=True, ramp=True,
               power=R_SLIT, extras=True, energies=None, temp=True, ism_t=True, clouds=3):
    c0 = envelope_cubic(lam, cnt)
    floor0 = max(np.min(cnt[(lam > 7595) & (lam < 7620)]), 0.0)
    w = 1 / np.sqrt(np.clip(cnt, 25, None))
    cache = {}
    n_base = 10 + len(rest)
    i_u, i_ism = n_base + 4, n_base + 5
    state = {"period": 10.0}
    energies = np.zeros(len(lines)) if energies is None else np.asarray(energies, float)

    def o2_at(off):
        if off not in cache:
            if len(cache) > 3:
                cache.clear()
            lf = fine_grid(lam, off)
            cache[off] = (lf, o2_basis(lf, lines))
        return cache[off]

    def model(p):
        off, floor, c, power, x_close, rv, b, t = p[0], p[1], p[2:6], p[6], p[7], p[8], p[9], p[10:n_base]
        moon, v_moon, fs, fc = p[n_base:n_base + 4]
        lf, basis = o2_at(off)
        tau = basis @ (o2_depths * warmth(energies, p[i_u]))
        own = stellar_tau(lf, rv, b, t, rest)
        for m in range(3):                          # up to three interstellar clouds in the target's light
            if p[i_ism + 3 * m] > 0:
                own = own + tau_ism(lf, p[i_ism + 3 * m], p[i_ism + 3 * m + 1], p[i_ism + 3 * m + 2])
        light = np.exp(-own) * cubic(lf, c)
        light = light + moon * np.exp(-stellar_tau(lf, v_moon, SUN_DOP, SUN_DEPTHS, rest))
        out = bin_pix(blur(mean_transmission(tau, x_open if ramp else x_close, x_close) * light, power)) + floor
        ph = 2 * np.pi * lam / state["period"]
        return out * (1 + fs * np.sin(ph) + fc * np.cos(ph))

    def resid(p):
        return (model(p) - cnt) * w

    start = np.r_[0.0, floor0, c0, power, x_open, 0.0, 0.09, np.full(len(rest), 0.5), 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 2.5, 0.0, 0.0, 2.5]
    best_rv, best_cost = 0.0, np.inf
    for rv in np.arange(-90.0, 90.1, 3.0):
        start[8] = rv
        cost = np.sum(resid(start) ** 2)
        if cost < best_cost:
            best_rv, best_cost = rv, cost
    start[8] = best_rv
    lo = np.r_[-0.3, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 3.0e4, 0.5, -150.0, 0.03, np.zeros(len(rest)),
               0.0, -3.0, -1e-9, -1e-9, -8.0, 0.0, -80.0, 0.5, 0.0, -80.0, 0.5, 0.0, -80.0, 0.5]
    hi = np.r_[0.3, np.inf, np.inf, np.inf, np.inf, np.inf, 1.6e5, 6.0, 150.0, 0.3, np.full(len(rest), 20.0),
               1e-9, 3.0, 1e-9, 1e-9, 8.0, 1e-9, 80.0, 8.0, 1e-9, 80.0, 8.0, 1e-9, 80.0, 8.0]
    if not temp:
        lo[i_u], hi[i_u] = -1e-9, 1e-9
    if not free_power:
        lo[6], hi[6] = power - 1e-6, power + 1e-6
    if not free_close:
        lo[7], hi[7] = x_open - 1e-9, x_open + 1e-9
    nfev = 300 if time.monotonic() - T0 < 30 else 60

    def sub_fit(active, x0, n_ev):
        # least squares over the active parameters only; the others stay at their values in x0
        idx = np.flatnonzero(active)
        full = x0.copy()

        def r_sub(q):
            full[idx] = q
            return resid(full)
        got = least_squares(r_sub, x0[idx], bounds=(lo[idx], hi[idx]), x_scale="jac", max_nfev=n_ev)
        out = x0.copy()
        out[idx] = got.x
        got.x = out
        return got

    base = np.zeros(start.size, bool)
    base[:n_base] = True
    base[i_u] = temp
    sol = sub_fit(base, start, nfev)
    if not extras:
        return sol.x
    # scattered moonlight: the solar spectrum (the notes' atlas) at a small velocity, added to the target's light;
    # kept only when it lowers chi-square by far more than the noise allows and carries at least 2 % of the light
    lo_m, hi_m = lo.copy(), hi.copy()
    hi_m[n_base] = 0.5 * abs(c0[0])
    trial = sol.x.copy()
    trial[n_base] = 0.05 * abs(c0[0])
    lo_keep, hi_keep = lo.copy(), hi.copy()
    lo[:], hi[:] = lo_m, hi_m
    with_moon = base.copy()
    with_moon[n_base:n_base + 2] = True
    moon_fit = sub_fit(with_moon, trial, 25)
    lo[:], hi[:] = lo_keep, hi_keep
    if 2 * (sol.cost - moon_fit.cost) > 400.0 and moon_fit.x[n_base] > 0.02 * abs(c0[0]):
        sol, lo, hi = moon_fit, lo_m, hi_m
        base = with_moon
    # interstellar K I in the target's own light, cloud by cloud: the deepest residual along the 7698.96 A member's
    # track starts the first Gaussian doublet (7664.91 A at twice the 7698.96 A depth), each further cloud is tried 3
    # km/s to either side of the last one found, the side where the 7698.96 A residual is deeper with the longer fit
    # (those found restarted at most 2 km/s wide); all clouds found so far are
    # refitted with it, and it is kept when it lowers chi-square by far more than the noise allows and is at least
    # 0.03 deep.
    # A sightline that crosses clouds a few km/s apart, narrower than the instrumental profile and saturated, needs
    # them all: one broad doublet in their place misjudges the 7664.91 A line inside the K I window.
    for m in range(clouds if ism_t else 0):
        k = i_ism + 3 * m
        if m == 0:
            v_i = ism_track(lam, cnt / model(sol.x) - 1, sol.x[0], lines, rest, sol.x[8])
            starts = [] if v_i is None else [v_i]
        else:                                       # a further cloud next to those found, on either side
            found = [sol.x[i_ism + 3 * j + 1] for j in range(m)]
            # the side where the 7698.96 A member's residual is deeper, 3 km/s from the last cloud found
            ratio = cnt / model(sol.x) - 1

            def dip(v):
                i = int(np.argmin(np.abs(lam - sol.x[0] - 7698.96 * (1 + v / C_KMS))))
                return float(np.mean(ratio[max(i - 1, 0):i + 2]))
            starts = sorted((found[-1] - 3.0, found[-1] + 3.0), key=dip)
        best = None
        for v_i in starts:
            trial = sol.x.copy()
            trial[k:k + 3] = 0.3, v_i, 1.5
            for j in range(m):                      # the clouds already found start narrower
                trial[i_ism + 3 * j + 2] = min(trial[i_ism + 3 * j + 2], 2.0)
            lo_i, hi_i = lo.copy(), hi.copy()
            hi_i[k] = 5.0
            lo_keep, hi_keep = lo.copy(), hi.copy()
            lo[:], hi[:] = lo_i, hi_i
            with_ism = base.copy()
            with_ism[k:k + 3] = True
            got = sub_fit(with_ism, trial, 25 if v_i == starts[0] else 15)
            lo[:], hi[:] = lo_keep, hi_keep
            if best is None or got.cost < best[0].cost:
                best = (got, lo_i, hi_i, with_ism)
        if best is not None and 2 * (sol.cost - best[0].cost) > 400.0 and best[0].x[k] > 0.03:
            sol, lo, hi, base = best
        else:
            break
    # residual fringing: a sinusoid in wavelength (period 7-18 A) multiplying the counts, found by a periodogram of
    # the ratio data / model and kept when its amplitude exceeds 1.5 %
    ratio = cnt / model(sol.x) - 1
    best = (0.0, 10.0, 0.0, 0.0)
    for period in np.arange(7.0, 18.0, 0.02):
        ph = 2 * np.pi * lam / period
        A = np.vstack([np.sin(ph), np.cos(ph), np.ones_like(lam)]).T
        coef = np.linalg.lstsq(A, ratio, rcond=None)[0]
        amp = math.hypot(coef[0], coef[1])
        if amp > best[0]:
            best = (amp, period, coef[0], coef[1])
    amp, period, fs, fc = best
    if amp > 0.015:
        state["period"] = period
        start2 = sol.x.copy()
        start2[n_base + 2], start2[n_base + 3] = fs, fc
        lo[n_base + 2:n_base + 4], hi[n_base + 2:n_base + 4] = -0.2, 0.2
        active = base.copy()
        active[n_base + 2:n_base + 4] = True
        sol = sub_fit(active, start2, 25)
    return sol.x


def reduce_exposure(folder, free_power=True, free_close=True, full=False, ramp=True, shared_power=False, ism=True,
                    extras=True, temp=True, ism_t=True, clouds=3):
    global T0
    T0 = time.monotonic()                     # the time guard counts from the start of this exposure
    folder = Path(folder)
    head = json.loads((folder / "exposure.json").read_text())
    lines_doc = json.loads((folder / "lines.json").read_text())
    o2 = lines_doc["telluric_o2_rest_air_a"]
    rest = lines_doc["stellar_rest_air_a"]
    lam_s, cnt_s = read_order(folder / "standard_order.csv")
    lam_t, cnt_t = read_order(folder / "target_order.csv")
    global SUN_DEPTHS, SUN_DOP
    SUN_DEPTHS = lines_doc.get("solar_peak_optical_depths", SUN_DEPTHS)
    SUN_DOP = lines_doc.get("solar_doppler_width_a", SUN_DOP)
    ps = fit_standard(lam_s, cnt_s, head["standard_airmass"], o2, free_power, ism)
    p = fit_target(lam_t, cnt_t, head["target_airmass"], o2, ps[7:], rest, free_power and not shared_power,
                   free_close, ramp, ps[6] if shared_power else R_SLIT, extras,
                   lines_doc.get("telluric_o2_lower_state_energy_cm1"), temp, ism_t and extras, clouds)
    if full:
        return ps, p
    off, power, rv, b, t = p[0], p[6], p[8], p[9], p[10:10 + len(rest)]
    lf = fine_grid(lam_t, off)
    clean = blur(np.exp(-stellar_tau(lf, rv, b, t, rest)), power)
    centre = rest[0] * (1 + rv / C_KMS)
    window = np.abs(lf - centre) <= 0.9
    return {"radial_velocity_kms": float(rv), "ew_7665_ma": float((1 - clean[window]).sum() * FINE * 1e3)}


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(reduce_exposure(sys.argv[1])))
