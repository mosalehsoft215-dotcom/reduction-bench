"""Reference analysis of the ring occultation.

What reaches a detector at the bottom of a dip has four sources, and the optical depth of the ring is
only as good as the three that are not the direct beam.

1. Scale. stations.csv gives the tracks as angles. They become km with the distance of the body at the
   reference epoch, interpolated between the two tabulated ranges (the range grows by 21 km/s).
2. Light that the ring never covers, in units of the primary star (unocculted()):
   - the companion, 0.31 arcsec away: its track passes 1500 km from the ring's centre;
   - the centaur. Its mean magnitude follows from the absolute magnitude, the two distances and the
     phase angle. Its light curve has two maxima per rotation and the tabulated maximum is in the time
     of the body, so the event is placed on it after taking off the light time (55 min, a quarter of
     a brightness cycle): the centaur was near minimum;
   - the ring in reflected sunlight, scaled from the opening angle of 2025 to the one fitted here.
   The catalogue magnitude of the star belongs to the pair; the ratios wanted are to the primary alone.
3. Each light curve is divided by a 6 s running median and gets one noise level from its scatter. The two
   deep dips of each station give sixteen points in the sky; a conic through them starts the ring's
   centre, radius, opening angle and position angle.
4. The direct beam (ring_model.star_flux) is fitted to all stations at once in windows of 5 s around the
   dips: centre, radius, opening angle, position angle, mean width, width modulation, longitude of
   periapse, and the mean and modulation of the line-of-sight optical depth of the direct beam. The
   transparency of the sky wanders by a few percent within seconds, so each window has its own cubic
   baseline, which multiplies everything in the aperture.
5. The particles are 6 to 30 m in radius: what one of them diffracts leaves it within about
   lambda / (2 r), and at the distance of the body that is 13 to 66 km on the ground, a few times the
   width of the ring and not more. A station inside the shadow of the ring therefore stands in the
   diffracted light of the neighbouring particles: 3 to 6 percent of the primary at the bottom of a
   dip whose direct beam is 5 to 7 percent. ring_model.halo computes it from the ring being fitted:
   the stated sizes, the distance, and the particle count on each line of sight, all orders of repeated
   diffraction. Fit and halo are iterated until they agree.
6. From the direct beam to the particles. The ring is many particles thick, so the direct beam falls as
   exp(-tau_los) with tau_los = 2 tau / sin B for the geometric normal optical depth tau: each particle
   takes its shadow and as much again by diffraction out of the direct beam.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import median_filter
from scipy.optimize import least_squares

from ring_model import AU_KM, halo, near_edges, patterns, ring_plane, star_flux

# ---- numbers of event_notes.txt ---------------------------------------------------------------
RANGE_AU = (6.598214, 6.610574)               # geocentric distance at 0 h UT, 13 and 14 September
EPOCH_DAY = 13 + (13 * 3600 + 10 * 60) / 86400.0      # reference epoch, day of September (UT)
R_SUN_AU, PHASE_DEG = 7.5482, 1.912
H_BODY, H_RING, RING_ANGLE_DEG, PHASE_COEFF = 6.330, 8.214, 21.5, 0.061
ROT_PERIOD_H, ROT_AMPLITUDE, ROT_MAXIMUM_DAY = 7.0040, 0.200, 11.39900
I_PAIR, DELTA_I = 13.752, 2.50
LIGHT_S_PER_AU = 499.004784
R_MIN_KM, R_MAX_KM, SLOPE = 0.006, 0.030, 3.0  # particle radii, dn/dr ~ r^-3
MAS = np.pi / 180 / 3600 / 1000
# ---- choices of the analysis ------------------------------------------------------------------
HALF = 2.5                                    # half-width of a fit window, s
NEAR_KM = 4.0                                 # diffraction integral within this distance of an edge, geometric optics beyond
ORDERS = 8                                    # diffractions kept (the ninth carries less than 1e-4 of the diffracted light)
DEGREE = 3                                    # baseline polynomial of a window
_g = np.polynomial.legendre.leggauss(16)
QUICK = ((_g[0] + 1) / 2, _g[1] / 2)
KEYS = ("semi_major_axis_km", "opening_angle_deg", "position_angle_deg", "width_periapse_km", "width_apoapse_km",
        "tau_periapse", "tau_apoapse")
#                  xi_c  eta_c    a    B     P     W0     q   theta_p  ln tau_los   h
SCALE = np.array([0.5, 0.5, 0.5, 0.2, 0.2, 0.3, 0.05, 0.1, 0.05, 0.05])
LOWER = np.array([-300, -300, 200, 5, -90, 0.5, -0.95, -20, -4.0, -0.95])
UPPER = np.array([300, 300, 700, 85, 270, 60, 0.95, 20, 2.5, 0.95])


def distance_km():
    """distance of the body at the reference epoch"""
    return (RANGE_AU[0] + (EPOCH_DAY - 13.0) * (RANGE_AU[1] - RANGE_AU[0])) * AU_KM


def unocculted(b_deg):
    """centaur, its ring in reflected sunlight and the companion, in units of the primary star.
    b_deg: opening angle of the ring at the event (the fitted one)"""
    dist = distance_km() / AU_KM
    geo = 5 * np.log10(R_SUN_AU * dist) + PHASE_COEFF * PHASE_DEG
    at_body = EPOCH_DAY - dist * LIGHT_S_PER_AU / 86400.0
    cycles = (at_body - ROT_MAXIMUM_DAY) * 24.0 / (ROT_PERIOD_H / 2)        # brightness cycles since the maximum
    body = H_BODY + geo - ROT_AMPLITUDE / 2 * np.cos(2 * np.pi * cycles)
    ring = H_RING + geo - 2.5 * np.log10(np.sin(np.radians(b_deg)) / np.sin(np.radians(RING_ANGLE_DEG)))
    companion = 10 ** (-0.4 * DELTA_I)
    primary = I_PAIR + 2.5 * np.log10(1 + companion)
    return 10 ** (-0.4 * (body - primary)) + 10 ** (-0.4 * (ring - primary)) + companion


def load(d, dist_km):
    d = Path(d)
    band = np.loadtxt(d / "bandpass.csv", delimiter=",", skiprows=1)
    lam, wt = band[:, 0] * 1e-9, band[:, 1] / band[:, 1].sum()
    k = dist_km * MAS
    out = []
    for line in (d / "stations.csv").read_text().splitlines()[1:]:
        c = line.split(",")
        st = dict(name=c[0], cycle=float(c[1]), exposure=float(c[2]), xi0=float(c[3]) * k, eta0=float(c[4]) * k,
                  vxi=float(c[5]) * k, veta=float(c[6]) * k)
        a = np.loadtxt(d / f"lightcurve_{c[0]}.csv", delimiter=",", skiprows=1)
        st["t"] = a[:, 0]
        n = int(round(6.0 / st["cycle"])) | 1
        st["f"] = a[:, 1] / median_filter(a[:, 1], size=n, mode="nearest")
        st["sigma"] = 1.4826 * np.median(np.abs(st["f"] - np.median(st["f"])))
        out.append(st)
    return out, lam, wt


def deep_dips(st):
    """times of the two deepest dips (the star behind the ring)"""
    k = max(int(round(0.3 / st["cycle"])), 3)
    s = np.convolve(1 - st["f"], np.ones(k) / k, "same")
    i = int(np.argmax(s))
    far = np.abs(st["t"] - st["t"][i]) > 3.0
    j = int(np.flatnonzero(far)[np.argmax(s[far])])
    out = []
    for c in (i, j):                                        # middle of the stretch below half depth
        low = np.abs(st["t"] - st["t"][c]) < 1.5
        low &= s > 0.5 * s[c]
        out.append(st["t"][low].mean() + st["exposure"] / 2)
    return sorted(out)


def conic_start(stations):
    pts = np.array([(st["xi0"] + st["vxi"] * t, st["eta0"] + st["veta"] * t) for st in stations for t in deep_dips(st)])
    x, y = pts[:, 0], pts[:, 1]
    A, Bc, C, D, E = np.linalg.lstsq(np.c_[x * x, x * y, y * y, x, y], np.ones(x.size), rcond=None)[0]
    M = np.array([[A, Bc / 2], [Bc / 2, C]])
    c = np.linalg.solve(2 * M, [-D, -E])
    val, vec = np.linalg.eigh(M / (1 + c @ M @ c))
    major, minor = 1 / np.sqrt(val[0]), 1 / np.sqrt(val[1])
    pa = np.degrees(np.arctan2(vec[0, 0], vec[1, 0])) % 180
    return [c[0], c[1], major, np.degrees(np.arcsin(minor / major)), pa]


def make_ring(p):
    mean, h, thp = np.exp(p[8]), p[9], p[7]
    return dict(centre=(p[0], p[1]), a=p[2], sin_b=np.sin(np.radians(p[3])), pa=np.radians(p[4]), w0=p[5], q=p[6],
                theta_p=thp, tau_los=lambda th: mean * (1 + h * np.cos(th - thp)))


def crossings(st, ring):
    """times at which the star crosses the mid-line of the ring"""
    t = np.arange(st["t"][0], st["t"][-1], 0.01)
    rho, _ = ring_plane(st["xi0"] + st["vxi"] * t, st["eta0"] + st["veta"] * t, ring["centre"], ring["sin_b"], ring["pa"])
    return t[np.flatnonzero(np.diff(np.sign(rho - ring["a"])))]


def windows(stations, centres, half):
    """per station: indices of the samples inside the windows and the window number of each"""
    out = []
    for st, cs in zip(stations, centres):
        lab = np.full(st["t"].size, -1)
        for n, c in enumerate(sorted(cs)):
            lab[np.abs(st["t"] + st["exposure"] / 2 - c) < half] = n
        idx = np.flatnonzero(lab >= 0)
        out.append((idx, lab[idx]))
    return out


def diffracted(stations, win, ring, table, share=None, orders=ORDERS, step=0.1):
    """the light diffracted by the particles at the samples of the windows (mid-exposure), in units of the primary"""
    out = []
    for st, (idx, lab) in zip(stations, win):
        h = np.zeros(idx.size)
        for n in np.unique(lab):
            s = lab == n
            t = st["t"][idx][s] + st["exposure"] / 2
            g = np.arange(t[0] - 2 * step, t[-1] + 3 * step, step)
            hg = halo(st["xi0"] + st["vxi"] * g, st["eta0"] + st["veta"] * g, ring, table, share=share, orders=orders)
            h[s] = CubicSpline(g, hg)(t)
        out.append(h)
    return out


def edge_rows(stations, win, ring):
    return [near_edges(st["t"][idx], st, ring, NEAR_KM) for st, (idx, _) in zip(stations, win)]


def model(stations, win, ring, lam, wt, dist_km, const, halos=None, scale=1.0, quick=False, rows=None):
    """everything in the aperture over its value outside the event.
    const: light the ring never covers, in units of the primary. quick: one wavelength and a 16-point
    exposure average, for the first steps only. rows: see ring_model.star_flux."""
    out = []
    for i, (st, (idx, _)) in enumerate(zip(stations, win)):
        r = None if rows is None else rows[i]
        if quick:
            m = star_flux(st["t"][idx], st, ring, np.array([lam @ wt]), np.ones(1), dist_km, nodes=QUICK, rows=r)
        else:
            m = star_flux(st["t"][idx], st, ring, lam, wt, dist_km, rows=r)
        if halos is not None:
            m = m + scale * halos[i]
        out.append((m + const) / (1 + const))
    return out


def residuals(stations, win, models):
    """weighted residuals; the baseline polynomial of each window multiplies the model and is solved exactly"""
    out = []
    for st, (idx, lab), m in zip(stations, win, models):
        t, f = st["t"][idx], st["f"][idx]
        for n in np.unique(lab):
            s = lab == n
            tt = (t[s] - t[s].mean()) / HALF
            G = m[s][:, None] * np.vander(tt, DEGREE + 1)
            out.append((f[s] - G @ np.linalg.lstsq(G, f[s], rcond=None)[0]) / st["sigma"])
    return np.concatenate(out)


def fit_ring(stations, win, lam, wt, dist_km, const, p0, free, halos=None, extra=None, quick=False, tol=1e-8):
    """least squares over the entries of p0 listed in `free`. extra (calibrate.py only): 'scale' adds a factor
    on the diffracted light, 'const' fits the uncovered light instead of taking it from the photometry."""
    p0 = np.array(p0, float)
    x0, xs, xl, xh = list(p0[free]), list(SCALE[free]), list(LOWER[free]), list(UPPER[free])
    if extra == "scale":
        x0, xs, xl, xh = x0 + [1.0], xs + [0.1], xl + [0.0], xh + [5.0]
    elif extra == "const":
        x0, xs, xl, xh = x0 + [const], xs + [0.02], xl + [0.0], xh + [3.0]
    rows = edge_rows(stations, win, make_ring(p0))            # fixed during the fit

    def unpack(x):
        p = p0.copy()
        p[free] = x[:len(free)]
        return p, (x[-1] if extra == "scale" else 1.0), (x[-1] if extra == "const" else const)

    def fun(x):
        p, sc, co = unpack(x)
        return residuals(stations, win, model(stations, win, make_ring(p), lam, wt, dist_km, co, halos, sc, quick, rows))

    r = least_squares(fun, x0, x_scale=xs, bounds=(xl, xh), max_nfev=60, xtol=tol, ftol=tol, gtol=tol)
    p, sc, co = unpack(r.x)
    chi2 = 2 * r.cost / (r.fun.size - r.x.size - (DEGREE + 1) * sum(len(np.unique(w[1])) for w in win))
    return p, r, chi2, sc, co


def halo_setup(kind, lam, wt, dist_km):
    """pattern table, share of the light per order, number of orders, for the reference ('full') and for the
    cruder treatments that calibrate.py compares with it"""
    def table(r_min=R_MIN_KM, r_max=R_MAX_KM, slope=SLOPE):
        return patterns(r_min, r_max, slope, lam, wt, dist_km)
    if kind == "full":
        return table(), None, ORDERS
    if kind == "single":                                    # one diffraction only
        return table(), None, 1
    one = lambda rho_pk: (rho_pk[0], rho_pk[1][:1])
    if kind == "lumped":                                    # everything diffracted and not stopped, in the single pattern
        return one(table()), (lambda x: (np.exp(-x) - np.exp(-2 * x))[None, :]), None
    raise ValueError(kind)


def solve(d, halo_kind="full", extra=None, const=None, dist_km=None, rounds=8):
    """const, dist_km, halo_kind and extra are varied by calibrate.py; the reference uses the defaults"""
    dist_km = distance_km() if dist_km is None else dist_km
    stations, lam, wt = load(d, dist_km)
    every = list(range(10))
    p = np.array(conic_start(stations) + [12.0, 0.0, 0.0, 0.5, 0.0])
    fixed = const
    const = unocculted(p[3]) if fixed is None else fixed      # with the opening angle of the conic, for the first steps
    fit = lambda win, *a, **k: fit_ring(stations, win, lam, wt, dist_km, const, *a, **k)
    win = windows(stations, [deep_dips(st) for st in stations], HALF)
    p = fit(win, p, [0, 1, 2, 3, 4, 5, 8], quick=True, tol=1e-5)[0]                 # no modulation first
    win = windows(stations, [crossings(st, make_ring(p)) for st in stations], HALF)
    best = None
    for th in np.radians([0, 90, 180, 270]):                                        # then width and depth modulation
        q0 = p.copy()
        q0[6], q0[7], q0[9] = 0.2, th, 0.2
        cand = fit(win, q0, every, quick=True, tol=1e-5)
        if best is None or cand[1].cost < best[1].cost:
            best = cand
    if fixed is None:
        const = unocculted(best[0][3])                        # with the fitted opening angle (it moves by 0.01 deg at most later)
    p, _, chi2, scale, used = fit(win, best[0], every)                              # direct beam alone, full model
    plain = chi2
    peak = 0.0
    if halo_kind is not None:
        table, share, orders = halo_setup(halo_kind, lam, wt, dist_km)
        for n in range(rounds):
            halos = diffracted(stations, win, make_ring(p), table, share, orders)
            old = p.copy()
            p = fit(win, p, [5, 6, 7, 8, 9], halos, tol=1e-6)[0]
            if abs(p[8] - old[8]) < 3e-4 and abs(p[9] - old[9]) < 3e-4 and abs(p[5] - old[5]) < 3e-3:
                break
        halos = diffracted(stations, win, make_ring(p), table, share, orders)
        p, _, chi2, scale, used = fit(win, p, every, halos, extra)
        peak = scale * max(h.max() for h in halos)
    if p[6] < 0:                                         # periapse is where the ring is narrowest
        p[6], p[7], p[9] = -p[6], p[7] + np.pi, -p[9]
    tau = np.exp(p[8]) * np.sin(np.radians(p[3])) / 2    # mean geometric normal optical depth
    values = [p[2], p[3], p[4] % 180, p[5] * (1 - p[6]), p[5] * (1 + p[6]), tau * (1 + p[9]), tau * (1 - p[9])]
    return dict(values=dict(zip(KEYS, map(float, values))), p=p, chi2=float(chi2), chi2_direct_beam_only=float(plain),
                centre=(float(p[0]), float(p[1])), tau_mean=float(tau), halo_peak=float(peak), halo_scale=float(scale),
                const=float(used), dist_km=float(dist_km))


def main(data_dir, out_path):
    r = solve(data_dir)
    out = {k: round(v, 5) for k, v in r["values"].items()}
    out.update(mean_tau=round(r["tau_mean"], 5), ring_centre_east_north_km=[round(x, 3) for x in r["centre"]],
               periapse_longitude_deg=round(float(np.degrees(r["p"][7]) % 360), 3),
               distance_au=round(r["dist_km"] / AU_KM, 6), uncovered_light_over_primary=round(r["const"], 5),
               largest_diffracted_light_over_primary=round(r["halo_peak"], 4), reduced_chi2=round(r["chi2"], 4),
               reduced_chi2_direct_beam_only=round(r["chi2_direct_beam_only"], 4))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
